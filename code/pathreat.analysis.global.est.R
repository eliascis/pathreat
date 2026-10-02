##########################################
### pathreat.analysis.global.est.R #######
### Global effects estimation ############
### Memory-efficient version #############
##########################################

library(dplyr)
library(fixest)
library(texreg)

source("code/pathreat.analysis.config.R")

######################
### load data ########
######################
{
cat("Loading matched data...\n")
if (!"d.mbase" %in% ls()) {
  d.mbase <- load_matched_data()
}
cat(sprintf("  Loaded %d observations\n", nrow(d.mbase)))
pair_index <- build_matched_pair_index(d.mbase)
cat(sprintf("  Validated %s complete matched pairs\n",
            format(pair_index$n_pairs, big.mark = ",")))
}

############################
### helper functions #######
############################

# extract coefficients and stats from estimation list (compact storage)
extract_coef <- function(est_list, estimate_type, ctrl_means = NULL) {
  do.call(rbind, lapply(names(est_list), function(v) {
    e <- est_list[[v]]
    ct <- e$coeftable["treat", ]
      estimate <- regression_estimate(e)
    data.frame(
      variable = v,
      estimate_type = estimate_type,
      coef = ct["Estimate"],
      se = ct["Std. Error"],
      t_stat = ct["t value"],
      pval = ct["Pr(>|t|)"],
      ci_low = estimate$ci_low,
      ci_high = estimate$ci_high,
      n_obs = e$nobs,
      control_mean = if (!is.null(ctrl_means)) ctrl_means[v] else NA,
      row.names = NULL
    )
  }))
}

# run single regression on the canonical outcome-specific pair sample
run_reg <- function(dep, data, prefix = "") {
  base_dep <- sub("^(sd\\.|nr\\.|pa\\.)", "", dep)
  sample_mask <- matched_pair_sample_mask(
    data,
    outcome = dep,
    pair_index = pair_index,
    filter_outcome = base_dep
  )
  n_pairs <- attr(sample_mask, "n_pairs")

  fe_spec <- "country_rast + biome_raster"
  clustvar <- c("country_rast")
  f <- as.formula(paste(dep, "~ treat |", fe_spec))

  if (anyNA(data$country_rast[sample_mask]) ||
      anyNA(data$biome_raster[sample_mask])) {
    stop(sprintf("%s has missing fixed-effect labels in its eligible pairs", dep))
  }
  n_treated <- sum(sample_mask & data$treat == 1)
  n_control <- sum(sample_mask & data$treat == 0)
  if (n_treated != n_pairs || n_control != n_pairs) {
    stop(sprintf("%s canonical sample is not pair balanced", dep))
  }

  e <- feols(f, data = data, cluster = clustvar, subset = sample_mask)
  if (e$nobs != 2L * n_pairs) {
    stop(sprintf("%s regression dropped rows from the canonical pair sample", dep))
  }

  list(
    fit = e,
    control_mean = mean(data[[dep]][sample_mask & data$treat == 0]),
    n_pairs = n_pairs,
    n_countries = length(unique(data$country_rast[sample_mask]))
  )
}

############################
### 1. UNSCALED estimates ##
############################
{
cat("\n=== Block 1: Unscaled estimates ===\n")

# work on a copy
d <- d.mbase

unscaled_runs <- lapply(deplist, function(dep) run_reg(dep, d))
names(unscaled_runs) <- deplist
e.unscaled <- lapply(unscaled_runs, `[[`, "fit")
names(e.unscaled) <- deplist

# control means and country counts use the exact regression samples
ctrl_means <- vapply(unscaled_runs, `[[`, numeric(1), "control_mean")

h <- htmlreg(
  e.unscaled,
  omit.coef = "(Intercept)",
  custom.model.names = names(e.unscaled),
  custom.gof.rows = list("Control Mean" = ctrl_means),
  stars = c(0.01, 0.05, 0.1),
  digits = 3,
  table = F
)

cat(file = paste0(paths$results_dir, "pathreat.est.newest.html"),
    c("<b>Unscaled estimates (original coding, any_threat = presence/absence)</b>", h), append = F)
cat(file = paste0(paths$results_dir, "pathreat.global.est.html"),
    c("<b>Unscaled estimates (original coding, any_threat = presence/absence)</b>", h), append = F)

# extract results
est_results <- extract_coef(e.unscaled, "unscaled", ctrl_means)
}

############################
### LaTeX Table A.2 ########
############################
{
cat("  Generating LaTeX Table A.2...\n")

# helper functions
fmt_num <- function(x) {
  if (is.na(x)) return("")
  formatC(x, format = "f", digits = 3)
}

sig_stars <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.01) return("***")
  if (p < 0.05) return("**")
  if (p < 0.1)  return("*")
  return("")
}

fmt_coef <- function(coef, pval) {
  if (is.na(coef)) return("")
  val <- fmt_num(coef)
  stars <- sig_stars(pval)
  if (nchar(stars) > 0) paste0(val, "\\sym{", stars, "}") else val
}

fmt_se <- function(se) {
  if (is.na(se)) return("")
  paste0("(", fmt_num(se), ")")
}

fmt_int <- function(x) formatC(x, format = "d", big.mark = ",")

# panel layout: 4 panels (5+5+4+2 vars)
panels <- list(
  A = c("built", "cropland", "planted", "pasture", "oil"),
  B = c("mining", "renewables", "roads", "powerlines", "fires"),
  C = c("swu", "dams", "light", "def0120_parea"),
  D = c("any_threat", "threat_composite")
)

# convert deplist_scale_labels to LaTeX column headers (name row + unit row)
# long names get \makecell{Line1\\Line2} for multi-line headers
header_breaks <- c(
  "Forest plantation"          = "Forest\\\\plantation",
  "No. of cattle"              = "No.\\ of\\\\cattle",
  "Oil \\& gas operation"      = "Oil \\& gas\\\\operation",
  "Mining operations"          = "Mining\\\\operations",
  "Renewable energy presence"  = "Renewable energy\\\\presence",
  "Road presence"              = "Road\\\\presence",
  "Power lines presence"       = "Power lines\\\\presence",
  "Burned area presence"       = "Burned area\\\\presence",
  "Surface water change"       = "Surface water\\\\change",
  "Dams presence"              = "Dams\\\\presence",
  "Forest loss 2001--20"       = "Forest loss\\\\2001--20",
  "Threat composite"           = "Threat\\\\composite"
)

make_header <- function(label) {
  # convert units to LaTeX
  label <- unit_label_tex(label)
  label <- gsub("%", "\\\\%", label)
  label <- gsub("&", "\\\\&", label)
  label <- gsub("2001-2020", "2001--20", label)
  # split on first [ to get name and unit
  if (grepl("\\[", label)) {
    parts <- regmatches(label, regexec("^(.+?)\\s*(\\[.+\\])$", label))[[1]]
    name_part <- trimws(parts[2])
    unit_part <- parts[3]
  } else {
    name_part <- label
    unit_part <- ""
  }
  # apply multi-line breaks for long names (normalize whitespace for matching)
  name_part <- gsub("\\s+", " ", trimws(name_part))
  if (name_part %in% names(header_breaks)) {
    name_part <- paste0("\\makecell{", header_breaks[name_part], "}")
  }
  list(name = name_part, unit = unit_part)
}

# collect stats per variable
var_stats <- lapply(deplist, function(v) {
  e <- e.unscaled[[v]]
  ct <- coeftable(e)["treat", ]
  n_countries <- unscaled_runs[[v]]$n_countries
  list(
    coef = ct["Estimate"],
    se = ct["Std. Error"],
    pval = ct["Pr(>|t|)"],
    ctrl_mean = ctrl_means[v],
    n_obs = e$nobs,
    n_countries = n_countries,
    ar2 = fitstat(e, "ar2")$ar2
  )
})
names(var_stats) <- deplist

# build LaTeX lines
lines <- c()
table_scale_labels <- deplist_scale_labels
table_scale_labels["light"] <- "Light intensity [DN]"

for (p in seq_along(panels)) {
  panel_name <- names(panels)[p]
  panel_vars <- panels[[panel_name]]
  n_vars <- length(panel_vars)
  n_empty <- 5 - n_vars  # pad if fewer than 5

  # panel header
  lines <- c(lines, paste0("% Panel ", panel_name))
  if (p == 1) {
    lines <- c(lines, paste0("%\\textit{Panel ", panel_name, ":} & & & & & \\\\"))
    lines <- c(lines, "%\\addlinespace[.5ex]")
  } else {
    lines <- c(lines, paste0("%\\multicolumn{6}{l}{\\textit{Panel ", panel_name, ":}} \\\\"))
    lines <- c(lines, "%\\addlinespace[.5ex]")
    lines <- c(lines, "\\midrule")
    lines <- c(lines, "\\addlinespace[2.5ex]")
  }

  # column headers (name row + unit row)
  headers <- lapply(panel_vars, function(v) make_header(table_scale_labels[v]))
  name_cells <- sapply(headers, function(h) paste0("\\mco{", h$name, "}"))
  unit_cells <- sapply(headers, function(h) if (nchar(h$unit) > 0) paste0("\\mco{", h$unit, "}") else "")
  # pad with empty cells
  if (n_empty > 0) {
    name_cells <- c(name_cells, rep("", n_empty))
    unit_cells <- c(unit_cells, rep("", n_empty))
  }
  lines <- c(lines, paste0("& ", paste(name_cells, collapse = " & "), " \\\\"))
  lines <- c(lines, paste0("& ", paste(unit_cells, collapse = " & "), " \\\\"))

  # column numbers (continuous across panels, below names/units)
  col_offset <- (p - 1) * 5
  num_cells <- sapply(seq_len(n_vars), function(i) paste0("\\mco{(", col_offset + i, ")}"))
  if (n_empty > 0) num_cells <- c(num_cells, rep("", n_empty))
  lines <- c(lines, paste0("& ", paste(num_cells, collapse = " & "), " \\\\"))
  lines <- c(lines, "\\midrule")

  # coefficient row
  coef_cells <- sapply(panel_vars, function(v) fmt_coef(var_stats[[v]]$coef, var_stats[[v]]$pval))
  if (n_empty > 0) coef_cells <- c(coef_cells, rep("", n_empty))
  lines <- c(lines, paste0("Protected area & ", paste(coef_cells, collapse = " & "), " \\\\"))

  # SE row
  se_cells <- sapply(panel_vars, function(v) fmt_se(var_stats[[v]]$se))
  if (n_empty > 0) se_cells <- c(se_cells, rep("", n_empty))
  lines <- c(lines, paste0("& ", paste(se_cells, collapse = " & "), " \\\\"))

  lines <- c(lines, "\\midrule")

  # control mean
  cm_cells <- sapply(panel_vars, function(v) fmt_num(var_stats[[v]]$ctrl_mean))
  if (n_empty > 0) cm_cells <- c(cm_cells, rep("", n_empty))
  lines <- c(lines, paste0("Control mean & ", paste(cm_cells, collapse = " & "), " \\\\"))

  # N observations (in braces for S-column protection)
  nobs_cells <- sapply(panel_vars, function(v) paste0("{", fmt_int(var_stats[[v]]$n_obs), "}"))
  if (n_empty > 0) nobs_cells <- c(nobs_cells, rep("", n_empty))
  lines <- c(lines, paste0("Num.\\ obs. & ", paste(nobs_cells, collapse = " & "), " \\\\"))

  # N countries
  ncountry_cells <- sapply(panel_vars, function(v) var_stats[[v]]$n_countries)
  if (n_empty > 0) ncountry_cells <- c(ncountry_cells, rep("", n_empty))
  lines <- c(lines, paste0("Num.\\ countries & ", paste(ncountry_cells, collapse = " & "), " \\\\"))

  # Adj. R²
  r2_cells <- sapply(panel_vars, function(v) paste0("\\mco{", formatC(var_stats[[v]]$ar2, format = "f", digits = 3), "}"))
  if (n_empty > 0) r2_cells <- c(r2_cells, rep("", n_empty))
  r2_line <- paste0("{Adj.\\ $R^2$} & ", paste(r2_cells, collapse = " & "), " \\\\")

  if (p < length(panels)) {
    lines <- c(lines, r2_line)
  } else {
    # last panel: strip trailing \\
    lines <- c(lines, sub(" *\\\\\\\\$", "", r2_line))
  }
}

# write table
tex_file <- "pub/tables/tab.ATT_by_threat_table.tex"
writeLines(lines, tex_file)
cat("  Written:", tex_file, "\n")


rm(d, e.unscaled, h, ctrl_means, unscaled_runs); gc()
}

###########################################
### 2. GROUP MEANS with clustered SEs #####
###########################################
{
cat("\n=== Block 2: Group means with clustered SEs (unscaled) ===\n")

# work on a copy
d <- d.mbase

# function to compute group mean with clustered SE
get_group_stats <- function(dep, data, treat_val, sample_mask,
                            cluster_var = "country_rast") {
  # run intercept-only regression with clustering
  f <- as.formula(paste(dep, "~ 1"))
  group_mask <- sample_mask & data$treat == treat_val
  e <- feols(f, data = data, cluster = cluster_var, subset = group_mask)

  data.frame(
    mean = coef(e)["(Intercept)"],
    se = e$se["(Intercept)"],
    n = e$nobs
  )
}

# compute for all variables
group_means_list <- lapply(deplist, function(v) {
  cat("  Processing:", v, "\n")
  sample_mask <- matched_pair_sample_mask(d, v, pair_index)
  n_pairs <- attr(sample_mask, "n_pairs")

  ctrl <- get_group_stats(v, d, treat_val = 0, sample_mask = sample_mask)
  treat <- get_group_stats(v, d, treat_val = 1, sample_mask = sample_mask)
  if (ctrl$n != n_pairs || treat$n != n_pairs) {
    stop(sprintf("%s group means do not use one row per pair member", v))
  }

  data.frame(
    variable = v,
    control_mean = ctrl$mean,
    control_se = ctrl$se,
    control_n = ctrl$n,
    treated_mean = treat$mean,
    treated_se = treat$se,
    treated_n = treat$n,
    row.names = NULL
  )
})

group_means <- do.call(rbind, group_means_list)
rownames(group_means) <- NULL

rm(d); gc()
}

###########################################################
### 3. FE-ADJUSTED LEVELS with country-pairs bootstrap ####
###########################################################
{
cat("\n=== Block 3: FE-adjusted levels with country-pairs bootstrap ===\n")

# fixed production settings; the named order preserves the validated prototype seeds
fe_level_n_boot <- 499L
fe_level_seed_base <- 20260712L
fe_level_outcome_order <- c(
  "built",
  "cropland",
  "planted",
  "pasture",
  "def0120_parea",
  "oil",
  "mining",
  "renewables",
  "roads",
  "powerlines",
  "fires",
  "swu",
  "dams",
  "light",
  "threat_composite",
  "any_threat"
)
fe_level_seed_map <- setNames(
  fe_level_seed_base + seq_along(fe_level_outcome_order),
  fe_level_outcome_order
)
fe_level_max_attempts <- fe_level_n_boot +
  max(50L, ceiling(0.20 * fe_level_n_boot))
fe_level_att_tolerance <- 1e-6
fe_level_identity_tolerance <- 1e-8

if (!setequal(fe_level_outcome_order, deplist)) {
  stop("FE-adjusted level outcome order does not match deplist")
}

est_unscaled_reference <- est_results[
  est_results$estimate_type == "unscaled" &
    est_results$variable %in% fe_level_outcome_order,
]

if (!setequal(est_unscaled_reference$variable, fe_level_outcome_order)) {
  stop("Production estimates do not contain every FE-adjusted level outcome")
}

fe_level_weighted_mean <- function(x, w) {
  sum(x * w) / sum(w)
}

fe_level_collapse_outcome <- function(variable, data) {
  keep <- matched_pair_sample_mask(data, variable, pair_index)
  if (anyNA(data$country_rast[keep]) || anyNA(data$biome_raster[keep])) {
    stop(sprintf("%s has missing FE labels in its canonical sample", variable))
  }

  d_variable <- data[keep, c(
    "country_rast",
    "biome_raster",
    "treat",
    variable
  )]
  names(d_variable)[names(d_variable) == variable] <- "outcome"

  cells <- d_variable %>%
    group_by(country_rast, biome_raster, treat) %>%
    summarise(
      cell_n = n(),
      cell_mean = mean(outcome),
      .groups = "drop"
    ) %>%
    data.frame()

  cells$variable <- variable
  cells
}

fe_level_fit_cells <- function(cells, weight_name) {
  weight_formula <- as.formula(paste0("~", weight_name))
  feols(
    cell_mean ~ treat | country_rast + biome_raster,
    data = cells,
    weights = weight_formula,
    fixef.rm = "none",
    warn = FALSE,
    notes = FALSE
  )
}

fe_level_get_margins <- function(fit, cells, weight_name) {
  target <- cells[cells$treat == 1, ]
  if (nrow(target) == 0) {
    stop("No protected observations in the FE-level standardization target")
  }

  target_control <- target
  target_protected <- target
  target_control$treat <- 0
  target_protected$treat <- 1

  pred_control <- predict(fit, newdata = target_control)
  pred_protected <- predict(fit, newdata = target_protected)
  target_weights <- target[[weight_name]]

  if (any(!is.finite(pred_control)) || any(!is.finite(pred_protected))) {
    stop("Non-finite FE-level prediction")
  }

  adjusted_control <- fe_level_weighted_mean(pred_control, target_weights)
  adjusted_protected <- fe_level_weighted_mean(pred_protected, target_weights)
  att <- unname(coef(fit)["treat"])

  if (!is.finite(att)) {
    stop("FE-level treatment coefficient is not finite")
  }

  c(
    adjusted_control = adjusted_control,
    adjusted_protected = adjusted_protected,
    adjusted_gap = adjusted_protected - adjusted_control,
    att = att,
    identity_error = abs((adjusted_protected - adjusted_control) - att)
  )
}

fe_level_get_raw_stats <- function(cells) {
  control <- cells[cells$treat == 0, ]
  protected <- cells[cells$treat == 1, ]

  raw_control <- fe_level_weighted_mean(control$cell_mean, control$cell_n)
  raw_protected <- fe_level_weighted_mean(
    protected$cell_mean,
    protected$cell_n
  )

  c(
    raw_control = raw_control,
    raw_protected = raw_protected,
    raw_gap = raw_protected - raw_control
  )
}

fe_level_run_boot_draw <- function(cells, countries) {
  draw <- sample(countries, length(countries), replace = TRUE)
  multiplicity <- tabulate(
    match(draw, countries),
    nbins = length(countries)
  )

  boot_cells <- cells
  boot_cells$boot_mult <- multiplicity[
    match(boot_cells$country_rast, countries)
  ]
  boot_cells <- boot_cells[boot_cells$boot_mult > 0, ]
  boot_cells$boot_n <- boot_cells$cell_n * boot_cells$boot_mult

  fit <- suppressWarnings(fe_level_fit_cells(boot_cells, "boot_n"))
  margins <- fe_level_get_margins(fit, boot_cells, "boot_n")

  if (margins["identity_error"] > fe_level_identity_tolerance) {
    stop(sprintf(
      "Bootstrap FE fit failed the margin identity check: %.12g",
      margins["identity_error"]
    ))
  }

  data.frame(
    adjusted_control = margins["adjusted_control"],
    adjusted_protected = margins["adjusted_protected"],
    adjusted_gap = margins["adjusted_gap"],
    att = margins["att"],
    identity_error = margins["identity_error"],
    error = NA_character_,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

fe_level_safe_boot_draw <- function(cells, countries) {
  tryCatch(
    fe_level_run_boot_draw(cells, countries),
    error = function(e) {
      data.frame(
        adjusted_control = NA_real_,
        adjusted_protected = NA_real_,
        adjusted_gap = NA_real_,
        att = NA_real_,
        identity_error = NA_real_,
        error = conditionMessage(e),
        stringsAsFactors = FALSE
      )
    }
  )
}

fe_level_run_outcome <- function(variable, cells) {
  set.seed(fe_level_seed_map[[variable]])

  production <- est_unscaled_reference[
    est_unscaled_reference$variable == variable,
  ]
  production_means <- group_means[group_means$variable == variable, ]

  full_fit <- fe_level_fit_cells(cells, "cell_n")
  full_margins <- fe_level_get_margins(full_fit, cells, "cell_n")
  raw_stats <- fe_level_get_raw_stats(cells)

  production_att_diff <- abs(full_margins["att"] - production$coef)
  n_obs <- sum(cells$cell_n)
  n_obs_diff <- n_obs - production$n_obs
  raw_control_diff <- abs(
    raw_stats["raw_control"] - production_means$control_mean
  )
  raw_protected_diff <- abs(
    raw_stats["raw_protected"] - production_means$treated_mean
  )

  if (production_att_diff > fe_level_att_tolerance) {
    stop(sprintf(
      "%s: collapsed ATT differs from production by %.12g",
      variable,
      production_att_diff
    ))
  }
  if (full_margins["identity_error"] > fe_level_identity_tolerance) {
    stop(sprintf(
      "%s: full-sample margin identity error is %.12g",
      variable,
      full_margins["identity_error"]
    ))
  }
  if (n_obs_diff != 0) {
    stop(sprintf(
      "%s: collapsed sample has %s observations but production has %s",
      variable,
      format(n_obs, scientific = FALSE),
      format(production$n_obs, scientific = FALSE)
    ))
  }
  if (max(raw_control_diff, raw_protected_diff) > fe_level_att_tolerance) {
    stop(sprintf(
      "%s: collapsed raw means do not reproduce production group means",
      variable
    ))
  }

  countries <- sort(unique(cells$country_rast))
  draws <- lapply(seq_len(fe_level_n_boot), function(i) {
    fe_level_safe_boot_draw(cells, countries)
  })
  draws <- do.call(rbind, draws)
  attempts <- fe_level_n_boot

  is_valid <- complete.cases(draws[, c(
    "adjusted_control",
    "adjusted_protected",
    "adjusted_gap",
    "att",
    "identity_error"
  )])

  while (sum(is_valid) < fe_level_n_boot &&
         attempts < fe_level_max_attempts) {
    replacements_needed <- min(
      fe_level_n_boot - sum(is_valid),
      fe_level_max_attempts - attempts
    )
    replacement_draws <- lapply(seq_len(replacements_needed), function(i) {
      fe_level_safe_boot_draw(cells, countries)
    })
    replacement_draws <- do.call(rbind, replacement_draws)
    draws <- rbind(draws, replacement_draws)
    attempts <- attempts + replacements_needed
    is_valid <- complete.cases(draws[, c(
      "adjusted_control",
      "adjusted_protected",
      "adjusted_gap",
      "att",
      "identity_error"
    )])
  }

  draws <- draws[is_valid, ]
  if (nrow(draws) < fe_level_n_boot) {
    stop(sprintf(
      "%s: obtained only %d valid draws after %d attempts",
      variable,
      nrow(draws),
      attempts
    ))
  }
  draws <- draws[seq_len(fe_level_n_boot), ]
  draws$variable <- variable
  draws$replicate <- seq_len(fe_level_n_boot)

  if (max(draws$identity_error) > fe_level_identity_tolerance) {
    stop(sprintf(
      "%s: maximum bootstrap margin identity error is %.12g",
      variable,
      max(draws$identity_error)
    ))
  }

  quantile_values <- function(x) {
    unname(quantile(x, probs = c(0.025, 0.975), type = 6))
  }

  control_ci <- quantile_values(draws$adjusted_control)
  protected_ci <- quantile_values(draws$adjusted_protected)
  gap_ci <- quantile_values(draws$adjusted_gap)

  summary_row <- data.frame(
    variable = variable,
    raw_control_mean = raw_stats["raw_control"],
    raw_protected_mean = raw_stats["raw_protected"],
    raw_gap = raw_stats["raw_gap"],
    raw_control_se = production_means$control_se,
    raw_protected_se = production_means$treated_se,
    adjusted_control_mean = full_margins["adjusted_control"],
    adjusted_protected_mean = full_margins["adjusted_protected"],
    adjusted_gap = full_margins["adjusted_gap"],
    bootstrap_control_se = sd(draws$adjusted_control),
    bootstrap_control_ci_low = control_ci[1],
    bootstrap_control_ci_high = control_ci[2],
    bootstrap_protected_se = sd(draws$adjusted_protected),
    bootstrap_protected_ci_low = protected_ci[1],
    bootstrap_protected_ci_high = protected_ci[2],
    bootstrap_gap_se = sd(draws$adjusted_gap),
    bootstrap_gap_ci_low = gap_ci[1],
    bootstrap_gap_ci_high = gap_ci[2],
    production_att = production$coef,
    production_se = production$se,
    production_pval = production$pval,
    collapsed_att = full_margins["att"],
    production_att_abs_diff = production_att_diff,
    full_identity_error = full_margins["identity_error"],
    max_boot_identity_error = max(draws$identity_error),
    raw_control_abs_diff = raw_control_diff,
    raw_protected_abs_diff = raw_protected_diff,
    n_obs = n_obs,
    n_control = sum(cells$cell_n[cells$treat == 0]),
    n_protected = sum(cells$cell_n[cells$treat == 1]),
    n_countries = length(countries),
    n_biomes = length(unique(cells$biome_raster)),
    n_cells = nrow(cells),
    bootstrap_reps_requested = fe_level_n_boot,
    bootstrap_reps_valid = nrow(draws),
    bootstrap_attempts = attempts,
    bootstrap_invalid = attempts - nrow(draws),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  cat(sprintf(
    "  %-18s ATT=% .6g; adjusted means=(% .6g, % .6g); valid=%d\n",
    variable,
    summary_row$production_att,
    summary_row$adjusted_control_mean,
    summary_row$adjusted_protected_mean,
    summary_row$bootstrap_reps_valid
  ))

  list(summary = summary_row, draws = draws)
}

cat("  Collapsing outcomes to country-biome-treatment cells...\n")
fe_level_cells <- lapply(fe_level_outcome_order, function(variable) {
  cat("    Collapsing", variable, "\n")
  fe_level_collapse_outcome(variable, d.mbase)
})
names(fe_level_cells) <- fe_level_outcome_order

old_fixest_nthreads <- fixest::getFixest_nthreads()
fixest::setFixest_nthreads(1)

cat(sprintf(
  "  Running %d valid country-pairs bootstrap draws per outcome...\n",
  fe_level_n_boot
))

fe_level_n_workers <- min(no_cluster, length(fe_level_outcome_order))
fe_level_results <- tryCatch(
  {
    if (.Platform$OS.type != "windows" && fe_level_n_workers > 1) {
      parallel::mclapply(
        fe_level_outcome_order,
        function(variable) {
          fe_level_run_outcome(variable, fe_level_cells[[variable]])
        },
        mc.cores = fe_level_n_workers,
        mc.preschedule = FALSE,
        mc.set.seed = FALSE
      )
    } else {
      lapply(fe_level_outcome_order, function(variable) {
        fe_level_run_outcome(variable, fe_level_cells[[variable]])
      })
    }
  },
  finally = fixest::setFixest_nthreads(old_fixest_nthreads)
)
names(fe_level_results) <- fe_level_outcome_order

fe_level_failed <- vapply(
  fe_level_results,
  inherits,
  logical(1),
  what = "try-error"
)
if (any(fe_level_failed)) {
  errors <- vapply(
    fe_level_results[fe_level_failed],
    as.character,
    character(1)
  )
  stop(paste(
    "FE-adjusted level estimation failed:",
    paste(errors, collapse = "\n")
  ))
}

fe_level_summary_list <- lapply(fe_level_results, `[[`, "summary")
fe_level_summary <- do.call(rbind, fe_level_summary_list)
rownames(fe_level_summary) <- NULL
fe_level_summary <- fe_level_summary[
  match(deplist, fe_level_summary$variable),
]
rownames(fe_level_summary) <- NULL

fe_level_draws <- lapply(fe_level_results, `[[`, "draws")
names(fe_level_draws) <- fe_level_outcome_order

fe_level_settings <- list(
  estimand = paste(
    "Both counterfactual means standardized over protected pixels",
    "in each estimation sample"
  ),
  fixed_effects = c("country_rast", "biome_raster"),
  bootstrap = "Country-pairs bootstrap with replacement",
  bootstrap_reps = fe_level_n_boot,
  seed_base = fe_level_seed_base,
  seed_map = fe_level_seed_map,
  significance_source = "Production country-clustered FE p-value",
  att_tolerance = fe_level_att_tolerance,
  identity_tolerance = fe_level_identity_tolerance,
  outcome_order = deplist,
  generated_at = format(Sys.time(), tz = "Europe/Berlin", usetz = TRUE)
)

fe_level_bundle <- list(
  settings = fe_level_settings,
  summary = fe_level_summary,
  bootstrap_draws = fe_level_draws
)

fe_level_csv_file <- file.path(
  paths$results_dir,
  "pathreat.global.fe-adjusted-levels.est.csv"
)
write.csv(fe_level_summary, fe_level_csv_file, row.names = FALSE)
saveRDS(fe_level_bundle, paths$est_global_fe_levels)

cat("  Saved FE-adjusted level results to:", paths$est_global_fe_levels, "\n")
cat("  Saved FE-adjusted level summary to:", fe_level_csv_file, "\n")

rm(fe_level_cells, fe_level_results)
invisible(gc())
}

############################
### Save results ###########
############################
{
# save treatment effects to .Rds file
results_file <- paste0(paths$results_dir, "pathreat.global.est.Rds")
saveRDS(est_results, results_file)
cat("\nSaved treatment effect results to:", results_file, "\n")
cat("  Dimensions:", nrow(est_results), "rows x", ncol(est_results), "columns\n")
cat("  Estimate types:", paste(unique(est_results$estimate_type), collapse = ", "), "\n")

# save group means with clustered SEs to separate .Rds file
group_means_file <- paste0(paths$results_dir, "pathreat.global.group_means.est.Rds")
saveRDS(group_means, group_means_file)
cat("\nSaved group means to:", group_means_file, "\n")
cat("  Dimensions:", nrow(group_means), "rows x", ncol(group_means), "columns\n")
}

cat("\n=== Global effects estimation complete ===\n")
cat("Output HTML: ", paths$results_dir, "pathreat.est.newest.html\n")
cat("Output HTML: ", paths$results_dir, "pathreat.global.est.html\n")
cat("Output Rds:  ", paths$results_dir, "pathreat.global.est.Rds\n")
cat("Note: any_threat uses presence/absence coding (binary)\n")
