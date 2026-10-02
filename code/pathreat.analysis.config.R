#####################################
### pathreat.analysis.config.R ####
### Shared parameters and settings ##
#####################################

#########################
### paths ###############
#########################

paths <- list(
  # inputs
  data_unmatched        = "data/store/pathreat.data.merge.unmatched.fst",
  data_outcome_filters   = "data/store/pathreat.data.outcome-sample-filters.Rds",
  data_matched          = "data/store/pathreat.data.merge.matched.fst",
  data_summary   = "data/store/Data_summary.xlsx",
  # outputs
  results_dir    = "results/",
  figures_dir    = "pub/figures/",
  # estimation results
  est_global      = "results/pathreat.global.est.Rds",
  est_group_means = "results/pathreat.global.group_means.est.Rds",
  est_global_fe_levels = "results/pathreat.global.fe-adjusted-levels.est.Rds",
  est_bycountry   = "results/pathreat.bycountry.est.Rds",
  est_bybiome     = "results/pathreat.bybiome.est.Rds",
  est_PA_type_group_means = "results/pathreat.PA-type.group_means.est.Rds",
  est_sensitivity_rosenbaum = "results/pathreat.sensitivity.rosenbaum.Rds",
  est_sfa         = "results/pathreat.sfa.est.Rds",
  est_sfa_country = "results/pathreat.sfa.country.est.Rds",
  # cached aggregations
  country_violin_means = "data/store/pathreat.data.global.violins.country-means.pair-complete.Rds"
)

# output directories (created on demand so a fresh clone without pub/ content runs)
for (d in c(paths$results_dir, paths$figures_dir, "pub/figures/presentation", "pub/tables")) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
rm(d)

#################################
### figure source data helper ###
#################################

# Writes the data frame handed to ggplot() for a manuscript figure to
# data/store/<figure root>.csv, so every figure has a same-named Source Data
# file (e.g. pub/figures/fig.delta_map.jpg -> data/store/fig.delta_map.csv).
save_figure_data <- function(df, fig_file, dir = "data/store") {
  if (inherits(df, "sf")) {
    df <- sf::st_drop_geometry(df)
  }
  df <- as.data.frame(df)
  root <- sub("\\.[A-Za-z0-9]+$", "", basename(fig_file))
  out_file <- file.path(dir, paste0(root, ".csv"))
  write.csv(df, out_file, row.names = FALSE)
  cat(sprintf("Saved figure data: %s (%d rows, %d cols)\n", out_file, nrow(df), ncol(df)))
  invisible(out_file)
}

#################################
### load matched data helper ####
#################################

load_matched_data <- function() {
  d <- fst::read_fst(paths$data_matched)
  cat(sprintf("Loaded matched data: %s obs (%d countries)\n",
              format(nrow(d), big.mark = ","), length(unique(d$country_rast))))
  d
}

######################
### variable lists ###
######################

# matching covariates
mlist <- c(
  "elevation",
  "pop_2000",
  "annual_total_precipitation",
  "slope",
  "temperature",
  "access",
  "cropsuit2",
  # "treecover_2000"
  "forest2000_parea"
)

# matching covariate labels (from Data_summary.xlsx)
data_summary <- readxl::read_excel(paths$data_summary)
mlist_labels_df <- data_summary[data_summary$varnames %in% mlist, c("varnames", "varlabels", "order")]
mlist_labels_df <- mlist_labels_df[order(mlist_labels_df$order), ]
mlist_labels <- setNames(mlist_labels_df$varlabels, mlist_labels_df$varnames)
mlist_labels <- mlist_labels[mlist]  # reorder to match mlist

# PA-level identifier columns (for SFA and PA-level analyses)
pa_cols <- c("wdpaid", "iucn_class", "size_class", "size",
             "country_rast", "biome_raster", "x", "y")

# 14 individual threat variables (raw, no summary variables)
threat_list_raw <- c(
  "built", "cropland", "planted", "pasture",
  "oil", "mining", "renewables",
  "roads", "powerlines",
  "def0120_parea",
  "fires", "swu", "dams",
  "light"
)

# threat categories (IUCN classification)
threat_categories <- list(
  tc_residential  = c("built"),
  tc_agriculture  = c("cropland", "planted", "pasture"),
  tc_energy       = c("oil", "mining", "renewables"),
  tc_transport    = c("roads", "powerlines"),
  tc_biological   = c("def0120_parea"),
  tc_modification = c("fires", "swu", "dams"),
  tc_pollution    = c("light")
)

# composite category index names and labels
composite_catlist <- names(threat_categories)
names(composite_catlist) <- composite_catlist
composite_catlist_labels <- c(
  tc_residential  = "Residential & commercial dev.",
  tc_agriculture  = "Agriculture & aquaculture",
  tc_energy       = "Energy production & mining",
  tc_transport    = "Transportation & service corridors",
  tc_biological   = "Biological resource use",
  tc_modification = "Natural system modifications",
  tc_pollution    = "Pollution"
)

# outcome variables (threats + summary)
deplist <- c(
  "built",
  "cropland",
  "planted",
  "pasture",
  "oil",
  "mining",
  "renewables",
  "roads",
  "powerlines",
  "fires",
  "swu",
  "dams",
  "light",
  "def0120_parea",
  "any_threat",
  "threat_composite"
)
names(deplist) <- deplist

# outcome variable labels (from Data_summary.xlsx)
deplist_labels_df <- data_summary[!is.na(data_summary$threat.no), ]
deplist_labels <- setNames(deplist_labels_df$threat.label, deplist_labels_df$variable)

# outcome variable labels with units for figure/table axes (from Data_summary.xlsx)
deplist_scale_labels <- setNames(deplist_labels_df$scale.label, deplist_labels_df$variable)

# unit codes in scale.label are ASCII ([sqm], [sqkm]); convert per output medium
unit_label_tex <- function(x) {
  x <- gsub("[sqkm]", "[km\\textsuperscript{2}]", x, fixed = TRUE)
  x <- gsub("[sqm]", "[m\\textsuperscript{2}]", x, fixed = TRUE)
  x
}
unit_label_plot <- function(x) {
  x <- gsub("[sqkm]", "[km\u00B2]", x, fixed = TRUE)
  x <- gsub("[sqm]", "[m\u00B2]", x, fixed = TRUE)
  x
}

# thresholds for presence/absence coding of any_threat and for the country
# prevalence weights of the threat composite.
# Uniform top-quartile rule: a pixel counts as threatened when it sits in the
# upper quartile of the threat's global distribution - Q75 for the thirteen
# indicators where higher means worse, Q25 for surface water abstraction, whose
# direction is reversed (lower values = greater water loss). A literal Q75 for
# swu would flag two thirds of observed pixels and is not the intended analogue.
# For nine of the fourteen indicators Q75 is zero, so the rule coincides with
# the fixed "> 0" convention used before 2026-08-26; only cropland (Q75 = 1) and
# def0120_parea (Q75 = 0.0027) change coding. See ideas/20260826_q75_threshold.md.
# value: resolved threshold, filled in by the consuming script
# quantile_prob: quantile of the GLOBAL merged pixel data
#   (data/store/pathreat.data.merge.fst) used as the cutoff
# dir: ">" means value > threshold = presence; "<" means value < threshold = presence
dep_thresholds <- data.frame(
  variable      = c("built", "cropland", "planted", "pasture", "oil", "mining",
                     "renewables", "roads", "powerlines", "fires", "swu", "dams",
                     "light", "def0120_parea"),
  value         = rep(NA_real_, 14),
  quantile_prob = c(rep(0.75, 10), 0.25, rep(0.75, 3)),
  dir           = c(rep(">", 10), "<", rep(">", 3)),
  stringsAsFactors = FALSE
)
rownames(dep_thresholds) <- dep_thresholds$variable

##################################
### treated-side target filters ##
##################################
# per-outcome treated-side target expressions (built by
# pathreat.data.outcome-sample-filters.R)
# balance filter (avg |SMD| <= 0.10) is pre-applied in merge.matched.fst
outcome_target_filters_available <- file.exists(paths$data_outcome_filters)
if (outcome_target_filters_available) {
  dep_sample_filters <- readRDS(paths$data_outcome_filters)
} else {
  cat("Note: outcome filters not found (run pathreat.data.outcome-sample-filters.R first)\n")
  dep_sample_filters <- setNames(replicate(length(deplist), NULL), deplist)
}

########################################
### canonical matched-pair samples ####
########################################

base_outcome_name <- function(outcome) {
  sub("^(sd\\.|nr\\.|pa\\.)", "", outcome)
}

get_outcome_target_filter <- function(outcome) {
  if (!outcome_target_filters_available) {
    stop(paste(
      "Outcome target filters are unavailable.",
      "Run code/pathreat.data.outcome-sample-filters.R first."
    ))
  }
  base_outcome <- base_outcome_name(outcome)
  if (!base_outcome %in% names(dep_sample_filters)) {
    stop(sprintf("No treated-side target-filter definition for outcome: %s",
                 base_outcome))
  }
  dep_sample_filters[[base_outcome]]
}

build_matched_pair_index <- function(data,
                                     pair_id = "matched_pair_id",
                                     treatment = "treat") {
  required <- c(pair_id, treatment)
  missing_required <- setdiff(required, names(data))
  if (length(missing_required) > 0) {
    stop(sprintf(
      "Matched data are missing required columns: %s",
      paste(missing_required, collapse = ", ")
    ))
  }

  pair_ids <- data[[pair_id]]
  treat <- data[[treatment]]

  if (anyNA(pair_ids)) {
    stop("matched_pair_id contains missing values")
  }
  if (anyNA(treat) || !all(treat %in% c(0, 1))) {
    stop("treat must be observed and coded 0/1 for every matched row")
  }

  treated_rows <- which(treat == 1)
  control_rows_unordered <- which(treat == 0)
  treated_ids <- pair_ids[treated_rows]
  control_ids <- pair_ids[control_rows_unordered]

  if (anyDuplicated(treated_ids)) {
    stop("At least one matched_pair_id has more than one protected row")
  }
  if (anyDuplicated(control_ids)) {
    stop("At least one matched_pair_id has more than one control row")
  }
  if (length(treated_ids) != length(control_ids)) {
    stop("Protected and control row counts differ in the matched data")
  }

  control_order <- match(treated_ids, control_ids)
  if (anyNA(control_order) || length(treated_ids) * 2L != nrow(data)) {
    stop("Every matched_pair_id must contain exactly one protected and one control row")
  }

  structure(
    list(
      pair_ids = treated_ids,
      treated_rows = treated_rows,
      control_rows = control_rows_unordered[control_order],
      n_pairs = length(treated_ids),
      n_rows = nrow(data),
      pair_id = pair_id,
      treatment = treatment
    ),
    class = "pathreat_matched_pair_index"
  )
}

matched_pair_sample_mask <- function(data,
                                     outcome,
                                     pair_index = NULL,
                                     treated_filter = NULL,
                                     filter_outcome = outcome) {
  if (is.null(pair_index)) {
    pair_index <- build_matched_pair_index(data)
  }
  if (!inherits(pair_index, "pathreat_matched_pair_index") ||
      pair_index$n_rows != nrow(data)) {
    stop("pair_index is not valid for the supplied matched data")
  }

  current_pair_ids <- data[[pair_index$pair_id]]
  current_treatment <- data[[pair_index$treatment]]
  index_matches_data <-
    all(current_pair_ids[pair_index$treated_rows] == pair_index$pair_ids) &&
    all(current_pair_ids[pair_index$control_rows] == pair_index$pair_ids) &&
    all(current_treatment[pair_index$treated_rows] == 1) &&
    all(current_treatment[pair_index$control_rows] == 0)
  if (!index_matches_data) {
    stop("pair_index row positions do not match the supplied data order")
  }
  if (!outcome %in% names(data)) {
    stop(sprintf("Outcome column not found: %s", outcome))
  }

  treated_rows <- pair_index$treated_rows
  control_rows <- pair_index$control_rows

  if (is.null(treated_filter)) {
    treated_filter <- get_outcome_target_filter(filter_outcome)
  }

  if (is.null(treated_filter)) {
    target_eligible <- rep(TRUE, pair_index$n_pairs)
  } else if (is.character(treated_filter) || is.language(treated_filter)) {
    filter_call <- if (is.character(treated_filter)) {
      parse(text = treated_filter)[[1]]
    } else {
      treated_filter
    }
    filter_variables <- all.vars(filter_call)
    missing_filter_variables <- setdiff(filter_variables, names(data))
    if (length(missing_filter_variables) > 0) {
      stop(sprintf(
        "Treated-side target filter uses missing columns: %s",
        paste(missing_filter_variables, collapse = ", ")
      ))
    }
    filter_data <- lapply(filter_variables, function(variable) {
      data[[variable]][treated_rows]
    })
    names(filter_data) <- filter_variables
    target_eligible <- eval(filter_call, envir = filter_data)
  } else if (is.logical(treated_filter)) {
    if (length(treated_filter) == nrow(data)) {
      target_eligible <- treated_filter[treated_rows]
    } else {
      target_eligible <- treated_filter
    }
  } else {
    stop("treated_filter must be NULL, an expression, a string, or a logical vector")
  }

  if (!is.logical(target_eligible) ||
      length(target_eligible) != pair_index$n_pairs) {
    stop("Treated-side target filter must return one logical value per matched pair")
  }
  target_eligible[is.na(target_eligible)] <- FALSE

  outcome_values <- data[[outcome]]
  observed <- function(x) {
    if (is.numeric(x)) {
      is.finite(x)
    } else {
      !is.na(x)
    }
  }
  pair_eligible <- target_eligible &
    observed(outcome_values[treated_rows]) &
    observed(outcome_values[control_rows])

  row_mask <- rep(FALSE, nrow(data))
  row_mask[treated_rows[pair_eligible]] <- TRUE
  row_mask[control_rows[pair_eligible]] <- TRUE
  attr(row_mask, "n_pairs") <- sum(pair_eligible)
  attr(row_mask, "outcome") <- outcome
  row_mask
}

#########################
### parallel settings ###
#########################

# Number of worker processes for parallel sections. Override with the
# PATHREAT_CORES environment variable; otherwise use all but one detected
# core, capped at 8.
detected_cores <- suppressWarnings(parallel::detectCores())
if (is.na(detected_cores) || detected_cores < 1) {
  detected_cores <- 1L
}

no_cluster <- suppressWarnings(as.integer(Sys.getenv("PATHREAT_CORES", unset = "")))
if (is.na(no_cluster) || no_cluster < 1L) {
  no_cluster <- max(1L, min(8L, detected_cores - 1L))
}

#########################
### helper functions ####
#########################

f.histoden <- function(x, scale = 10^0, logplot = F, asinh = F) {
  print(summary(x))
  x <- x[!is.na(x)]
  if (logplot == F & asinh == F) {
    hist((x/scale), freq = F)
    lines(density((x/scale)), col = "orange")
  }
  if (logplot == T) {
    hist(log(x/scale + 1), freq = F)
    lines(density(log(x/scale) + 1), col = "orange")
  }
  if (asinh == T) {
    hist(asinh(x/scale), freq = F)
    lines(density(asinh(x/scale)), col = "orange")
  }
}

##################################################
### Canonical regression estimates and levels ####
##################################################

# The fitted model owns its covariance and degrees of freedom. Consumers
# retain their established formulas, samples, weights and covariance choices.
regression_estimate <- function(fit, term = "treat", level = 0.95) {
  ct <- fixest::coeftable(fit)[term, ]
  ci <- stats::confint(fit, parm = term, level = level)
  data.frame(
    coef = unname(ct[1]), se = unname(ct[2]), t_stat = unname(ct[3]),
    pval = unname(ct[4]), ci_low = unname(ci[1, 1]),
    ci_high = unname(ci[1, 2]), n_obs = stats::nobs(fit)
  )
}

# Average both treatment counterfactuals over the protected observations.
# This is deliberately model-based even where pair balance makes the levels
# numerically identical to the observed group means.
regression_protected_levels <- function(fit, data, treatment = "treat") {
  target <- data[data[[treatment]] == 1, , drop = FALSE]
  if (!nrow(target)) stop("No protected observations for standardization")
  target[[treatment]] <- 0
  control <- stats::predict(fit, newdata = target)
  target[[treatment]] <- 1
  protected <- stats::predict(fit, newdata = target)
  if (any(!is.finite(c(control, protected)))) {
    stop("Non-finite standardized regression levels")
  }
  c(adjusted_control = mean(control), adjusted_protected = mean(protected))
}

# QR uses explicit FE indicators when iterative absorption fails. Compute the
# treatment sandwich variance from this stable fit, retaining only the original
# model's FE parameter count, test degrees of freedom and small-sample rules.
# In particular, country indicators nested within country clusters must not
# acquire a different correction merely because they are written explicitly.
regression_qr_inference <- function(fit, reference_fit, data, n_countries) {
  reference <- if (n_countries > 1L) {
    summary(reference_fit, cluster = ~country_rast)
  } else summary(reference_fit, vcov = "hetero")
  correction <- attributes(stats::vcov(reference, attr = TRUE))
  n <- stats::nobs(fit)
  k <- correction$df.K
  df <- correction$df.t
  if (length(k) != 1L || length(df) != 1L ||
      !is.finite(k) || !is.finite(df) || n <= k || df <= 0) {
    stop("Invalid degrees of freedom for QR regression inference")
  }
  retained <- fit$qr$pivot[seq_len(fit$rank)]
  x <- stats::model.matrix(fit)[, retained, drop = FALSE]
  treatment_column <- match("treat", colnames(x))
  if (is.na(treatment_column)) stop("QR regression did not identify treat")
  r <- qr.R(fit$qr)[seq_len(fit$rank), seq_len(fit$rank), drop = FALSE]
  bread <- chol2inv(r)
  score <- as.vector(x %*% bread[, treatment_column]) * stats::residuals(fit)
  if (n_countries > 1L) {
    variance <- sum(rowsum(score, data$country_rast)^2)
    if (isTRUE(correction$ssc$G.adj)) {
      variance <- variance * n_countries / (n_countries - 1)
    }
    if (isTRUE(correction$ssc$K.adj)) variance <- variance * (n - 1) / (n - k)
  } else {
    variance <- sum(score^2)
    if (isTRUE(correction$ssc$K.adj)) variance <- variance * n / (n - k)
  }
  beta <- unname(stats::coef(fit)["treat"])
  se <- sqrt(variance)
  critical <- stats::qt(0.975, df = df)
  data.frame(
    coef = beta, se = se, t_stat = beta / se,
    pval = 2 * stats::pt(-abs(beta / se), df = df),
    ci_low = beta - critical * se, ci_high = beta + critical * se, n_obs = n
  )
}

# Complete-pair subgroup model shared by PA and species producers. The caller
# still owns membership and eligibility. Singleton OLS identifies the effect,
# but cannot supply residual-based uncertainty.
fit_pair_subgroup <- function(dd) {
  n_pairs <- sum(dd$treat == 1)
  if (nrow(dd) != 2L * n_pairs || n_pairs < 1L) {
    stop("Subgroup regression is not pair-complete")
  }
  n_countries <- length(unique(dd$country_rast))
  n_biomes <- length(unique(dd$biome_raster))
  fe <- c(if (n_countries > 1L) "country_rast",
          if (n_biomes > 1L) "biome_raster")
  fe_spec <- if (length(fe)) {
    paste(c(if (n_countries > 1L) "country", if (n_biomes > 1L) "biome"),
          collapse = "+")
  } else "none"
  constant <- length(unique(dd$threat_composite)) == 1L
  reference_fit <- NULL
  out <- data.frame(
    coef = NA_real_, se_delta = NA_real_, t_stat = NA_real_, p_value = NA_real_,
    ci_low = NA_real_, ci_high = NA_real_, n_obs = nrow(dd),
    n_biomes = n_biomes, fe_spec = fe_spec,
    se_type = if (n_countries > 1L) "cluster_country" else "hetero",
    status = "estimated", fit_engine = "fixest", inference_reason = NA_character_,
    adjusted_control = NA_real_, adjusted_protected = NA_real_
  )
  if (n_pairs == 1L || constant) {
    rhs <- c("treat", if (length(fe)) paste0("factor(", fe, ")"))
    fit <- stats::lm(stats::reformulate(rhs, response = "threat_composite"),
                     data = dd, na.action = na.fail)
    out$fit_engine <- "lm"
    out$status <- if (n_pairs == 1L) "single_pair" else "constant_outcome"
    out$inference_reason <- if (n_pairs == 1L) {
      "Saturated single-pair regression has no residual degrees of freedom"
    } else "Constant outcome: inferential statistics suppressed"
  } else {
    fml <- stats::as.formula(paste("threat_composite ~ treat",
      if (length(fe)) paste("|", paste(fe, collapse = " + "))))
    fit <- fixest::feols(fml, data = dd)
    if (isFALSE(fit$convStatus)) {
      reference_fit <- fit
      rhs <- c("treat", paste0("factor(", fe, ")"))
      fit <- stats::lm(stats::reformulate(rhs, response = "threat_composite"),
                       data = dd, na.action = na.fail)
      out$fit_engine <- "lm_qr_fallback"
      out$inference_reason <- paste(
        "Absorbed FE solver did not converge.",
        "Explicit-indicator QR regression retains the original covariance correction."
      )
    }
  }
  if (stats::nobs(fit) != nrow(dd)) stop("Subgroup regression removed observations")
  out$coef <- unname(stats::coef(fit)["treat"])
  if (!is.finite(out$coef)) stop("Subgroup regression did not identify treat")
  levels <- regression_protected_levels(fit, dd)
  out$adjusted_control <- levels["adjusted_control"]
  out$adjusted_protected <- levels["adjusted_protected"]
  if (abs(diff(levels) - out$coef) > 1e-8 * max(1, abs(out$coef))) {
    stop("Standardized level gap does not equal the regression coefficient")
  }
  if (out$status != "estimated") return(out)
  inference <- tryCatch({
    values <- if (!is.null(reference_fit)) {
      regression_qr_inference(fit, reference_fit, dd, n_countries)
    } else {
      fitted_summary <- if (n_countries > 1L) {
        summary(fit, cluster = ~country_rast)
      } else summary(fit, vcov = "hetero")
      regression_estimate(fitted_summary)
    }
    if (any(!is.finite(unlist(values))) || values$se <= 0 ||
        values$pval < 0 || values$pval > 1) {
      stop("Non-finite or degenerate regression inference")
    }
    values
  }, error = function(err) err)
  if (inherits(inference, "error")) {
    out$status <- "no_inference"
    out$inference_reason <- conditionMessage(inference)
  } else {
    out[1, c("se_delta", "t_stat", "p_value", "ci_low", "ci_high")] <-
      inference[1, c("se", "t_stat", "pval", "ci_low", "ci_high")]
  }
  out
}
