##########################################
### pathreat.analysis.sensitivity.R ####
### Sensitivity analysis #################
### E-values (VanderWeele & Ding, 2017) ##
### Rosenbaum bounds (Rosenbaum, 2002) ###
##########################################

library(dplyr)
library(EValue)

source("code/pathreat.analysis.config.R")

######################
### load data ########
######################
{
cat("Loading matched data...\n")
if (!"d.mbase" %in% ls()) {
  d.mbase <- load_matched_data()
}
cat(sprintf("  Loaded %s observations\n", format(nrow(d.mbase), big.mark = ",")))
}

pair_index <- build_matched_pair_index(d.mbase)
global_estimates <- readRDS(paths$est_global) %>%
  filter(estimate_type == "unscaled")
if (anyDuplicated(global_estimates$variable)) stop("Duplicated global ATT outcomes")

############################
### compute E-values #######
############################
{
cat("\n=== E-value sensitivity analysis ===\n")
cat("Method: Normalized regression ATT -> approximate RR -> E-value\n")
cat("Reference: VanderWeele & Ding (2017), Annals of Internal Medicine\n\n")

# use deplist and labels from config (stays in sync with main estimation)
sens_threats <- deplist
sens_labels  <- deplist_labels

evalue_results <- lapply(sens_threats, function(dv) {
  # dv<-"any_threat"
  cat("  ", dv, "\n")

  sample_mask <- matched_pair_sample_mask(d.mbase, dv, pair_index)
  d <- d.mbase[
    sample_mask,
    c("matched_pair_id", "treat", dv),
    drop = FALSE
  ]
  n_treated <- sum(d$treat == 1)
  n_control <- sum(d$treat == 0)
  if (n_treated != n_control || nrow(d) != 2L * n_treated) {
    stop(sprintf("%s: E-value sample is not pair-complete", dv))
  }

  y_t <- d[[dv]][d$treat == 1]
  y_c <- d[[dv]][d$treat == 0]
  if (length(y_t) != length(y_c) ||
      any(!is.finite(y_t)) || any(!is.finite(y_c))) {
    stop(sprintf("%s: E-value outcomes are not complete on both pair sides", dv))
  }

  # normalized mean difference: d = (mean_T - mean_C) / sqrt(var_T + var_C)
  # Preserve the sensitivity normalization and point-only E-value calculation.
  # The global regression now supplies the numerator; pair means are an audit.
  estimate <- global_estimates[global_estimates$variable == dv, , drop = FALSE]
  if (nrow(estimate) != 1L || estimate$n_obs != nrow(d) ||
      !is.finite(estimate$coef) ||
      abs(estimate$coef - (mean(y_t) - mean(y_c))) >
        1e-8 * max(1, abs(estimate$coef))) {
    stop("Global ATT disagrees with the E-value matched sample: ", dv)
  }
  norm_diff <- estimate$coef / sqrt(var(y_t) + var(y_c))

  # E-value via EValue package
  # evalues.OLS standardizes est by sd; passing sd=1 with already-standardized d
  # gives approximate RR = exp(0.91 * d) and corresponding E-value
  ev <- evalues.OLS(est = norm_diff, se = NA, sd = 1)

  data.frame(
    variable = dv,
    norm_diff = round(norm_diff, 3),
    risk_ratio = round(as.numeric(ev["RR", "point"]), 3),
    evalue = round(as.numeric(ev["E-values", "point"]), 3),
    n_treated = length(y_t),
    n_control = length(y_c),
    row.names = NULL
  )
})

evalue_results <- do.call(rbind, evalue_results)
rownames(evalue_results) <- NULL

cat("\nE-value results:\n")
print(evalue_results)
}

################################
### compute Rosenbaum bounds ###
################################
{
cat("\n=== Rosenbaum bounds sensitivity analysis ===\n")
cat("Method: Sign test with sensitivity parameter Gamma\n")
cat("Reference: Rosenbaum (2002), Observational Studies\n\n")

# sign test bounds for a given gamma (Rosenbaum, 2002)
sign_test_bounds <- function(n_positive, n_total, gamma) {
  if (gamma == 1) {
    pval <- 2 * pbinom(min(n_positive, n_total - n_positive), n_total, 0.5)
    return(c(pval, pval))
  }
  p_plus  <- gamma / (1 + gamma)
  p_minus <- 1 / (1 + gamma)
  pval_upper <- pbinom(n_positive, n_total, p_plus)
  pval_lower <- 1 - pbinom(n_positive - 1, n_total, p_minus)
  c(pval_upper, pval_lower)
}

GAMMA_SEQ <- seq(1, 3, by = 0.1)

rosenbaum_results <- lapply(sens_threats, function(dv) {
  cat("  ", dv, "\n")

  sample_mask <- matched_pair_sample_mask(d.mbase, dv, pair_index)
  d <- d.mbase[
    sample_mask,
    c("matched_pair_id", "treat", dv),
    drop = FALSE
  ]
  n_eligible_pairs <- sum(d$treat == 1)
  if (n_eligible_pairs != sum(d$treat == 0) ||
      nrow(d) != 2L * n_eligible_pairs) {
    stop(sprintf("%s: Rosenbaum sample is not pair-complete", dv))
  }

  # pair treated and control by matched_pair_id
  treated <- d %>% filter(treat == 1) %>% select(matched_pair_id, outcome_t = all_of(dv))
  control <- d %>% filter(treat == 0) %>% select(matched_pair_id, outcome_c = all_of(dv))
  pairs <- inner_join(treated, control, by = "matched_pair_id") %>%
    mutate(diff = outcome_t - outcome_c)
  if (nrow(pairs) != n_eligible_pairs || any(!is.finite(pairs$diff))) {
    stop(sprintf("%s: paired-difference construction lost eligible pairs", dv))
  }
  pairs <- pairs %>% filter(diff != 0)

  n_total <- nrow(pairs)
  n_positive <- sum(pairs$diff > 0)

  if (n_total == 0) {
    cat("    WARNING: No valid non-zero paired differences\n")
    return(data.frame(variable = dv, n_pairs = 0,
                      critical_gamma_05 = NA, critical_gamma_01 = NA,
                      row.names = NULL))
  }

  # test across gamma sequence
  pvals <- t(sapply(GAMMA_SEQ, function(g) sign_test_bounds(n_positive, n_total, g)))

  # find critical gamma: smallest gamma where upper p-value >= alpha
  crit_05_idx <- which(pvals[, 1] >= 0.05)[1]
  crit_01_idx <- which(pvals[, 1] >= 0.01)[1]
  critical_gamma_05 <- if (!is.na(crit_05_idx)) GAMMA_SEQ[crit_05_idx] else Inf
  critical_gamma_01 <- if (!is.na(crit_01_idx)) GAMMA_SEQ[crit_01_idx] else Inf

  cat(sprintf("    Pairs: %s | Critical Gamma(0.05): %.2f\n",
              format(n_total, big.mark = ","), critical_gamma_05))

  data.frame(
    variable = dv,
    n_pairs = n_total,
    critical_gamma_05 = critical_gamma_05,
    critical_gamma_01 = critical_gamma_01,
    row.names = NULL
  )
})

rosenbaum_results <- do.call(rbind, rosenbaum_results)
rownames(rosenbaum_results) <- NULL

cat("\nRosenbaum bounds results:\n")
print(rosenbaum_results)
}

######################
### save results #####
######################
{
# save E-value results
evalue_file <- paste0(paths$results_dir, "pathreat.sensitivity.evalues.Rds")
saveRDS(evalue_results, evalue_file)

# save Rosenbaum results
rosenbaum_file <- paste0(paths$results_dir, "pathreat.sensitivity.rosenbaum.Rds")
saveRDS(rosenbaum_results, rosenbaum_file)

# save Rosenbaum summary CSV
rosenbaum_csv <- paste0(paths$results_dir, "pathreat.sensitivity.rosenbaum_summary.csv")
write.csv(rosenbaum_results, rosenbaum_csv, row.names = FALSE)
}

###################################
### combined LaTeX table ##########
###################################
{
# merge E-value and Rosenbaum results
combined <- merge(evalue_results, rosenbaum_results[, c("variable", "critical_gamma_05", "critical_gamma_01")],
                  by = "variable", all.x = TRUE)
# restore deplist order
combined <- combined[match(deplist, combined$variable), ]

# formatting helpers
fmt_num <- function(x) {
  s <- sprintf("%.3f", abs(x))
  ifelse(x < 0, paste0("-", s), s)
}

fmt_gamma <- function(x) {
  ifelse(is.infinite(x) | x > 3, "{>3.0}", sprintf("%.1f", x))
}

# variable labels with units from config (scale.label from Data_summary.xlsx)
sens_scale_labels <- unit_label_tex(deplist_scale_labels)

# escape special LaTeX characters
tex_label <- function(x) { x <- gsub("%", "\\\\%", x); gsub("&", "\\\\&", x) }

tex_lines <- c(
  paste0("Variable & {Norm.\\ ATT} & {Risk ratio} & {E-value}",
         " & {$\\Gamma$ ($\\alpha{=}0.05$)} & {$\\Gamma$ ($\\alpha{=}0.01$)} \\\\"),
  "\\midrule"
)

for (i in seq_len(nrow(combined))) {
  tex_lines <- c(tex_lines,
    paste0(tex_label(sens_scale_labels[combined$variable[i]]), " & ",
           fmt_num(combined$norm_diff[i]), " & ",
           sprintf("%.3f", combined$risk_ratio[i]), " & ",
           sprintf("%.3f", combined$evalue[i]), " & ",
           fmt_gamma(combined$critical_gamma_05[i]), " & ",
           fmt_gamma(combined$critical_gamma_01[i]), " \\\\"))
}

## Note: \begin{tabular}, \toprule, \bottomrule, \\, \end{tabular} are in pathreat_03.tex
## Strip trailing \\ from last line (main file provides \input{table}\\)
tex_lines[length(tex_lines)] <- sub(" *\\\\\\\\$", "", tex_lines[length(tex_lines)])

writeLines(tex_lines, "pub/tables/tab.sensitivity.tex")
}

cat("\n=== Sensitivity analysis complete ===\n")
cat("Output Rds:    ", evalue_file, "\n")
cat("Output Rds:    ", rosenbaum_file, "\n")
cat("Output CSV:    ", rosenbaum_csv, "\n")
cat("Output LaTeX:  pub/tables/tab.sensitivity.tex\n")
