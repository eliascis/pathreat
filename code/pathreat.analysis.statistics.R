##################################################
### pathreat.analysis.statistics.R ###############
### Manuscript summary statistics ################
##################################################
#
# Collects headline statistics used in the manuscript abstract and main text.
# This script reads existing analysis outputs only; it does not re-estimate
# effects, regenerate tables, or modify figures.
#
# Inputs:
#   results/pathreat.matched-sample.summary.csv
#   results/pathreat.global.est.Rds
#   results/pathreat.bycountry.est.Rds
#   data/store/pathreat.analysis.byPA.est.Rds
#   data/store/pathreat.analysis.byPA-country-membership.Rds
#   results/pathreat.sfa.est.Rds
#   results/pathreat.sfa.country.est.Rds
#   results/pathreat.sfa.aggregates.Rds
#   results/pathreat.byspecies.est.Rds
#   results/pathreat.bytaxa.est.Rds
#   results/pathreat.analysis.hotspots.Rds
#
# Outputs:
#   results/pathreat.analysis.statistics.csv
#   results/pathreat.analysis.statistics.txt
#

library(dplyr)

source("code/pathreat.analysis.config.R")

fmt_num <- function(x, digits = 1) {
  formatC(x, format = "f", digits = digits, big.mark = ",")
}

fmt_int <- function(x) {
  format(as.integer(round(x)), scientific = FALSE, big.mark = ",", trim = TRUE)
}

add_stat <- function(rows, key, value, label, unit = "", digits = NA_integer_) {
  if (length(value) != 1L) {
    stop("Statistic must be scalar: ", key)
  }

  rows[[length(rows) + 1L]] <- data.frame(
    key = key,
    value = as.numeric(value),
    unit = unit,
    digits = digits,
    label = label,
    stringsAsFactors = FALSE
  )
  rows
}

require_columns <- function(x, required, label) {
  if (!is.data.frame(x)) {
    stop(label, " must be a data frame.")
  }

  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(label, " is missing required fields: ", paste(missing, collapse = ", "))
  }
}

read_rds_checked <- function(path, label) {
  if (!file.exists(path)) {
    stop("Missing ", label, ": ", path)
  }

  tryCatch(
    readRDS(path),
    error = function(err) {
      stop("Cannot read ", label, " (", path, "): ", conditionMessage(err))
    }
  )
}

is_count_vector <- function(x) {
  is.numeric(x) &&
    all(is.finite(x)) &&
    all(x >= 0) &&
    all(abs(x - round(x)) <= 1e-8)
}

assert_close <- function(observed, expected, label, tolerance = 1e-9) {
  if (length(observed) != 1L || length(expected) != 1L ||
      !is.finite(observed) || !is.finite(expected)) {
    stop(label, " must compare two finite scalars.")
  }

  scale <- max(1, abs(observed), abs(expected))
  if (abs(observed - expected) > tolerance * scale) {
    stop(
      label, " is inconsistent: observed ", format(observed, digits = 16),
      ", expected ", format(expected, digits = 16), "."
    )
  }
}

cat("=== Loading headline result objects ===\n")

matched_summary_path <- file.path(
  paths$results_dir,
  "pathreat.matched-sample.summary.csv"
)
if (!file.exists(matched_summary_path)) {
  stop("Missing matched-sample audit: ", matched_summary_path)
}
matched_summary <- tryCatch(
  read.csv(
    matched_summary_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  ),
  error = function(err) {
    stop(
      "Cannot read matched-sample audit (", matched_summary_path, "): ",
      conditionMessage(err)
    )
  }
)

global_est <- read_rds_checked(paths$est_global, "global estimates")
bycountry_est <- read_rds_checked(paths$est_bycountry, "by-country estimates")
pa_est <- read_rds_checked(
  "data/store/pathreat.analysis.byPA.est.Rds",
  "PA-level estimates"
)
pa_country_membership <- read_rds_checked(
  "data/store/pathreat.analysis.byPA-country-membership.Rds",
  "PA-country membership audit"
)
sfa_est <- read_rds_checked(paths$est_sfa, "PA-level SFA estimates")
sfa_country_est <- read_rds_checked(
  paths$est_sfa_country,
  "country-level SFA estimates"
)
sfa_aggregates <- read_rds_checked(
  file.path(paths$results_dir, "pathreat.sfa.aggregates.Rds"),
  "SFA aggregates"
)
sp_est <- read_rds_checked(
  file.path(paths$results_dir, "pathreat.byspecies.est.Rds"),
  "species estimates"
)
taxon_bundle <- read_rds_checked(
  file.path(paths$results_dir, "pathreat.bytaxa.est.Rds"),
  "taxon estimates"
)
hotspot_est <- read_rds_checked(
  file.path(paths$results_dir, "pathreat.analysis.hotspots.Rds"),
  "hotspot estimates"
)

if (!all(c("settings", "counts", "effects", "levels") %in% names(taxon_bundle))) {
  stop("Taxon bundle is missing required components.")
}

cat(sprintf("  Global estimates: %s rows\n", format(nrow(global_est), big.mark = ",")))
cat(sprintf("  Matched-sample audit: %s rows\n",
            format(nrow(matched_summary), big.mark = ",")))
cat(sprintf("  By-country estimates: %s rows\n",
            format(nrow(bycountry_est), big.mark = ",")))
cat(sprintf("  PA-level estimates: %s rows\n",
            format(nrow(pa_est), big.mark = ",")))
cat(sprintf("  SFA estimates: %s PAs, %s countries\n",
            format(nrow(sfa_est), big.mark = ","),
            format(nrow(sfa_country_est), big.mark = ",")))
cat(sprintf("  Species effects: %s rows\n", format(nrow(sp_est), big.mark = ",")))
cat(sprintf("  Hotspot estimates: %s rows\n",
            format(nrow(hotspot_est), big.mark = ",")))
cat(sprintf("  Centralized taxon effects: %s rows\n",
            format(nrow(taxon_bundle$effects), big.mark = ",")))

stats <- list()

#################################
### Matched-sample audit ########
#################################

matched_count_fields <- c(
  "observations",
  "matched_pairs",
  "countries",
  "protected_areas",
  "treated_pixels",
  "treated_cells_25km",
  "indonesia_observations",
  "indonesia_matched_pairs",
  "indonesia_protected_areas",
  "indonesia_treated_pixels",
  "indonesia_cells_25km",
  "manifest_expected_strata"
)
matched_logical_fields <- c(
  "indonesia_passes_balance",
  "indonesia_in_stage"
)
matched_required_fields <- c(
  "stage",
  matched_count_fields,
  "indonesia_avg_abs_smd",
  "balance_threshold",
  matched_logical_fields,
  "manifest_contract_version",
  "manifest_input_md5"
)

require_columns(
  matched_summary,
  matched_required_fields,
  "Matched-sample audit"
)
if (nrow(matched_summary) != 2L || anyDuplicated(matched_summary$stage) ||
    !setequal(matched_summary$stage, c("pre_balance", "post_balance"))) {
  stop(
    "Matched-sample audit must contain exactly one pre_balance row and ",
    "one post_balance row."
  )
}

matched_summary <- matched_summary[
  match(c("pre_balance", "post_balance"), matched_summary$stage),
  ,
  drop = FALSE
]
rownames(matched_summary) <- NULL

if (anyNA(matched_summary[, matched_required_fields, drop = FALSE])) {
  stop("Matched-sample audit contains missing required values.")
}
for (field in matched_count_fields) {
  if (!is_count_vector(matched_summary[[field]])) {
    stop("Matched-sample audit field is not a nonnegative count: ", field)
  }
}
for (field in matched_logical_fields) {
  if (!is.logical(matched_summary[[field]])) {
    stop("Matched-sample audit field must be logical: ", field)
  }
}
if (!is.numeric(matched_summary$indonesia_avg_abs_smd) ||
    any(!is.finite(matched_summary$indonesia_avg_abs_smd)) ||
    any(matched_summary$indonesia_avg_abs_smd < 0)) {
  stop("Indonesia average absolute SMD must be finite and nonnegative.")
}
if (!is.numeric(matched_summary$balance_threshold) ||
    any(!is.finite(matched_summary$balance_threshold)) ||
    any(matched_summary$balance_threshold <= 0)) {
  stop("Country-balance thresholds must be finite and positive.")
}
if (!is.character(matched_summary$manifest_contract_version) ||
    any(!nzchar(matched_summary$manifest_contract_version))) {
  stop("Manifest contract versions must be nonempty strings.")
}
if (!is.character(matched_summary$manifest_input_md5) ||
    any(!grepl("^[[:xdigit:]]{32}$", matched_summary$manifest_input_md5))) {
  stop("Manifest input fingerprints must be 32-character MD5 strings.")
}

constant_audit_fields <- c(
  "indonesia_avg_abs_smd",
  "balance_threshold",
  "indonesia_passes_balance",
  "manifest_contract_version",
  "manifest_expected_strata",
  "manifest_input_md5"
)
for (field in constant_audit_fields) {
  if (length(unique(matched_summary[[field]])) != 1L) {
    stop("Matched-sample audit field must agree across stages: ", field)
  }
}
if (any(abs(matched_summary$balance_threshold - 0.10) > 1e-12)) {
  stop("Matched-sample audit must retain the approved 0.10 balance threshold.")
}
if (any(matched_summary$manifest_contract_version != "pathreat_matching_v3") ||
    any(matched_summary$manifest_expected_strata != 668L)) {
  stop("Matched-sample audit must attest the 668-stratum v3 matching contract.")
}

pre_balance <- matched_summary[matched_summary$stage == "pre_balance", ]
post_balance <- matched_summary[matched_summary$stage == "post_balance", ]

stage_count_fields <- setdiff(
  matched_count_fields,
  "manifest_expected_strata"
)
if (any(
  unlist(post_balance[, stage_count_fields, drop = FALSE], use.names = FALSE) >
    unlist(pre_balance[, stage_count_fields, drop = FALSE], use.names = FALSE)
)) {
  stop("Post-balance counts cannot exceed their pre-balance counterparts.")
}

indonesia_count_fields <- c(
  "indonesia_observations",
  "indonesia_matched_pairs",
  "indonesia_protected_areas",
  "indonesia_treated_pixels",
  "indonesia_cells_25km"
)
for (i in seq_len(nrow(matched_summary))) {
  row <- matched_summary[i, , drop = FALSE]

  if (row$observations != 2 * row$matched_pairs) {
    stop(row$stage, ": observations must equal twice the matched-pair count.")
  }
  if (row$matched_pairs != row$treated_pixels) {
    stop(row$stage, ": matched pairs must equal distinct treated pixels.")
  }
  if (row$protected_areas > row$treated_pixels ||
      row$treated_cells_25km > row$treated_pixels) {
    stop(row$stage, ": PA/cell counts exceed treated pixels.")
  }
  if (row$indonesia_observations != 2 * row$indonesia_matched_pairs ||
      row$indonesia_matched_pairs != row$indonesia_treated_pixels) {
    stop(row$stage, ": Indonesia observation, pair, and pixel counts disagree.")
  }
  if (row$indonesia_protected_areas > row$indonesia_treated_pixels ||
      row$indonesia_cells_25km > row$indonesia_treated_pixels) {
    stop(row$stage, ": Indonesia PA/cell counts exceed treated pixels.")
  }
  if (row$indonesia_observations > row$observations ||
      row$indonesia_matched_pairs > row$matched_pairs ||
      row$indonesia_protected_areas > row$protected_areas ||
      row$indonesia_treated_pixels > row$treated_pixels ||
      row$indonesia_cells_25km > row$treated_cells_25km) {
    stop(row$stage, ": Indonesia counts exceed full-sample counts.")
  }

  indonesia_counts <- unlist(
    row[, indonesia_count_fields, drop = FALSE],
    use.names = FALSE
  )
  if (row$indonesia_in_stage && any(indonesia_counts <= 0)) {
    stop(row$stage, ": retained Indonesia counts must all be positive.")
  }
  if (!row$indonesia_in_stage && any(indonesia_counts != 0)) {
    stop(row$stage, ": excluded Indonesia counts must all be zero.")
  }
}

if (!pre_balance$indonesia_in_stage) {
  stop("Indonesia must be present before the country-balance filter.")
}
expected_indonesia_pass <-
  pre_balance$indonesia_avg_abs_smd <= pre_balance$balance_threshold
if (!identical(
      pre_balance$indonesia_passes_balance,
      expected_indonesia_pass
    )) {
  stop("Indonesia balance status disagrees with its average absolute SMD.")
}
if (!identical(
      post_balance$indonesia_in_stage,
      pre_balance$indonesia_passes_balance
    )) {
  stop("Indonesia post-balance retention disagrees with its balance status.")
}
if (post_balance$indonesia_in_stage) {
  for (field in indonesia_count_fields) {
    if (post_balance[[field]] != pre_balance[[field]]) {
      stop("Retained Indonesia count changed at the balance filter: ", field)
    }
  }
}

sample_stat_specs <- list(
  observations = c("observations", "observations"),
  pairs = c("matched_pairs", "pairs"),
  countries = c("countries", "countries"),
  protected_areas = c("protected_areas", "PAs"),
  treated_pixels = c("treated_pixels", "pixels"),
  cells_25km = c("treated_cells_25km", "cells")
)
for (stage_name in c("pre_balance", "post_balance")) {
  stage_row <- matched_summary[matched_summary$stage == stage_name, ]
  stage_label <- gsub("_", "-", stage_name)

  for (stat_name in names(sample_stat_specs)) {
    field <- sample_stat_specs[[stat_name]][1]
    unit <- sample_stat_specs[[stat_name]][2]
    stats <- add_stat(
      stats,
      paste0("matched_", stage_name, "_", stat_name),
      stage_row[[field]],
      paste0("Matched-sample ", stage_label, " ", gsub("_", " ", stat_name)),
      unit,
      0
    )
  }

  indonesia_sample_stats <- setdiff(names(sample_stat_specs), "countries")
  for (stat_name in indonesia_sample_stats) {
    field <- paste0("indonesia_", sample_stat_specs[[stat_name]][1])
    if (stat_name == "cells_25km") {
      field <- "indonesia_cells_25km"
    }
    stats <- add_stat(
      stats,
      paste0("indonesia_", stage_name, "_", stat_name),
      stage_row[[field]],
      paste0("Indonesia ", stage_label, " ", gsub("_", " ", stat_name)),
      sample_stat_specs[[stat_name]][2],
      0
    )
  }
}

stats <- add_stat(
  stats,
  "indonesia_avg_abs_smd",
  pre_balance$indonesia_avg_abs_smd,
  "Indonesia average absolute standardized mean difference before country filtering",
  "average absolute SMD",
  3
)
stats <- add_stat(
  stats,
  "indonesia_balance_threshold",
  pre_balance$balance_threshold,
  "Country-balance threshold applied to Indonesia average absolute SMD",
  "average absolute SMD",
  2
)
stats <- add_stat(
  stats,
  "indonesia_passes_balance",
  pre_balance$indonesia_passes_balance,
  "Indicator that Indonesia passes the country-balance threshold",
  "indicator",
  0
)
stats <- add_stat(
  stats,
  "indonesia_retained_post_balance",
  post_balance$indonesia_in_stage,
  "Indicator that Indonesia is retained in the post-balance matched sample",
  "indicator",
  0
)

#################################
### Global threat estimates #####
#################################

require_columns(
  global_est,
  c("variable", "estimate_type", "coef", "pval", "n_obs", "control_mean"),
  "Global estimates"
)

global_unscaled <- global_est %>%
  filter(estimate_type == "unscaled")

composite <- global_unscaled %>%
  filter(variable == "threat_composite")

any_threat <- global_unscaled %>%
  filter(variable == "any_threat")

individual <- global_unscaled %>%
  filter(variable %in% threat_list_raw)

if (nrow(composite) != 1L || nrow(any_threat) != 1L) {
  stop("Global results must contain one unscaled composite and any-threat row.")
}
if (nrow(individual) != length(threat_list_raw) ||
    anyDuplicated(individual$variable) ||
    !setequal(individual$variable, threat_list_raw)) {
  stop("Global results do not contain exactly one unscaled row per threat.")
}
individual_values <- individual[
  ,
  c("coef", "pval", "n_obs", "control_mean"),
  drop = FALSE
]
if (any(!is.finite(unlist(individual_values, use.names = FALSE))) ||
    any(individual$pval < 0 | individual$pval > 1) ||
    any(individual$n_obs <= 0)) {
  stop("At least one individual-threat estimate has invalid headline values.")
}
if (any(!is.finite(unlist(
  composite[, c("coef", "pval", "n_obs", "control_mean")],
  use.names = FALSE
))) || composite$pval < 0 || composite$pval > 1 ||
    composite$control_mean == 0) {
  stop("Global composite result contains invalid headline values.")
}
if (any(!is.finite(unlist(
  any_threat[, c("coef", "pval", "n_obs", "control_mean")],
  use.names = FALSE
))) || any_threat$pval < 0 || any_threat$pval > 1) {
  stop("Global any-threat result contains invalid headline values.")
}
if (composite$n_obs != post_balance$observations ||
    any_threat$n_obs != post_balance$observations) {
  stop("Global estimation samples disagree with the post-balance audit.")
}

lower_p05 <- individual %>%
  filter(coef < 0, pval < 0.05)

lower_p10 <- individual %>%
  filter(coef < 0, pval < 0.10)

stats <- add_stat(
  stats,
  "global_composite_percent_change",
  100 * composite$coef / composite$control_mean,
  "Composite mapped-pressure index percentage difference inside PAs relative to matched controls",
  "%",
  1
)

stats <- add_stat(
  stats,
  "global_composite_coef",
  composite$coef,
  "Composite mapped-pressure index ATT in raw index units",
  "index units",
  3
)

stats <- add_stat(
  stats,
  "global_composite_p_value",
  composite$pval,
  "Composite mapped-pressure index p-value",
  "",
  3
)

stats <- add_stat(
  stats,
  "global_any_threat_percentage_points",
  100 * any_threat$coef,
  "Any-threat indicator percentage-point difference inside PAs relative to matched controls",
  "percentage points",
  1
)

stats <- add_stat(
  stats,
  "global_any_threat_p_value",
  any_threat$pval,
  "Any-threat indicator p-value",
  "",
  3
)

stats <- add_stat(
  stats,
  "individual_threats_lower_p05",
  nrow(lower_p05),
  "Number of individual mapped threat indicators significantly lower inside PAs at p < 0.05",
  "indicators",
  0
)

stats <- add_stat(
  stats,
  "individual_threats_lower_p10",
  nrow(lower_p10),
  "Number of individual mapped threat indicators lower inside PAs at p < 0.10",
  "indicators",
  0
)

stats <- add_stat(
  stats,
  "global_composite_matched_pairs",
  composite$n_obs / 2,
  "Matched 1-km protected-control pairs in the composite-effect estimation sample",
  "pairs",
  0
)

#################################
### Country composite audit #####
#################################

country_required <- c(
  "country_id",
  "country_name",
  "variable",
  "coef",
  "pval",
  "n_obs",
  "n_treat",
  "n_control",
  "control_mean",
  "status"
)
require_columns(bycountry_est, country_required, "By-country estimates")

country_composite <- bycountry_est[
  bycountry_est$variable == "threat_composite",
  country_required,
  drop = FALSE
]
if (nrow(country_composite) != post_balance$countries ||
    anyDuplicated(country_composite$country_id) ||
    anyDuplicated(country_composite$country_name)) {
  stop(
    "Composite by-country rows must identify every post-balance country ",
    "exactly once."
  )
}
if (anyNA(country_composite[, c(
  "country_id", "country_name", "n_obs", "n_treat", "n_control", "status"
), drop = FALSE])) {
  stop("Composite by-country results contain missing identifiers or counts.")
}
for (field in c("n_obs", "n_treat", "n_control")) {
  if (!is_count_vector(country_composite[[field]])) {
    stop("Composite by-country result has an invalid count field: ", field)
  }
}
if (any(
  country_composite$n_obs !=
    country_composite$n_treat + country_composite$n_control
) || any(country_composite$n_treat != country_composite$n_control)) {
  stop("Composite by-country samples are not pair-complete.")
}
if (sum(country_composite$n_obs) != post_balance$observations ||
    sum(country_composite$n_treat) != post_balance$matched_pairs) {
  stop("Composite by-country sample counts disagree with the matched audit.")
}

country_estimated <- country_composite$status == "estimated"
if (anyNA(country_estimated)) {
  stop("Composite by-country status contains missing values.")
}
estimated_values <- country_composite[
  country_estimated,
  c("coef", "pval", "control_mean"),
  drop = FALSE
]
if (nrow(estimated_values) > 0L &&
    any(!is.finite(unlist(estimated_values, use.names = FALSE)))) {
  stop("Estimated composite by-country rows contain nonfinite values.")
}
if (nrow(estimated_values) > 0L &&
    any(estimated_values$pval < 0 | estimated_values$pval > 1)) {
  stop("Estimated composite by-country p-values lie outside [0, 1].")
}

country_sig_reduction <- country_estimated &
  country_composite$pval < 0.05 & country_composite$coef < 0
country_sig_increase <- country_estimated &
  country_composite$pval < 0.05 & country_composite$coef > 0
country_null <- country_estimated & country_composite$pval >= 0.05
country_unclassified <- country_estimated & !(
  country_sig_reduction | country_sig_increase | country_null
)
if (any(country_unclassified)) {
  stop("An estimated composite country effect has no valid direction class.")
}
if (sum(country_sig_reduction) + sum(country_sig_increase) +
    sum(country_null) != sum(country_estimated)) {
  stop("Composite country-effect classifications are not exhaustive.")
}

indonesia_country <- country_composite$country_id == 102 |
  tolower(trimws(country_composite$country_name)) == "indonesia"
if (sum(indonesia_country) > 1L) {
  stop("Composite by-country results contain multiple Indonesia rows.")
}
if (any(indonesia_country) && !(
  country_composite$country_id[indonesia_country] == 102 &&
    tolower(trimws(country_composite$country_name[indonesia_country])) ==
      "indonesia"
)) {
  stop("Indonesia country ID and name disagree in by-country results.")
}
if (any(indonesia_country) != post_balance$indonesia_in_stage) {
  stop("Indonesia by-country presence disagrees with matched-sample retention.")
}

indonesia_country_estimated <- FALSE
indonesia_country_coef <- NA_real_
indonesia_country_pval <- NA_real_
indonesia_country_percent <- NA_real_
indonesia_country_direction <- NA_real_
if (any(indonesia_country)) {
  indonesia_row <- country_composite[indonesia_country, , drop = FALSE]

  if (indonesia_row$n_obs != post_balance$indonesia_observations ||
      indonesia_row$n_treat != post_balance$indonesia_matched_pairs) {
    stop("Indonesia by-country sample counts disagree with the matched audit.")
  }

  indonesia_country_estimated <- indonesia_row$status == "estimated"
  if (!indonesia_country_estimated) {
    stop("Retained Indonesia lacks a valid composite by-country estimate.")
  }

  indonesia_country_coef <- indonesia_row$coef
  indonesia_country_pval <- indonesia_row$pval
  if (indonesia_row$control_mean != 0) {
    indonesia_country_percent <-
      100 * indonesia_row$coef / indonesia_row$control_mean
  }
  indonesia_country_direction <- if (
    indonesia_row$pval >= 0.05
  ) {
    0
  } else if (indonesia_row$coef < 0) {
    -1
  } else {
    1
  }
}

stats <- add_stat(
  stats,
  "country_composite_estimated_n",
  sum(country_estimated),
  "Countries with an estimated composite mapped-pressure effect",
  "countries",
  0
)
stats <- add_stat(
  stats,
  "country_composite_significant_reduction_n",
  sum(country_sig_reduction),
  "Countries with a significant composite mapped-pressure reduction at p < 0.05",
  "countries",
  0
)
stats <- add_stat(
  stats,
  "country_composite_significant_increase_n",
  sum(country_sig_increase),
  "Countries with a significant composite mapped-pressure increase at p < 0.05",
  "countries",
  0
)
stats <- add_stat(
  stats,
  "country_composite_null_n",
  sum(country_null),
  "Countries without a significant composite mapped-pressure effect at p < 0.05",
  "countries",
  0
)
stats <- add_stat(
  stats,
  "country_composite_not_estimated_n",
  sum(!country_estimated),
  "Post-balance countries without an estimated composite mapped-pressure effect",
  "countries",
  0
)
stats <- add_stat(
  stats,
  "country_composite_indonesia_present",
  any(indonesia_country),
  "Indicator that Indonesia is present in composite by-country results",
  "indicator",
  0
)
stats <- add_stat(
  stats,
  "country_composite_indonesia_estimated",
  indonesia_country_estimated,
  "Indicator that Indonesia has an estimated composite by-country effect",
  "indicator",
  0
)
stats <- add_stat(
  stats,
  "country_composite_indonesia_coef",
  indonesia_country_coef,
  "Indonesia composite mapped-pressure ATT in raw index units",
  "index units",
  3
)
stats <- add_stat(
  stats,
  "country_composite_indonesia_percent_change",
  indonesia_country_percent,
  "Indonesia composite mapped-pressure percentage difference relative to matched controls",
  "%",
  1
)
stats <- add_stat(
  stats,
  "country_composite_indonesia_p_value",
  indonesia_country_pval,
  "Indonesia composite mapped-pressure p-value",
  "",
  3
)
stats <- add_stat(
  stats,
  "country_composite_indonesia_direction",
  indonesia_country_direction,
  "Indonesia effect direction (-1 significant reduction, 0 null, 1 significant increase)",
  "category",
  0
)

#################################
### PA and SFA audits ############
#################################

pa_required <- c("wdpaid", "n_pixels", "country_rast", "n_countries",
                 "coef", "mean_delta", "mean_pair_difference", "n_obs", "n_biomes",
                 "fe_spec", "se_type", "status", "fit_engine", "effectiveness")
require_columns(pa_est, pa_required, "PA-level estimates")
if (nrow(pa_est) == 0L || nrow(pa_est) != post_balance$protected_areas ||
    anyNA(pa_est[, pa_required, drop = FALSE]) ||
    anyDuplicated(pa_est$wdpaid)) {
  stop("PA-level results must identify every post-balance PA exactly once.")
}
for (field in c("n_pixels", "n_countries", "n_obs", "n_biomes")) {
  if (!is_count_vector(pa_est[[field]]) || any(pa_est[[field]] <= 0)) {
    stop("PA-level result has an invalid count field: ", field)
  }
}
if (sum(pa_est$n_pixels) != post_balance$matched_pairs) {
  stop("PA-level pair counts disagree with the post-balance matched sample.")
}

require_columns(pa_est, c("se_delta", "t_stat", "p_value", "ci_low", "ci_high",
                          "inference_reason"), "PA regression inference")
pa_has_inference <- pa_est$status == "estimated"
pa_inference_fields <- c("se_delta", "t_stat", "p_value", "ci_low", "ci_high")
expected_fe <- ifelse(pa_est$n_countries > 1,
                      ifelse(pa_est$n_biomes > 1, "country+biome", "country"),
                      ifelse(pa_est$n_biomes > 1, "biome", "none"))
if (any(!is.finite(pa_est$coef)) || any(pa_est$n_obs != 2 * pa_est$n_pixels) ||
    any(abs(pa_est$coef - pa_est$mean_pair_difference) > 1e-8) ||
    any(pa_est$mean_delta != pa_est$coef) || any(pa_est$fe_spec != expected_fe) ||
    any(pa_est$se_type != ifelse(pa_est$n_countries > 1, "cluster_country", "hetero")) ||
    any(!pa_est$status %in% c("estimated", "single_pair", "constant_outcome", "no_inference")) ||
    any(!is.finite(as.matrix(pa_est[pa_has_inference, pa_inference_fields]))) ||
    any(!is.na(as.matrix(pa_est[!pa_has_inference, pa_inference_fields]))) ||
    any(pa_est$p_value[pa_has_inference] < 0 | pa_est$p_value[pa_has_inference] > 1) ||
    any(pa_est$se_delta[pa_has_inference] <= 0) ||
    any(is.na(pa_est$inference_reason[!pa_has_inference])) ||
    any((pa_est$effectiveness == "no_inference") != !pa_has_inference)) {
  stop("PA regression coefficient, specification, or inference audit failed")
}

pa_country_required <- c("wdpaid", "country_rast", "n_pixels_country")
require_columns(
  pa_country_membership,
  pa_country_required,
  "PA-country membership audit"
)
if (nrow(pa_country_membership) == 0L ||
    anyNA(pa_country_membership[, pa_country_required, drop = FALSE]) ||
    anyDuplicated(pa_country_membership[c("wdpaid", "country_rast")]) ||
    !is_count_vector(pa_country_membership$n_pixels_country) ||
    any(pa_country_membership$n_pixels_country <= 0) ||
    sum(pa_country_membership$n_pixels_country) !=
      post_balance$matched_pairs ||
    n_distinct(pa_country_membership$wdpaid) !=
      post_balance$protected_areas) {
  stop("PA-country membership audit is incomplete or malformed.")
}

pa_country_check <- pa_country_membership %>%
  group_by(wdpaid) %>%
  summarise(
    n_pixels = sum(n_pixels_country),
    country_rast = min(country_rast),
    n_countries = n(),
    .groups = "drop"
  )
if (!setequal(pa_est$wdpaid, pa_country_check$wdpaid)) {
  stop("PA-level results and the PA-country membership audit disagree.")
}
pa_country_index <- match(pa_est$wdpaid, pa_country_check$wdpaid)
for (field in c("n_pixels", "country_rast", "n_countries")) {
  if (any(pa_est[[field]] != pa_country_check[[field]][pa_country_index])) {
    stop("PA-level result disagrees with country membership field: ", field)
  }
}

indonesia_pa_country <- pa_country_membership$country_rast == 102
if (sum(indonesia_pa_country) !=
      post_balance$indonesia_protected_areas ||
    sum(pa_country_membership$n_pixels_country[indonesia_pa_country]) !=
      post_balance$indonesia_matched_pairs) {
  stop("Indonesia PA-country membership disagrees with the matched audit.")
}

stats <- add_stat(
  stats,
  "pa_results_n",
  nrow(pa_est),
  "PA records in the canonical regression output",
  "PAs",
  0
)

stats <- add_stat(stats, "pa_coefficients_n", sum(is.finite(pa_est$coef)),
                  "PAs with fitted regression coefficients", "PAs", 0)
stats <- add_stat(stats, "pa_inference_n", sum(pa_has_inference),
                  "PAs with available regression inference", "PAs", 0)
stats <- add_stat(stats, "pa_no_inference_n", sum(!pa_has_inference),
                  "PAs without regression inference", "PAs", 0)
for (status_name in c("estimated", "single_pair", "constant_outcome", "no_inference")) {
  stats <- add_stat(stats, paste0("pa_status_", status_name, "_n"),
                    sum(pa_est$status == status_name),
                    paste("PA estimation status:", status_name), "PAs", 0)
}

sfa_pa_required <- c(
  "wdpaid",
  "country_rast",
  "n_countries",
  "n_pixels",
  "delta",
  "tc_protected",
  "frontier_level",
  "frontier",
  "u_hat",
  "size"
)
require_columns(sfa_est, sfa_pa_required, "PA-level SFA estimates")
if (nrow(sfa_est) == 0L || anyNA(sfa_est[, sfa_pa_required, drop = FALSE]) ||
    anyDuplicated(sfa_est$wdpaid)) {
  stop("PA-level SFA results must contain unique, complete PA rows.")
}
if (any(!is.finite(unlist(
  sfa_est[, setdiff(sfa_pa_required, "wdpaid"), drop = FALSE],
  use.names = FALSE
))) || any(sfa_est$n_pixels <= 0) || any(sfa_est$size <= 0)) {
  stop("PA-level SFA results contain invalid numeric values.")
}
if (!all(sfa_est$wdpaid %in% pa_est$wdpaid)) {
  stop("SFA results contain PAs absent from PA-level results.")
}

pa_sfa_index <- match(sfa_est$wdpaid, pa_est$wdpaid)
for (field in c("country_rast", "n_countries", "n_pixels")) {
  if (any(sfa_est[[field]] != pa_est[[field]][pa_sfa_index])) {
    stop("SFA and PA-level results disagree on field: ", field)
  }
}

sfa_country_required <- c(
  "country_rast",
  "n_pas",
  "n_pixels",
  "delta_mean",
  "frontier_mean"
)
require_columns(
  sfa_country_est,
  sfa_country_required,
  "Country-level SFA estimates"
)
if (nrow(sfa_country_est) == 0L ||
    anyNA(sfa_country_est[, sfa_country_required, drop = FALSE]) ||
    anyDuplicated(sfa_country_est$country_rast)) {
  stop("Country-level SFA results must contain unique, complete country rows.")
}
for (field in c("n_pas", "n_pixels")) {
  if (!is_count_vector(sfa_country_est[[field]]) ||
      any(sfa_country_est[[field]] <= 0)) {
    stop("Country-level SFA result has an invalid count field: ", field)
  }
}
if (any(!is.finite(sfa_country_est$delta_mean)) ||
    any(!is.finite(sfa_country_est$frontier_mean)) ||
    any(abs(sfa_country_est$frontier_mean) <= .Machine$double.eps)) {
  stop("Country-level SFA achieved-share inputs must be finite and nonzero.")
}
sfa_country_achieved_share <-
  sfa_country_est$delta_mean / sfa_country_est$frontier_mean
if (any(!is.finite(sfa_country_achieved_share))) {
  stop("Country-level SFA achieved shares must be finite.")
}

sfa_country_check <- pa_country_membership %>%
  filter(wdpaid %in% sfa_est$wdpaid) %>%
  group_by(country_rast) %>%
  summarise(
    n_pas = n_distinct(wdpaid),
    n_pixels = sum(n_pixels_country),
    .groups = "drop"
  )
if (!setequal(sfa_country_est$country_rast, sfa_country_check$country_rast)) {
  stop("Country-level and PA-level SFA country sets disagree.")
}
sfa_country_index <- match(
  sfa_country_est$country_rast,
  sfa_country_check$country_rast
)
if (any(
  sfa_country_est$n_pas != sfa_country_check$n_pas[sfa_country_index]
) || any(
  sfa_country_est$n_pixels != sfa_country_check$n_pixels[sfa_country_index]
)) {
  stop("Country-level SFA PA/pixel counts disagree with PA-level SFA results.")
}
if (sum(sfa_country_est$n_pixels) != sum(sfa_est$n_pixels)) {
  stop("Country-level SFA pixel counts do not preserve the global sample.")
}

indonesia_sfa_country <- sfa_country_est$country_rast == 102
if (post_balance$indonesia_in_stage) {
  if (sum(indonesia_sfa_country) != 1L ||
      sfa_country_est$n_pas[indonesia_sfa_country] !=
        post_balance$indonesia_protected_areas ||
      sfa_country_est$n_pixels[indonesia_sfa_country] !=
        post_balance$indonesia_matched_pairs) {
    stop("Indonesia country-level SFA counts disagree with the matched audit.")
  }
} else if (any(indonesia_sfa_country)) {
  stop("Country-level SFA output retains balance-excluded Indonesia.")
}

sfa_aggregate_required <- c(
  "n_pas",
  "total_pa_area_km2",
  "aggregate_achieved_threat_reduction_threat_index_km2",
  "aggregate_underperformance_threat_index_km2",
  "aggregate_total_potential_threat_reduction_threat_index_km2",
  "share_achieved",
  "share_underperformance",
  "aggregate_net_residual_threat_index_km2",
  "aggregate_positive_residual_threat_index_km2",
  "area_weighted_mean_underperformance",
  "area_weighted_mean_net_residual",
  "area_weighted_mean_positive_residual",
  "median_underperformance",
  "median_net_residual",
  "share_negative_residual"
)
require_columns(sfa_aggregates, sfa_aggregate_required, "SFA aggregates")
if (nrow(sfa_aggregates) != 1L ||
    any(!is.finite(unlist(
      sfa_aggregates[, sfa_aggregate_required, drop = FALSE],
      use.names = FALSE
    )))) {
  stop("SFA aggregates must contain one complete finite row.")
}

sfa_epsilon <- sfa_est$tc_protected - sfa_est$frontier_level
if (max(abs(sfa_epsilon - (sfa_est$delta - sfa_est$frontier))) > 1e-10) {
  stop("SFA level- and delta-space residuals disagree.")
}
sfa_total_area <- sum(sfa_est$size)
sfa_achieved <- -sum(sfa_est$size * sfa_est$delta)
sfa_underperformance <- sum(sfa_est$size * sfa_est$u_hat)
sfa_total_potential <- sfa_achieved + sfa_underperformance
sfa_net_residual <- sum(sfa_est$size * sfa_epsilon)
sfa_positive_residual <- sum(sfa_est$size * pmax(sfa_epsilon, 0))
if (!is.finite(sfa_total_area) || sfa_total_area <= 0 ||
    !is.finite(sfa_total_potential) || sfa_total_potential <= 0) {
  stop("SFA total area and potential reduction must be finite and positive.")
}

sfa_recomputed <- c(
  n_pas = nrow(sfa_est),
  total_pa_area_km2 = sfa_total_area,
  aggregate_achieved_threat_reduction_threat_index_km2 = sfa_achieved,
  aggregate_underperformance_threat_index_km2 = sfa_underperformance,
  aggregate_total_potential_threat_reduction_threat_index_km2 =
    sfa_total_potential,
  share_achieved = sfa_achieved / sfa_total_potential,
  share_underperformance = sfa_underperformance / sfa_total_potential,
  aggregate_net_residual_threat_index_km2 = sfa_net_residual,
  aggregate_positive_residual_threat_index_km2 = sfa_positive_residual,
  area_weighted_mean_underperformance = sfa_underperformance / sfa_total_area,
  area_weighted_mean_net_residual = sfa_net_residual / sfa_total_area,
  area_weighted_mean_positive_residual = sfa_positive_residual / sfa_total_area,
  median_underperformance = median(sfa_est$u_hat),
  median_net_residual = median(sfa_epsilon),
  share_negative_residual = mean(sfa_epsilon < 0)
)
for (field in names(sfa_recomputed)) {
  assert_close(
    sfa_aggregates[[field]],
    sfa_recomputed[[field]],
    paste0("SFA aggregate ", field)
  )
}
assert_close(
  sfa_aggregates$share_achieved + sfa_aggregates$share_underperformance,
  1,
  "SFA achieved and underperformance shares"
)

stats <- add_stat(
  stats,
  "sfa_pa_n",
  nrow(sfa_est),
  "Protected areas in the stochastic-frontier results",
  "PAs",
  0
)
stats <- add_stat(
  stats,
  "sfa_country_n",
  nrow(sfa_country_est),
  "Countries in the stochastic-frontier country results",
  "countries",
  0
)

sfa_stat_specs <- list(
  total_pa_area_km2 = c("total_pa_area_km2", "km2", 1),
  aggregate_achieved_threat_reduction = c(
    "aggregate_achieved_threat_reduction_threat_index_km2",
    "threat-index km2",
    1
  ),
  aggregate_underperformance = c(
    "aggregate_underperformance_threat_index_km2",
    "threat-index km2",
    1
  ),
  aggregate_total_potential_reduction = c(
    "aggregate_total_potential_threat_reduction_threat_index_km2",
    "threat-index km2",
    1
  ),
  aggregate_net_residual = c(
    "aggregate_net_residual_threat_index_km2",
    "threat-index km2",
    1
  ),
  aggregate_positive_residual = c(
    "aggregate_positive_residual_threat_index_km2",
    "threat-index km2",
    1
  ),
  area_weighted_mean_underperformance = c(
    "area_weighted_mean_underperformance",
    "index units",
    4
  ),
  area_weighted_mean_net_residual = c(
    "area_weighted_mean_net_residual",
    "index units",
    4
  ),
  area_weighted_mean_positive_residual = c(
    "area_weighted_mean_positive_residual",
    "index units",
    4
  ),
  median_underperformance = c(
    "median_underperformance",
    "index units",
    3
  ),
  median_net_residual = c("median_net_residual", "index units", 3)
)
for (stat_name in names(sfa_stat_specs)) {
  spec <- sfa_stat_specs[[stat_name]]
  stats <- add_stat(
    stats,
    paste0("sfa_", stat_name),
    sfa_aggregates[[spec[1]]],
    paste0("SFA ", gsub("_", " ", stat_name)),
    spec[2],
    as.integer(spec[3])
  )
}
stats <- add_stat(
  stats,
  "sfa_share_achieved_percent",
  100 * sfa_aggregates$share_achieved,
  "Area-weighted share of total potential threat reduction achieved",
  "%",
  1
)
stats <- add_stat(
  stats,
  "sfa_share_underperformance_percent",
  100 * sfa_aggregates$share_underperformance,
  "Area-weighted share of total potential threat reduction attributed to underperformance",
  "%",
  1
)
stats <- add_stat(
  stats,
  "sfa_negative_residual_percent",
  100 * sfa_aggregates$share_negative_residual,
  "Share of PAs with a negative untruncated SFA frontier residual",
  "%",
  1
)
stats <- add_stat(
  stats,
  "sfa_country_achieved_share_median_percent",
  100 * median(sfa_country_achieved_share),
  "Median country ratio of observed to frontier-implied threat reduction",
  "%",
  1
)
stats <- add_stat(
  stats,
  "sfa_country_achieved_share_p25_percent",
  100 * unname(quantile(sfa_country_achieved_share, 0.25)),
  "25th percentile of the country achieved-share ratio",
  "%",
  1
)
stats <- add_stat(
  stats,
  "sfa_country_achieved_share_p75_percent",
  100 * unname(quantile(sfa_country_achieved_share, 0.75)),
  "75th percentile of the country achieved-share ratio",
  "%",
  1
)

#################################
### Species coverage stats ######
#################################

species_required <- c(
  "species_id",
  "taxon",
  "share_protected",
  "delta_s",
  "se_delta_s",
  "p_delta_s",
  "pct_effect_s",
  "control_mean_s",
  "status"
)
require_columns(sp_est, species_required, "Species estimates")
if (anyDuplicated(sp_est$species_id) || anyNA(sp_est$status)) {
  stop("Species estimates must have unique IDs and nonmissing statuses.")
}

canonical_species_numeric <- c(
  "delta_s",
  "se_delta_s",
  "p_delta_s",
  "pct_effect_s",
  "control_mean_s"
)
species_numeric_finite <- Reduce(
  `&`,
  lapply(
    canonical_species_numeric,
    function(field) is.finite(sp_est[[field]])
  )
)
species_estimated <- sp_est$status == "estimated"
if (any(species_estimated & !species_numeric_finite)) {
  stop("Canonical estimated-species rows contain nonfinite model values.")
}
species_eligible <- species_estimated & species_numeric_finite
if (!any(species_eligible)) {
  stop("No canonical species are eligible for taxon summaries.")
}

taxon_count_required <- c(
  "taxon",
  "n_species_estimated",
  "n_species_equal"
)
require_columns(taxon_bundle$counts, taxon_count_required, "Taxon counts")
taxon_effect_required <- c(
  "taxon",
  "weighting",
  "n_species",
  "pct_median"
)
require_columns(taxon_bundle$effects, taxon_effect_required, "Taxon effects")

effect_taxa <- taxon_bundle$settings$effect_taxon_order
if (!is.character(effect_taxa) || length(effect_taxa) == 0L ||
    anyNA(effect_taxa) || anyDuplicated(effect_taxa)) {
  stop("Taxon settings contain an invalid effect-taxon order.")
}
if (nrow(taxon_bundle$counts) != length(effect_taxa) ||
    anyDuplicated(taxon_bundle$counts$taxon) ||
    !setequal(taxon_bundle$counts$taxon, effect_taxa)) {
  stop("Taxon counts must contain exactly one row per configured taxon.")
}
for (field in c("n_species_estimated", "n_species_equal")) {
  if (!is_count_vector(taxon_bundle$counts[[field]])) {
    stop("Taxon count field is invalid: ", field)
  }
}
if (any(
  taxon_bundle$counts$n_species_estimated !=
    taxon_bundle$counts$n_species_equal
)) {
  stop("Equal-weight taxon counts disagree with estimated-species counts.")
}

eligible_taxa <- sp_est$taxon[species_eligible]
if (anyNA(eligible_taxa) || !all(eligible_taxa %in% effect_taxa)) {
  stop("Eligible canonical species contain missing or unsupported taxa.")
}
canonical_taxon_counts <- table(factor(eligible_taxa, levels = effect_taxa))
count_index <- match(effect_taxa, taxon_bundle$counts$taxon)
observed_taxon_counts <-
  taxon_bundle$counts$n_species_estimated[count_index]
if (!identical(
      as.integer(observed_taxon_counts),
      as.integer(canonical_taxon_counts)
    )) {
  stop("Per-taxon counts disagree with canonical eligible species.")
}
if (sum(observed_taxon_counts) != sum(species_eligible)) {
  stop("Summed taxon counts disagree with the canonical eligible-species count.")
}

taxon_equal <- taxon_bundle$effects[
  taxon_bundle$effects$weighting == "equal",
  taxon_effect_required,
  drop = FALSE
]
if (nrow(taxon_equal) != length(effect_taxa) ||
    anyDuplicated(taxon_equal$taxon) ||
    !setequal(taxon_equal$taxon, effect_taxa)) {
  stop("Equal-weight taxon effects must contain one row per configured taxon.")
}
taxon_effect_index <- match(effect_taxa, taxon_equal$taxon)
if (!is_count_vector(taxon_equal$n_species) || any(
  taxon_equal$n_species[taxon_effect_index] != observed_taxon_counts
) || any(!is.finite(taxon_equal$pct_median))) {
  stop("Equal-weight taxon effects disagree with canonical taxon counts.")
}
canonical_taxon_medians <- vapply(
  effect_taxa,
  function(taxon_name) {
    median(sp_est$pct_effect_s[
      species_eligible & sp_est$taxon == taxon_name
    ])
  },
  numeric(1)
)
observed_taxon_medians <- taxon_equal$pct_median[taxon_effect_index]
median_scale <- pmax(
  1,
  abs(canonical_taxon_medians),
  abs(observed_taxon_medians)
)
if (any(
  abs(observed_taxon_medians - canonical_taxon_medians) >
    1e-12 * median_scale
)) {
  stop("Equal-weight taxon medians disagree with canonical species effects.")
}

sp_cover <- sp_est %>%
  filter(is.finite(share_protected)) %>%
  mutate(
    share_unprotected = 1 - share_protected,
    no_pa_overlap = share_protected == 0
  ) %>%
  select(species_id, share_protected, share_unprotected, no_pa_overlap)
if (nrow(sp_cover) == 0L || any(
  sp_cover$share_protected < 0 | sp_cover$share_protected > 1
)) {
  stop("Species coverage shares must be finite probabilities.")
}

cat(sprintf("  Species coverage: %s rows\n", format(nrow(sp_cover), big.mark = ",")))

stats <- add_stat(
  stats,
  "species_coverage_n",
  nrow(sp_cover),
  "Threatened terrestrial vertebrate species in the PA-coverage analysis",
  "species",
  0
)

stats <- add_stat(
  stats,
  "species_no_pa_overlap_n",
  sum(sp_cover$no_pa_overlap, na.rm = TRUE),
  "Threatened terrestrial vertebrate species with no mapped PA overlap",
  "species",
  0
)

stats <- add_stat(
  stats,
  "species_no_pa_overlap_percent",
  mean(sp_cover$no_pa_overlap, na.rm = TRUE) * 100,
  "Share of threatened terrestrial vertebrate species with no mapped PA overlap",
  "%",
  1
)

stats <- add_stat(
  stats,
  "species_mean_range_protected_percent",
  mean(sp_cover$share_protected, na.rm = TRUE) * 100,
  "Mean share of threatened terrestrial vertebrate mapped range inside any PA",
  "%",
  1
)

stats <- add_stat(
  stats,
  "species_mean_range_unprotected_percent",
  mean(sp_cover$share_unprotected, na.rm = TRUE) * 100,
  "Mean share of threatened terrestrial vertebrate mapped range outside PAs",
  "%",
  1
)

#################################
### Species range-effect stats ##
#################################

sp_range <- sp_est[species_eligible, , drop = FALSE]

stats <- add_stat(
  stats,
  "species_range_effect_estimable_n",
  nrow(sp_range),
  "Threatened terrestrial vertebrate species with estimable range-based pressure contrasts",
  "species",
  0
)

stats <- add_stat(
  stats,
  "species_range_effect_median_percent",
  median(sp_range$pct_effect_s, na.rm = TRUE),
  "Median species-level range-based percentage difference in the composite mapped-pressure index",
  "%",
  1
)

stats <- add_stat(
  stats,
  "species_range_effect_negative_share_percent",
  100 * mean(sp_range$delta_s < 0),
  "Share of species with lower estimated composite mapped pressure inside protected areas",
  "%",
  1
)

taxon_range <- taxon_equal %>%
  transmute(
    taxon = taxon,
    n_species = n_species,
    median_pct_effect = pct_median
  )

for (i in seq_len(nrow(taxon_range))) {
  tx <- taxon_range$taxon[i]

  stats <- add_stat(
    stats,
    paste0("species_range_effect_", tx, "_n"),
    taxon_range$n_species[i],
    paste0("Threatened ", tx, " species with estimable range-based pressure contrasts"),
    "species",
    0
  )

  stats <- add_stat(
    stats,
    paste0("species_range_effect_", tx, "_median_percent"),
    taxon_range$median_pct_effect[i],
    paste0("Median ", tx, " range-based percentage difference in the composite mapped-pressure index"),
    "%",
    1
  )

  taxon_delta <- sp_range$delta_s[sp_range$taxon == tx]
  if (length(taxon_delta) != taxon_range$n_species[i]) {
    stop(
      "Canonical species rows for taxon ", tx,
      " do not match the centralized taxon count."
    )
  }

  stats <- add_stat(
    stats,
    paste0("species_range_effect_", tx, "_negative_share_percent"),
    100 * mean(taxon_delta < 0),
    paste0("Share of ", tx, " species with lower estimated composite mapped pressure inside protected areas"),
    "%",
    1
  )
}

#################################
### Hotspot heterogeneity #######
#################################

require_columns(
  hotspot_est,
  c("hotspot_id", "hotspot_name", "coef", "pval", "control_mean", "status"),
  "Hotspot estimates"
)

if (!all(hotspot_est$status == "estimated")) {
  stop("Hotspot estimates contain non-estimated rows.")
}

if (!any(hotspot_est$hotspot_id == 0)) {
  stop("Hotspot estimates are missing the non-hotspot reference group.")
}

hotspot_only <- hotspot_est %>%
  filter(
    hotspot_id != 0
  ) %>%
  mutate(
    pct_effect = 100 * coef / control_mean
  )

if (nrow(hotspot_only) < 1L) {
  stop("Hotspot estimates contain no named hotspot rows.")
}

if (any(!is.finite(hotspot_only$pct_effect))) {
  stop("Hotspot percentage effects contain nonfinite values.")
}

stats <- add_stat(
  stats,
  "hotspot_estimated_n",
  nrow(hotspot_only),
  "Biodiversity hotspots with estimable composite mapped-pressure contrasts",
  "hotspots",
  0
)

stats <- add_stat(
  stats,
  "hotspot_negative_n",
  sum(hotspot_only$coef < 0),
  "Biodiversity hotspots with lower estimated composite mapped pressure inside protected areas",
  "hotspots",
  0
)

stats <- add_stat(
  stats,
  "hotspot_significant_reduction_n",
  sum(hotspot_only$coef < 0 & hotspot_only$pval < 0.05),
  "Biodiversity hotspots with significant reductions at the 5 percent level",
  "hotspots",
  0
)

stats <- add_stat(
  stats,
  "hotspot_significant_reduction_p10_n",
  sum(hotspot_only$coef < 0 & hotspot_only$pval < 0.10),
  "Biodiversity hotspots with significant reductions at the 10 percent level",
  "hotspots",
  0
)

stats <- add_stat(
  stats,
  "hotspot_significant_increase_n",
  sum(hotspot_only$coef > 0 & hotspot_only$pval < 0.05),
  "Biodiversity hotspots with significant increases at the 5 percent level",
  "hotspots",
  0
)

stats <- add_stat(
  stats,
  "hotspot_median_percent",
  median(hotspot_only$pct_effect),
  "Median biodiversity-hotspot percentage difference in the composite mapped-pressure index",
  "%",
  1
)

hotspot_named <- c(
  sundaland = "Sundaland",
  wallacea = "Wallacea",
  indoburma = "Indo-Burma",
  cerrado = "Cerrado"
)

for (nm in names(hotspot_named)) {
  hotspot_row <- hotspot_only %>%
    filter(
      hotspot_name == hotspot_named[[nm]]
    )

  if (nrow(hotspot_row) != 1L) {
    stop("Expected exactly one hotspot row for ", hotspot_named[[nm]], ".")
  }

  stats <- add_stat(
    stats,
    paste0("hotspot_", nm, "_percent"),
    hotspot_row$pct_effect,
    paste0(
      hotspot_named[[nm]],
      " percentage difference in the composite mapped-pressure index"
    ),
    "%",
    1
  )
}

stats_df <- bind_rows(stats)
if (anyDuplicated(stats_df$key)) {
  stop("Headline-statistic keys must be unique.")
}

stats_csv <- file.path(paths$results_dir, "pathreat.analysis.statistics.csv")
write.csv(stats_df, stats_csv, row.names = FALSE)

stats_txt <- file.path(paths$results_dir, "pathreat.analysis.statistics.txt")

get_stat <- function(key) {
  index <- which(stats_df$key == key)
  if (length(index) != 1L) {
    stop("Headline-statistic key does not occur exactly once: ", key)
  }
  stats_df$value[index]
}

global_direction <- if (get_stat("global_composite_coef") < 0) {
  "lower"
} else {
  "higher"
}
any_threat_direction <- if (get_stat("global_any_threat_percentage_points") < 0) {
  "lower"
} else {
  "higher"
}

indonesia_balance_text <- if (post_balance$indonesia_in_stage) {
  "passes and is retained"
} else {
  "fails and is excluded"
}
indonesia_effect_text <- if (post_balance$indonesia_in_stage) {
  indonesia_effect_class <- if (indonesia_country_pval >= 0.05) {
    "not statistically significant at p < 0.05"
  } else if (indonesia_country_coef < 0) {
    "a statistically significant reduction"
  } else {
    "a statistically significant increase"
  }

  sprintf(
    paste0(
      "Indonesia composite country effect: ATT = %s, p = %s, ",
      "percentage difference = %s%% (%s)."
    ),
    fmt_num(indonesia_country_coef, 3),
    fmt_num(indonesia_country_pval, 3),
    fmt_num(indonesia_country_percent, 1),
    indonesia_effect_class
  )
} else {
  paste(
    "Indonesia composite country effect: absent because Indonesia is",
    "excluded by the unchanged country-balance rule."
  )
}

taxon_count_text <- paste(
  paste0(
    effect_taxa,
    " = ",
    fmt_int(observed_taxon_counts)
  ),
  collapse = "; "
)

txt <- c(
  "Headline statistics for manuscript text",
  "========================================",
  "",
  "Matched-sample audit",
  "--------------------",
  sprintf(
    paste0(
      "Pre-balance: %s observations; %s pairs; %s countries; %s PAs; ",
      "%s treated 1-km pixels; %s treated 25-km cells."
    ),
    fmt_int(get_stat("matched_pre_balance_observations")),
    fmt_int(get_stat("matched_pre_balance_pairs")),
    fmt_int(get_stat("matched_pre_balance_countries")),
    fmt_int(get_stat("matched_pre_balance_protected_areas")),
    fmt_int(get_stat("matched_pre_balance_treated_pixels")),
    fmt_int(get_stat("matched_pre_balance_cells_25km"))
  ),
  sprintf(
    paste0(
      "Post-balance: %s observations; %s pairs; %s countries; %s PAs; ",
      "%s treated 1-km pixels; %s treated 25-km cells."
    ),
    fmt_int(get_stat("matched_post_balance_observations")),
    fmt_int(get_stat("matched_post_balance_pairs")),
    fmt_int(get_stat("matched_post_balance_countries")),
    fmt_int(get_stat("matched_post_balance_protected_areas")),
    fmt_int(get_stat("matched_post_balance_treated_pixels")),
    fmt_int(get_stat("matched_post_balance_cells_25km"))
  ),
  sprintf(
    paste0(
      "Indonesia pre-balance: %s observations; %s pairs; %s PAs; ",
      "%s treated 1-km pixels; %s treated 25-km cells."
    ),
    fmt_int(get_stat("indonesia_pre_balance_observations")),
    fmt_int(get_stat("indonesia_pre_balance_pairs")),
    fmt_int(get_stat("indonesia_pre_balance_protected_areas")),
    fmt_int(get_stat("indonesia_pre_balance_treated_pixels")),
    fmt_int(get_stat("indonesia_pre_balance_cells_25km"))
  ),
  sprintf(
    paste0(
      "Indonesia post-balance: %s observations; %s pairs; %s PAs; ",
      "%s treated 1-km pixels; %s treated 25-km cells."
    ),
    fmt_int(get_stat("indonesia_post_balance_observations")),
    fmt_int(get_stat("indonesia_post_balance_pairs")),
    fmt_int(get_stat("indonesia_post_balance_protected_areas")),
    fmt_int(get_stat("indonesia_post_balance_treated_pixels")),
    fmt_int(get_stat("indonesia_post_balance_cells_25km"))
  ),
  sprintf(
    "Indonesia balance: average |SMD| = %s against threshold %s; %s.",
    fmt_num(get_stat("indonesia_avg_abs_smd"), 3),
    fmt_num(get_stat("indonesia_balance_threshold"), 2),
    indonesia_balance_text
  ),
  "",
  "Global and country effects",
  "--------------------------",
  sprintf(
    "Composite mapped-pressure index: %s%% %s inside PAs (ATT = %s, p = %s).",
    fmt_num(abs(get_stat("global_composite_percent_change")), 1),
    global_direction,
    fmt_num(get_stat("global_composite_coef"), 3),
    fmt_num(get_stat("global_composite_p_value"), 3)
  ),
  sprintf(
    "Any-threat indicator: %s percentage points %s inside PAs (p = %s).",
    fmt_num(abs(get_stat("global_any_threat_percentage_points")), 1),
    any_threat_direction,
    fmt_num(get_stat("global_any_threat_p_value"), 3)
  ),
  sprintf(
    "Individual indicators lower inside PAs: %d at p < 0.05; %d at p < 0.10.",
    as.integer(get_stat("individual_threats_lower_p05")),
    as.integer(get_stat("individual_threats_lower_p10"))
  ),
  sprintf(
    "Composite estimation sample: %s matched 1-km protected-control pairs.",
    fmt_int(get_stat("global_composite_matched_pairs"))
  ),
  sprintf(
    paste0(
      "Composite country results: %s estimated; %s significant reductions; ",
      "%s significant increases; %s null; %s not estimated."
    ),
    fmt_int(get_stat("country_composite_estimated_n")),
    fmt_int(get_stat("country_composite_significant_reduction_n")),
    fmt_int(get_stat("country_composite_significant_increase_n")),
    fmt_int(get_stat("country_composite_null_n")),
    fmt_int(get_stat("country_composite_not_estimated_n"))
  ),
  indonesia_effect_text,
  "",
  "PA and stochastic-frontier results",
  "----------------------------------",
  sprintf(
    "PA regression results: %s records, %s fitted coefficients, %s with inference, %s without inference.",
    fmt_int(get_stat("pa_results_n")), fmt_int(get_stat("pa_coefficients_n")),
    fmt_int(get_stat("pa_inference_n")), fmt_int(get_stat("pa_no_inference_n"))
  ),
  sprintf(
    "SFA sample: %s PAs across %s countries; summed PA area = %s km2.",
    fmt_int(get_stat("sfa_pa_n")),
    fmt_int(get_stat("sfa_country_n")),
    fmt_num(get_stat("sfa_total_pa_area_km2"), 1)
  ),
  sprintf(
    paste0(
      "SFA area-weighted totals: achieved reduction = %s threat-index km2 ",
      "(%s%%); underperformance = %s threat-index km2 (%s%%); ",
      "total potential reduction = %s threat-index km2."
    ),
    fmt_num(get_stat("sfa_aggregate_achieved_threat_reduction"), 1),
    fmt_num(get_stat("sfa_share_achieved_percent"), 1),
    fmt_num(get_stat("sfa_aggregate_underperformance"), 1),
    fmt_num(get_stat("sfa_share_underperformance_percent"), 1),
    fmt_num(get_stat("sfa_aggregate_total_potential_reduction"), 1)
  ),
  sprintf(
    paste0(
      "SFA diagnostics: median underperformance = %s; median net residual = %s; ",
      "%s%% of PAs have a negative net residual."
    ),
    fmt_num(get_stat("sfa_median_underperformance"), 3),
    fmt_num(get_stat("sfa_median_net_residual"), 3),
    fmt_num(get_stat("sfa_negative_residual_percent"), 1)
  ),
  sprintf(
    "Country achieved-share ratio: median %s%%; IQR %s--%s%%.",
    fmt_num(get_stat("sfa_country_achieved_share_median_percent"), 1),
    fmt_num(get_stat("sfa_country_achieved_share_p25_percent"), 1),
    fmt_num(get_stat("sfa_country_achieved_share_p75_percent"), 1)
  ),
  "",
  "Species and taxon results",
  "-------------------------",
  sprintf(
    "Species coverage analysis: %s threatened terrestrial vertebrate species.",
    fmt_int(get_stat("species_coverage_n"))
  ),
  sprintf(
    "No mapped PA overlap: %s species (%s%%).",
    fmt_int(get_stat("species_no_pa_overlap_n")),
    fmt_num(get_stat("species_no_pa_overlap_percent"), 1)
  ),
  sprintf(
    "Mean mapped range shares: %s%% protected and %s%% unprotected.",
    fmt_num(get_stat("species_mean_range_protected_percent"), 1),
    fmt_num(get_stat("species_mean_range_unprotected_percent"), 1)
  ),
  "",
  sprintf(
    "Range-based species effect sample: %s species; median composite pressure difference = %s%%.",
    fmt_int(get_stat("species_range_effect_estimable_n")),
    fmt_num(get_stat("species_range_effect_median_percent"), 1)
  ),
  paste0("Eligible species by taxon: ", taxon_count_text, ".")
)

writeLines(txt, stats_txt)

cat(sprintf("Saved: %s\n", stats_csv))
cat(sprintf("Saved: %s\n", stats_txt))
cat("\n")
cat(paste(txt, collapse = "\n"))
cat("\n")
