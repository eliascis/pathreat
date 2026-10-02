##################################################
### pathreat.analysis.byPA.est.R ##################
### Individual-PA regressions and classification ##
##################################################
# Fits the composite outcome on the unchanged complete-pair matched sample.
# Outputs: canonical PA estimates, PA-country membership, and migration audits.

library(dplyr)
source("code/pathreat.analysis.config.R")

##################################################
### Complete-pair sample and regression helpers ###
##################################################

build_pa_sample <- function(data) {
  pair_index <- build_matched_pair_index(data)
  sample_mask <- matched_pair_sample_mask(data, "threat_composite", pair_index)
  keep <- sample_mask[pair_index$treated_rows]
  tr <- pair_index$treated_rows[keep]
  cr <- pair_index$control_rows[keep]
  if (!length(tr)) stop("No eligible PA pairs")
  if (anyNA(data[tr, c("wdpaid", "country_rast", "biome_raster")]) ||
      anyNA(data[cr, c("country_rast", "biome_raster")]) ||
      any(data$country_rast[tr] != data$country_rast[cr]) ||
      any(data$biome_raster[tr] != data$biome_raster[cr])) {
    stop("PA pairs must have complete, identical stored country/biome labels")
  }
  difference <- data$threat_composite[tr] - data$threat_composite[cr]
  if (any(!is.finite(data$delta[tr])) ||
      any(abs(difference - data$delta[tr]) > 1e-8)) {
    stop("Stored delta disagrees with independently calculated pair differences")
  }
  sample <- data[c(tr, cr),
                 c("treat", "country_rast", "biome_raster", "threat_composite"),
                 drop = FALSE]
  # A reused control remains a separate row for every pair. Its own wdpaid
  # never determines membership in a PA regression.
  sample$wdpaid <- rep(data$wdpaid[tr], 2L)
  pairs <- data.frame(
    wdpaid = data$wdpaid[tr], country_rast = data$country_rast[tr],
    delta = difference, tc_control = data$threat_composite[cr]
  )
  list(sample = sample, pairs = pairs)
}

# Keep the public PA helper name for callers and regression tests.
# fixest rejects constant outcomes. Base OLS also retains the saturated
# single-pair coefficient without attempting residual-based inference.
# Estimate first, then request the species estimator's covariance with
# unchanged fixest defaults (including ssc), so failed inference cannot
# erase a valid coefficient. Estimation errors remain fatal.
fit_pa_regression <- function(dd) {
  fit_pair_subgroup(dd)
}

classify_pa_effectiveness <- function(d_pa) {
  valid <- d_pa$status == "estimated" & is.finite(d_pa$p_value)
  sig_neg <- d_pa$coef[valid & d_pa$p_value < 0.05 & d_pa$coef < 0]
  cutpoints <- if (length(sig_neg)) {
    unname(quantile(sig_neg, c(1/3, 2/3)))
  } else c(NA_real_, NA_real_)
  d_pa$effectiveness <- case_when(
    !valid ~ "no_inference",
    d_pa$p_value >= 0.05 ~ "not_significant",
    d_pa$coef >= 0 ~ "harmful",
    d_pa$coef <= cutpoints[1] ~ "high",
    d_pa$coef <= cutpoints[2] ~ "medium",
    TRUE ~ "low"
  )
  attr(d_pa, "tercile_cutpoints") <- cutpoints
  d_pa
}

##################################################
### Fit and audit one row per global WDPA ID ######
##################################################

run_pa_estimation <- function() {
  outfile_pa <- "data/store/pathreat.analysis.byPA.est.Rds"
  outfile_membership <- "data/store/pathreat.analysis.byPA-country-membership.Rds"
  comparison_file <- "results/pathreat.byPA.inference-comparison.csv"
  transition_file <- "results/pathreat.byPA.classification-transitions.csv"
  old <- if (file.exists(outfile_pa)) readRDS(outfile_pa) else NULL
  old_membership <- if (file.exists(outfile_membership)) readRDS(outfile_membership) else NULL
  columns <- c("matched_pair_id", "treat", "wdpaid", "delta", "country_rast",
               "biome_raster", "threat_composite")
  cat("Loading complete-pair PA sample...\n")
  data <- fst::read_fst(paths$data_matched, columns = columns)
  inputs <- build_pa_sample(data)
  rm(data)
  d_pa_country <- tibble::as_tibble(inputs$pairs) %>%
    count(wdpaid, country_rast, name = "n_pixels_country") %>%
    arrange(wdpaid, country_rast)
  d_pa <- inputs$pairs %>%
    group_by(wdpaid) %>%
    summarise(
      mean_pair_difference = mean(delta), sd_delta = sd(delta),
      mean_control = mean(tc_control), n_pixels = n(),
      n_countries = n_distinct(country_rast), country_rast = min(country_rast),
      .groups = "drop"
    )
  # Precompute row groups once, avoiding a complete-data scan for every PA.
  groups <- split(seq_len(nrow(inputs$sample)), inputs$sample$wdpaid)
  groups <- groups[as.character(d_pa$wdpaid)]
  cat(sprintf("Fitting %d PA regressions (%d complete pairs)...\n",
              nrow(d_pa), sum(d_pa$n_pixels)))
  fits <- lapply(seq_along(groups), function(k) {
    if (k %% 2000L == 0L) cat(sprintf("  PA %d / %d\n", k, length(groups)))
    tryCatch(fit_pa_regression(inputs$sample[groups[[k]], , drop = FALSE]),
             error = function(err) stop("PA ", d_pa$wdpaid[k], ": ", conditionMessage(err)))
  })
  fitted <- do.call(rbind, fits)
  d_pa <- bind_cols(d_pa, fitted) %>%
    mutate(mean_delta = coef)
  d_pa <- classify_pa_effectiveness(d_pa)
  if (anyDuplicated(d_pa$wdpaid) || any(!is.finite(d_pa$coef)) ||
      any(d_pa$n_obs != 2L * d_pa$n_pixels) ||
      any(abs(d_pa$coef - d_pa$mean_pair_difference) > 1e-8)) {
    stop("PA regression coefficient/complete-pair audit failed")
  }
  membership_check <- d_pa_country %>%
    group_by(wdpaid) %>%
    summarise(n_pixels = sum(n_pixels_country), n_countries = n(), .groups = "drop")
  stopifnot(identical(d_pa$wdpaid, membership_check$wdpaid),
            all(d_pa$n_pixels == membership_check$n_pixels),
            all(d_pa$n_countries == membership_check$n_countries))
  if (!is.null(old_membership) &&
      !isTRUE(all.equal(as.data.frame(d_pa_country), as.data.frame(old_membership),
                       check.attributes = FALSE))) {
    stop("PA-country membership changed: outputs not published")
  }
  if (!is.null(old)) {
    if (!setequal(old$wdpaid, d_pa$wdpaid)) stop("PA identifier set changed")
    old <- old[match(d_pa$wdpaid, old$wdpaid), ]
    if (any(old$n_pixels != d_pa$n_pixels) ||
        any(abs(old$mean_delta - d_pa$coef) > 1e-8)) {
      stop("PA pair counts or coefficients changed beyond tolerance")
    }
  }
  old_field <- function(name, default = NA_real_) {
    if (!is.null(old) && name %in% names(old)) old[[name]] else rep(default, nrow(d_pa))
  }
  comparison <- data.frame(
    wdpaid = d_pa$wdpaid, n_pixels = d_pa$n_pixels,
    old_coef = old_field("mean_delta"), new_coef = d_pa$coef,
    mean_pair_difference = d_pa$mean_pair_difference,
    old_se = old_field("se_delta"), new_se = d_pa$se_delta,
    old_t_stat = old_field("t_stat"), new_t_stat = d_pa$t_stat,
    old_p_value = old_field("p_value"), new_p_value = d_pa$p_value,
    old_effectiveness = old_field("effectiveness", NA_character_),
    new_effectiveness = d_pa$effectiveness,
    old_status = old_field("status", "legacy_paired_difference"), new_status = d_pa$status,
    old_significant = is.finite(old_field("p_value")) & old_field("p_value") < 0.05,
    new_significant = d_pa$status == "estimated" & is.finite(d_pa$p_value) & d_pa$p_value < 0.05,
    fe_spec = d_pa$fe_spec, se_type = d_pa$se_type, fit_engine = d_pa$fit_engine,
    ci_low = d_pa$ci_low, ci_high = d_pa$ci_high, inference_reason = d_pa$inference_reason
  )
  # Preserve the first legacy-to-regression audit across routine reruns.
  if (!is.null(old) && "coef" %in% names(old) && file.exists(comparison_file)) {
    first <- read.csv(comparison_file)
    if (setequal(first$wdpaid, comparison$wdpaid)) {
      ix <- match(comparison$wdpaid, first$wdpaid)
      old_columns <- grep("^old_", names(comparison), value = TRUE)
      comparison[old_columns] <- first[ix, old_columns]
    }
  }
  transitions <- comparison %>%
    count(old_effectiveness, new_effectiveness, name = "n_pas")
  attr(d_pa, "estimation") <- list(
    specification = "unweighted OLS; varying additive country/biome FE",
    fixest_version = as.character(utils::packageVersion("fixest")),
    ssc_defaults = fixest::ssc(), ssc_overrides = fixest::getFixest_ssc(),
    confidence_level = 0.95,
    matched_input_md5 = unname(tools::md5sum(paths$data_matched)),
    coefficient_tolerance = 1e-8
  )
  cat("\nEstimation statuses:\n")
  print(table(d_pa$status))
  cat("\nClassification transitions:\n")
  print(as.data.frame(transitions), row.names = FALSE)
  cat("\nSignificance transitions:\n")
  print(table(old = comparison$old_significant, new = comparison$new_significant))
  cutpoints <- attr(d_pa, "tercile_cutpoints")
  cat(if (all(is.na(cutpoints))) "No tercile thresholds defined.\n" else
      sprintf("Tercile thresholds: %.10f, %.10f\n", cutpoints[1], cutpoints[2]))
  cat(sprintf("Maximum coefficient / pair-difference gap: %.3g\n",
              max(abs(d_pa$coef - d_pa$mean_pair_difference))))

  # Stage and round-trip every artifact before atomic, same-directory renames.
  destinations <- c(outfile_pa, outfile_membership, comparison_file, transition_file)
  objects <- list(d_pa, d_pa_country, comparison, transitions)
  staged <- vapply(destinations, function(path) tempfile(".byPA-", dirname(path)), character(1))
  on.exit(unlink(staged), add = TRUE)
  for (k in seq_along(destinations)) {
    if (k <= 2L) {
      saveRDS(objects[[k]], staged[k])
      stopifnot(identical(readRDS(staged[k]), objects[[k]]))
    } else {
      write.csv(objects[[k]], staged[k], row.names = FALSE)
      stopifnot(nrow(read.csv(staged[k])) == nrow(objects[[k]]))
    }
  }
  for (k in seq_along(destinations)) {
    if (!file.rename(staged[k], destinations[k])) stop("Atomic publication failed: ", destinations[k])
  }
  invisible(d_pa)
}

# Source helpers without data loading for focused regression checks.
if (!isTRUE(getOption("pathreat.byPA.functions_only", FALSE))) {
  d_pa <- run_pa_estimation()
}

##
# head(d_pa)
