# Run from the project root: Rscript code/tests/pathreat.analysis.byPA.R
options(pathreat.byPA.functions_only = TRUE)
source("code/pathreat.analysis.byPA.est.R")

make_pairs <- function(country = rep(1, 8), biome = rep(1, 8)) {
  n <- length(country)
  control <- seq(0.1, 0.8, length.out = n)
  protected <- control + rep(c(-0.03, 0.01, -0.02, -0.04), length.out = n)
  data.frame(
    matched_pair_id = rep(seq_len(n), 2), treat = rep(c(1, 0), each = n),
    wdpaid = c(rep(42, n), rep(999, n)),
    country_rast = rep(country, 2), biome_raster = rep(biome, 2),
    threat_composite = c(protected, control), delta = rep(protected - control, 2),
    # Repeated pixel IDs represent reused donors, without collapsing their rows.
    pixel_id = c(seq_len(n), 100 + seq_len(n))
  )
}
check_fit <- function(country, biome, fe_spec, se_type) {
  data <- make_pairs(country, biome)
  sample <- build_pa_sample(data)$sample
  result <- fit_pa_regression(sample)
  stopifnot(result$fe_spec == fe_spec, result$se_type == se_type,
            result$n_obs == nrow(data), result$status == "estimated",
            abs(result$coef - mean(data$delta[data$treat == 1])) < 1e-8)
  fe <- c(if (length(unique(country)) > 1) "country_rast",
          if (length(unique(biome)) > 1) "biome_raster")
  fml <- as.formula(paste("threat_composite ~ treat",
                         if (length(fe)) paste("|", paste(fe, collapse = " + "))))
  direct <- if (se_type == "cluster_country") {
    fixest::feols(fml, sample, cluster = ~country_rast)
  } else fixest::feols(fml, sample, vcov = "hetero")
  stopifnot(isTRUE(all.equal(as.numeric(result[1, c("coef", "se_delta", "t_stat", "p_value")]),
                            as.numeric(fixest::coeftable(direct)["treat", ]))),
            isTRUE(all.equal(as.numeric(result[1, c("ci_low", "ci_high")]),
                             as.numeric(confint(direct, "treat")))))
  invisible(result)
}
one <- rep(1, 8)
countries <- rep(1:2, each = 4)
biomes <- rep(1:2, 4)
check_fit(one, one, "none", "hetero")
check_fit(one, biomes, "biome", "hetero")
check_fit(countries, one, "country", "cluster_country")
check_fit(countries, biomes, "country+biome", "cluster_country")

original <- make_pairs(countries, biomes)
# The same donor, with identical outcome and stored FE labels, is reused in
# pairs 1 and 3. Its own PA identifier differs from the protected member's.
original$pixel_id[11] <- original$pixel_id[9]
original$threat_composite[11] <- original$threat_composite[9]
original$delta[c(3, 11)] <- original$threat_composite[3] - original$threat_composite[11]
set.seed(18)
shuffled <- original[sample(nrow(original)), ]
x <- build_pa_sample(original)
y <- build_pa_sample(shuffled)
stopifnot(all(x$sample$wdpaid == 42), nrow(x$sample) == 16,
          abs(fit_pa_regression(x$sample)$coef - fit_pa_regression(y$sample)$coef) < 1e-8)
# Two protected WDPA identifiers each propagate to their complete pair rows.
multiple <- original
multiple$wdpaid[c(1, 3)] <- 43
multi <- build_pa_sample(multiple)
stopifnot(sum(multi$sample$wdpaid == 43) == 4,
          sum(multi$sample$wdpaid == 42) == 12,
          all(multi$sample$wdpaid %in% c(42, 43)))
# Missing control or protected outcomes remove both pair members.
missing <- original
missing$threat_composite[c(1, 10)] <- NA_real_
z <- build_pa_sample(missing)
stopifnot(nrow(z$sample) == 12, nrow(z$pairs) == 6,
          all(z$sample$wdpaid == 42), all(is.finite(z$sample$threat_composite)))
single <- build_pa_sample(original[c(1, 9), ])$sample
s <- fit_pa_regression(single)
stopifnot(s$status == "single_pair", s$fit_engine == "lm", is.na(s$p_value),
          abs(s$coef + 0.03) < 1e-8)
constant <- x$sample
constant$threat_composite <- 0.2
c <- fit_pa_regression(constant)
stopifnot(c$status == "constant_outcome", c$fit_engine == "lm",
          abs(c$coef) < 1e-8, is.na(c$se_delta), is.na(c$ci_low))
# Nonconstant outcome entirely absorbed by FE: coefficient remains available.
degenerate <- x$sample
degenerate$threat_composite <- as.numeric(degenerate$country_rast)
u <- suppressWarnings(fit_pa_regression(degenerate))
stopifnot(u$status == "no_inference", is.finite(u$coef),
          is.na(u$p_value), !is.na(u$inference_reason))
# Structural sample errors must fail rather than silently dropping rows.
bad <- original
bad$biome_raster[9] <- 99
stopifnot(inherits(try(build_pa_sample(bad), silent = TRUE), "try-error"))

dropped <- x$sample
dropped$threat_composite[1] <- NA_real_
stopifnot(inherits(try(suppressWarnings(fit_pa_regression(dropped)), silent = TRUE), "try-error"))

classes <- data.frame(coef = c(-3, -2, -1, 1, -2, 0),
                      p_value = c(.01, .01, .01, .01, .5, NA),
                      status = c(rep("estimated", 5), "single_pair"))
cl <- classify_pa_effectiveness(classes)
stopifnot(identical(cl$effectiveness,
                    c("high", "medium", "low", "harmful", "not_significant", "no_inference")))
empty <- classify_pa_effectiveness(classes[4:6, ])
stopifnot(all(is.na(attr(empty, "tercile_cutpoints"))),
          !any(empty$effectiveness %in% c("high", "medium", "low")))
cat("All PA regression checks passed.\n")

# Standardized levels must follow protected observations and preserve the ATT.
for (sample in list(x$sample, single, constant, degenerate)) {
  fitted <- suppressWarnings(fit_pair_subgroup(sample))
  stopifnot(abs(fitted$adjusted_protected - fitted$adjusted_control - fitted$coef) < 1e-8,
            abs(fitted$adjusted_control - mean(sample$threat_composite[sample$treat == 0])) < 1e-8,
            abs(fitted$adjusted_protected - mean(sample$threat_composite[sample$treat == 1])) < 1e-8)
}
# A two-country cluster interval must use model degrees of freedom, rather
# than the normal critical value. Keep the same fitted covariance and test.
fitted <- fixest::feols(threat_composite ~ treat | country_rast + biome_raster,
                        data = x$sample, cluster = ~country_rast)
estimate <- regression_estimate(fitted)
stopifnot(isTRUE(all.equal(unname(as.numeric(estimate[1, c("ci_low", "ci_high")])),
                          unname(as.numeric(confint(fitted, "treat"))))),
          abs((estimate$ci_high - estimate$coef) / estimate$se - 1.96) > 1)
# Extraction also respects a caller's nonuniform weights and supplied vcov.
weighted <- fixest::feols(threat_composite ~ treat | country_rast + biome_raster,
                          data = x$sample, weights = seq_len(nrow(x$sample)),
                          vcov = "hetero")
weighted_estimate <- regression_estimate(weighted)
stopifnot(isTRUE(all.equal(weighted_estimate$coef, unname(coef(weighted)["treat"]))),
          isTRUE(all.equal(as.numeric(weighted_estimate[1, c("ci_low", "ci_high")]),
                           as.numeric(confint(weighted, "treat")))))
cat("All shared regression and standardized-level checks passed.\n")

# Explicit-indicator QR covariance must preserve the absorbed model's
# small-sample correction, including FE nested within country clusters.
qr_country <- rep(1:3, each = 12)
qr_biome <- rep(rep(1:3, each = 4), 3)
qr_cases <- list(
  list(rep(1, 36), rep(1, 36)),
  list(rep(1, 36), qr_biome),
  list(qr_country, rep(1, 36)),
  list(qr_country, qr_biome),
  list(qr_country, qr_country)
)
for (labels in qr_cases) {
  sample <- build_pa_sample(make_pairs(labels[[1]], labels[[2]]))$sample
  sample$threat_composite[sample$treat == 1] <-
    sample$threat_composite[sample$treat == 1] + seq(-0.02, 0.02, length.out = 36)
  n_countries <- length(unique(sample$country_rast))
  fe <- c(if (n_countries > 1) "country_rast",
          if (length(unique(sample$biome_raster)) > 1) "biome_raster")
  formula <- as.formula(paste("threat_composite ~ treat",
    if (length(fe)) paste("|", paste(fe, collapse = " + "))))
  direct <- fixest::feols(formula, data = sample)
  reference <- if (n_countries > 1) {
    summary(direct, cluster = ~country_rast)
  } else summary(direct, vcov = "hetero")
  qr_fit <- lm(reformulate(c("treat", if (length(fe)) paste0("factor(", fe, ")")),
                           response = "threat_composite"), data = sample)
  qr_values <- regression_qr_inference(qr_fit, direct, sample, n_countries)
  direct_values <- regression_estimate(reference)
  stopifnot(max(abs(as.numeric(qr_values) - as.numeric(direct_values))) < 1e-8)
}
cat("All QR covariance and nested fixed-effect checks passed.\n")

# Force an absorption failure on an unbalanced FE layout. The normal public
# fitter must recover by QR, without depending on a particular species fixture.
local({
  previous_settings <- fixest::getFixest_estimation()
  on.exit(do.call(fixest::setFixest_estimation,
                  c(list(reset = TRUE), previous_settings)))
  sample <- build_pa_sample(make_pairs(
    c(1, 1, 1, 2, 2, 3, 3, 3), c(1, 1, 2, 2, 3, 1, 2, 3)
  ))$sample
  reference <- fixest::feols(
    threat_composite ~ treat | country_rast + biome_raster,
    data = sample, cluster = ~country_rast, fixef.tol = 1e-10
  )
  fixest::setFixest_estimation(fixef.iter = 1)
  recovered <- suppressWarnings(fit_pair_subgroup(sample))
  expected <- regression_estimate(reference)
  stopifnot(recovered$fit_engine == "lm_qr_fallback",
            recovered$status == "estimated", !is.na(recovered$inference_reason),
            max(abs(as.numeric(recovered[1, c("coef", "se_delta", "t_stat", "p_value",
                                             "ci_low", "ci_high")]) -
                    as.numeric(expected[1, c("coef", "se", "t_stat", "pval",
                                            "ci_low", "ci_high")]))) < 1e-8,
            abs(recovered$adjusted_control -
                  mean(sample$threat_composite[sample$treat == 0])) < 1e-10)
})
cat("Forced absorption failure recovered with the correct regression inference.\n")
