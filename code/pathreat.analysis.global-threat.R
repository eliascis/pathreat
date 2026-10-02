##############################################
### pathreat.analysis.global-threat.R ########
### Global threat_composite ATT ###############
### Biodiversity-weighted variants ############
##############################################
#
# Re-estimates the main global threat_composite ATT from
# pathreat.analysis.global.est.R, then reweights pixels by
# biodiversity importance using pixel-level species data.
#
# Weighting design:
# - Each pixel is weighted by its own species importance score,
#   computed in species-count.R and stored in merge.matched.fst.
# - Red List weighting uses pixel-level threatened-species richness.
# - TBL, ED, and EDGE weighting use pixel-level sums of the
#   corresponding species scores across all present species.
# - Species with NA phylogenetic scores (unmatched to EDGE dataset)
#   contribute to species counts but not to score sums.
#
# Inputs:
#   data/store/pathreat.data.merge.matched.fst
#
# Outputs:
#   results/pathreat.global-threat.est.Rds
#   results/pathreat.global-threat.est.html
#   results/pathreat.est.newest.html
#

suppressPackageStartupMessages({
  library(dplyr)
  library(fixest)
  library(fst)
  library(texreg)
})

source("code/pathreat.analysis.config.R")

main_dep <- "threat_composite"
out_rds <- file.path(paths$results_dir, "pathreat.global-threat.est.Rds")
out_html <- file.path(paths$results_dir, "pathreat.global-threat.est.html")
out_html_latest <- file.path(paths$results_dir, "pathreat.est.newest.html")


run_global_model <- function(data, weight_col = NULL) {
  fml <- threat_composite ~ treat | country_rast + biome_raster

  if (is.null(weight_col)) {
    fit <- feols(fml, data = data, cluster = ~country_rast)
  } else {
    weights <- data[[weight_col]]
    if (any(!is.finite(weights)) || any(weights <= 0)) {
      stop(sprintf("%s must be finite and strictly positive", weight_col))
    }
    fit <- feols(
      fml,
      data = data,
      weights = as.formula(paste0("~", weight_col)),
      cluster = ~country_rast
    )
  }
  if (fit$nobs != nrow(data)) {
    stop("Weighted-global model changed the canonical matched-pair sample")
  }
  fit
}

extract_fit_row <- function(model_name, weight_definition, fit, data,
                            weight_col = NULL) {
  ct <- coeftable(fit)["treat", ]
  estimate <- regression_estimate(fit)

  if (is.null(weight_col)) {
    w <- rep(1, nrow(data))
  } else {
    w <- data[[weight_col]]
    if (any(!is.finite(w)) || any(w <= 0)) {
      stop(sprintf("%s contains missing, non-finite, or non-positive weights", weight_col))
    }
  }

  treat_idx <- data$treat == 1
  ctrl_idx <- data$treat == 0
  if (sum(treat_idx) != sum(ctrl_idx) || nrow(data) != 2L * sum(treat_idx)) {
    stop("Weighted-global input sample is not pair-complete")
  }
  control_mean <- if (is.null(weight_col)) {
    mean(data$threat_composite[ctrl_idx], na.rm = TRUE)
  } else {
    weighted.mean(data$threat_composite[ctrl_idx], w[ctrl_idx], na.rm = TRUE)
  }

  data.frame(
    model = model_name,
    weight_definition = weight_definition,
    coef = ct["Estimate"],
    se = ct["Std. Error"],
    t_stat = ct["t value"],
    pval = ct["Pr(>|t|)"],
    ci_low = estimate$ci_low,
    ci_high = estimate$ci_high,
    n_obs = fit$nobs,
    n_pairs = sum(treat_idx, na.rm = TRUE),
    control_mean = control_mean,
    pct_of_control = ifelse(
      is.na(control_mean) || control_mean == 0,
      NA_real_,
      ct["Estimate"] / control_mean * 100
    ),
    mean_weight_treat = mean(w[treat_idx], na.rm = TRUE),
    median_weight_treat = stats::median(w[treat_idx], na.rm = TRUE),
    min_weight_treat = min(w[treat_idx], na.rm = TRUE),
    max_weight_treat = max(w[treat_idx], na.rm = TRUE),
    share_positive_weight_treat = mean(w[treat_idx] > 0, na.rm = TRUE),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


############################
### load matched data ######
############################
{
  cat("=== Loading matched data ===\n")

  target_filter <- get_outcome_target_filter(main_dep)
  target_filter_vars <- if (is.null(target_filter)) {
    character()
  } else {
    all.vars(parse(text = target_filter)[[1]])
  }
  needed_cols <- unique(c(
    "matched_pair_id",
    "treat",
    "country_rast",
    "biome_raster",
    "threat_composite",
    "n_redlist_species",
    "tbl_sum",
    "ed_sum",
    "edge_sum",
    target_filter_vars
  ))

  d <- read_fst(paths$data_matched, columns = needed_cols)
  cat(sprintf("  Loaded %s rows\n", format(nrow(d), big.mark = ",")))

  pair_index <- build_matched_pair_index(d)
  sample_mask <- matched_pair_sample_mask(d, main_dep, pair_index)
  d <- d[sample_mask, , drop = FALSE]
  n_treated <- sum(d$treat == 1)
  n_control <- sum(d$treat == 0)
  if (n_treated != n_control || nrow(d) != 2L * n_treated) {
    stop("Weighted-global canonical sample is not pair-complete")
  }

  cat(sprintf("  Estimation sample: %s rows, %s matched pairs\n",
              format(nrow(d), big.mark = ","),
              format(n_treated, big.mark = ",")))

  cat("\n  Pixel-level weight summaries (treated pixels):\n")
  treat_idx <- d$treat == 1
  for (wc in c("n_redlist_species", "tbl_sum", "ed_sum", "edge_sum")) {
    x <- d[[wc]][treat_idx]
    cat(sprintf("    %-20s mean=%.2f  median=%.2f  min=%.2f  max=%.2f  pct_positive=%.1f%%\n",
                wc,
                mean(x, na.rm = TRUE),
                stats::median(x, na.rm = TRUE),
                min(x, na.rm = TRUE),
                max(x, na.rm = TRUE),
                100 * mean(x > 0, na.rm = TRUE)))
  }
}


############################
### estimate models ########
############################
{
  cat("\n=== Estimating global ATT models ===\n")

  model_specs <- list(
    list(
      name = "Equal weights",
      weight_col = NULL,
      definition = "Replicates the main threat_composite ATT from global.est.R"
    ),
    list(
      name = "Threatened spp.",
      weight_col = "n_redlist_species",
      definition = "Pixel-level threatened-species richness"
    ),
    list(
      name = "TBL",
      weight_col = "tbl_sum",
      definition = "Pixel-level sum of species terminal branch length scores"
    ),
    list(
      name = "ED",
      weight_col = "ed_sum",
      definition = "Pixel-level sum of species evolutionary distinctiveness scores"
    ),
    list(
      name = "EDGE",
      weight_col = "edge_sum",
      definition = "Pixel-level sum of species EDGE2 scores"
    )
  )

  fit_list <- list()
  estimate_rows <- vector("list", length(model_specs))

  for (i in seq_along(model_specs)) {
    spec <- model_specs[[i]]
    cat(sprintf("  %-20s ", spec$name))
    fit <- run_global_model(d, spec$weight_col)
    fit_list[[i]] <- fit
    estimate_rows[[i]] <- extract_fit_row(
      model_name = spec$name,
      weight_definition = spec$definition,
      fit = fit,
      data = d,
      weight_col = spec$weight_col
    )
    cat(sprintf("coef = %.6f, se = %.6f, p = %.4f\n",
                estimate_rows[[i]]$coef,
                estimate_rows[[i]]$se,
                estimate_rows[[i]]$pval))
  }

  names(fit_list) <- vapply(model_specs, `[[`, character(1), "name")
  estimate_results <- bind_rows(estimate_rows)
}


#################################
### compare to stored global ####
#################################
{
  baseline_check <- NULL

  if (file.exists(paths$est_global)) {
    stored_global <- readRDS(paths$est_global)
    stored_row <- stored_global[
      stored_global$variable == main_dep &
        stored_global$estimate_type == "unscaled",
      ,
      drop = FALSE
    ]

    if (nrow(stored_row) == 1) {
      baseline_row <- estimate_results[estimate_results$model == "Equal weights", , drop = FALSE]
      baseline_check <- data.frame(
        coef_current = baseline_row$coef,
        coef_stored = stored_row$coef,
        coef_diff = baseline_row$coef - stored_row$coef,
        se_current = baseline_row$se,
        se_stored = stored_row$se,
        se_diff = baseline_row$se - stored_row$se,
        row.names = NULL
      )

      cat(sprintf("\nStored global.est comparison: coef diff = %.12f, se diff = %.12f\n",
                  baseline_check$coef_diff,
                  baseline_check$se_diff))
    }
  }
}


############################
### save outputs ###########
############################
{
  html <- htmlreg(
    fit_list,
    omit.coef = "(Intercept)",
    custom.model.names = estimate_results$model,
    custom.gof.rows = list(
      "Control Mean" = sprintf("%.3f", estimate_results$control_mean),
      "Effect (% control mean)" = sprintf("%.1f", estimate_results$pct_of_control),
      "Mean Treat Weight" = sprintf("%.2f", estimate_results$mean_weight_treat)
    ),
    stars = c(0.01, 0.05, 0.1),
    digits = 3,
    table = FALSE
  )

  cat(
    file = out_html,
    c("<b>Global threat_composite ATT with biodiversity-weighted variants</b>", html),
    append = FALSE
  )
  cat(
    file = out_html_latest,
    c("<b>Global threat_composite ATT with biodiversity-weighted variants</b>", html),
    append = FALSE
  )

  out_obj <- list(
    estimates = estimate_results,
    baseline_check = baseline_check
  )
  saveRDS(out_obj, out_rds)

  cat(sprintf("\nSaved HTML: %s\n", out_html))
  cat(sprintf("Saved HTML: %s\n", out_html_latest))
  cat(sprintf("Saved RDS:  %s\n", out_rds))
}

## Figure A.16 excluded — coefficient plot now in pathreat.analysis.bytaxa.fig.R or similar
# ############################
# ### coefficient figure #####
# ############################
# {
#   cat("\n=== Creating coefficient figure ===\n")
#
#   plot_data <- estimate_results %>%
#     mutate(
#       pct_ci_low = ci_low / control_mean * 100,
#       pct_ci_high = ci_high / control_mean * 100,
#       model = factor(model, levels = rev(model))
#     )
#
#   fig <- ggplot(plot_data, aes(x = pct_of_control,
#                                y = model)) +
#     geom_vline(xintercept = 0, linewidth = 0.4, color = "grey40") +
#     geom_linerange(aes(xmin = pct_ci_low, xmax = pct_ci_high),
#                    linewidth = 0.8,
#                    color = "grey30") +
#     geom_point(size = 2.8, color = "black") +
#     scale_x_continuous(breaks = seq(-70, 10, by = 10)) +
#     labs(
#       x = "Effect on composite threat index (% of control mean)",
#       y = "Weighting scheme"
#     ) +
#     theme_minimal(base_size = 12) +
#     theme(
#       panel.grid.major.y = element_blank(),
#       panel.grid.minor = element_blank(),
#       axis.text.y = element_text(size = 11),
#       axis.title.x = element_text(margin = margin(t = 10)),
#       plot.margin = margin(15, 20, 10, 10)
#     )
#
#   figpath <- file.path(paths$figures_dir, "fig.global-threat.est.jpg")
#   ggsave(figpath, fig, width = 7, height = 3.5, dpi = 300)
#   cat("Saved:", figpath, "\n")
# }

cat("\n=== Done ===\n")
