############################################
### pathreat.randomrob.analysis.R ##########
### Matching-stage robustness sweep ########
### Single (main-text) estimation: country +
### biome FE with country-clustered SEs.
############################################

library(fixest)
library(fst)
library(dplyr)
library(ggplot2)

source("code/pathreat.analysis.config.R")


############################
### config #################
############################
{
  randomrob_dir <- "data/store/randomrob"

  # matched-sample subsample fractions (pair-level resampling of the full
  # production matched dataset — analogous to what Models 7/8 of
  # robustness.est.R used to do, but with multiple seeds for stability)
  match_sub_fracs <- c(`MATCH-SUB-20` = 0.20,
                       `MATCH-SUB-10` = 0.10)
  seeds <- 1001:1005

  # row labels for figure/table (left-to-right ordering on the x-axis)
  row_order <- c("REF-PROD", "MATCH-SUB-20", "MATCH-SUB-10", "REF-SUB",
                 "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8")
  row_labels <- c(
    "REF-PROD"     = "Main text (full sample)",
    "MATCH-SUB-20" = "20% matched subsample",
    "MATCH-SUB-10" = "10% matched subsample",
    "REF-SUB"      = "1% unmatched, baseline matching",
    "M1"           = "1% unmatched, parsimonious covariates",
    "M2"           = "1% unmatched, expanded covariates",
    "M3"           = "1% unmatched, 1:2 matching ratio",
    "M4"           = "1% unmatched, 1:5 matching ratio",
    "M5"           = "1% unmatched, without replacement",
    "M6"           = "1% unmatched, propensity score",
    "M7"           = "1% unmatched, propensity score + caliper",
    "M8"           = "1% unmatched, coarsened exact matching"
  )
}


############################
### helper: run E0 on a ####
### matched-like df ########
############################

# E0: threat_composite ~ treat | country_rast + biome_raster, cluster=country_rast
run_e0 <- function(m) {
  e <- tryCatch(
    feols(threat_composite ~ treat | country_rast + biome_raster,
          data = m, cluster = "country_rast"),
    error = function(err) {
      message(sprintf("  E0 failed: %s", conditionMessage(err)))
      NULL
    }
  )
  if (is.null(e)) return(NULL)
  ct <- coeftable(e)
  if (!"treat" %in% rownames(ct)) return(NULL)
  row <- ct["treat", , drop = TRUE]
  estimate <- regression_estimate(e)
  data.frame(
    coef         = row[[1]],
    se           = row[[2]],
    t_stat       = row[[3]],
    pval         = row[[4]],
    ci_low       = estimate$ci_low,
    ci_high      = estimate$ci_high,
    n_obs        = e$nobs,
    control_mean = mean(m$threat_composite[m$treat == 0], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}


############################
### 1. matching variants ###
############################
{
  cat("=== running E0 over matched datasets from randomrob.data.R ===\n")

  match_files <- list.files(randomrob_dir,
                            pattern = "^match\\.seed\\d+\\.variant-.*\\.Rds$",
                            full.names = TRUE)
  if (length(match_files) == 0) {
    stop(sprintf("No matched files in %s — run randomrob.data.R first",
                 randomrob_dir))
  }
  cat(sprintf("  found %d matched files\n", length(match_files)))

  results_match <- lapply(match_files, function(f) {
    base <- basename(f)
    seed <- as.integer(sub(".*seed(\\d+)\\.variant-.*", "\\1", base))
    mv   <- sub(".*variant-([^.]+)\\.Rds$", "\\1", base)

    m <- readRDS(f)
    n_pairs <- length(unique(m$matched_pair_id))
    cat(sprintf("  [seed %d, %s] n=%s, pairs=%s\n",
                seed, mv,
                format(nrow(m), big.mark = ","),
                format(n_pairs, big.mark = ",")))

    res <- run_e0(m)
    if (is.null(res)) return(NULL)
    res$seed          <- seed
    res$match_variant <- mv
    res$n_pairs       <- n_pairs
    res
  })
  results_match <- results_match[!vapply(results_match, is.null, logical(1))]
  results_match <- do.call(rbind, results_match)
}


############################
### 2. matched-sample #####
### subsample variants ####
############################
{
  cat("\n=== running E0 on random subsamples of the production matched sample ===\n")

  # load the production matched dataset once; select only columns needed for E0
  cat("  loading merge.matched.fst...\n")
  d_match <- read_fst(paths$data_matched,
                      columns = c("matched_pair_id", "treat",
                                  "country_rast", "biome_raster",
                                  "threat_composite"))
  cat(sprintf("  loaded %s rows\n", format(nrow(d_match), big.mark = ",")))

  pair_index <- build_matched_pair_index(d_match)
  sample_mask <- matched_pair_sample_mask(
    d_match,
    "threat_composite",
    pair_index
  )
  d_match <- d_match[sample_mask, , drop = FALSE]
  pair_index <- build_matched_pair_index(d_match)
  if (nrow(d_match) != 2L * pair_index$n_pairs) {
    stop("Production matched sample is not pair-complete")
  }

  results_matchsub <- do.call(rbind, lapply(names(match_sub_fracs), function(lab) {
    frac <- match_sub_fracs[[lab]]
    do.call(rbind, lapply(seeds, function(s) {
      set.seed(s)
      n_sub_pairs <- max(1L, round(pair_index$n_pairs * frac))
      sampled_pair_ids <- sample(pair_index$pair_ids, size = n_sub_pairs)
      m <- d_match[d_match$matched_pair_id %in% sampled_pair_ids, , drop = FALSE]
      n_treated <- sum(m$treat == 1)
      n_control <- sum(m$treat == 0)
      if (n_treated != n_sub_pairs || n_control != n_sub_pairs ||
          nrow(m) != 2L * n_sub_pairs) {
        stop(sprintf("%s / seed %d: matched subsample splits pairs", lab, s))
      }
      cat(sprintf("  [seed %d, %s] n=%s, pairs=%s\n",
                  s, lab, format(nrow(m), big.mark = ","),
                  format(n_sub_pairs, big.mark = ",")))
      res <- run_e0(m)
      if (is.null(res)) return(NULL)
      if (res$n_obs != nrow(m)) {
        stop(sprintf("%s / seed %d: E0 changed the canonical sample", lab, s))
      }
      res$seed          <- s
      res$match_variant <- lab
      res$n_pairs       <- n_sub_pairs
      res
    }))
  }))

  rm(d_match, pair_index, sample_mask)
  gc(verbose = FALSE)
}


############################
### 3. reference row: ######
###    REF-PROD ############
############################
{
  prod <- readRDS(paths$est_global)
  prod <- prod[prod$variable == "threat_composite", ]
  if (nrow(prod) != 1) {
    stop("Expected exactly one threat_composite row in global.est.Rds")
  }
  prod_n_pairs <- prod$n_obs / 2
  if (!is.finite(prod_n_pairs) || prod_n_pairs != floor(prod_n_pairs)) {
    stop("REF-PROD does not report an even pair-complete observation count")
  }
  ref_prod <- data.frame(
    seed          = NA_integer_,
    match_variant = "REF-PROD",
    coef          = prod$coef,
    se            = prod$se,
    t_stat        = prod$t_stat,
    pval          = prod$pval,
    ci_low        = prod$ci_low,
    ci_high       = prod$ci_high,
    n_obs         = prod$n_obs,
    control_mean  = prod$control_mean,
    n_pairs       = as.integer(prod_n_pairs),
    stringsAsFactors = FALSE
  )
}


############################
### 4. combine + save #####
############################
{
  col_order <- c("seed", "match_variant", "coef", "se", "t_stat", "pval",
                 "ci_low", "ci_high", "n_obs", "control_mean", "n_pairs")
  results_all <- rbind(
    results_match[, col_order],
    results_matchsub[, col_order],
    ref_prod[, col_order]
  )
  rownames(results_all) <- NULL

  saveRDS(results_all,
          file.path(paths$results_dir, "pathreat.randomrob.est.Rds"),
          compress = FALSE)
  cat(sprintf("\n  saved: %spathreat.randomrob.est.Rds (%d rows)\n",
              paths$results_dir, nrow(results_all)))
}


############################
### 5. figure ##############
############################
{
  cat("\n=== building figure ===\n")

  # relabel M0 → REF-SUB so its range is computed over seeds
  plot_df <- results_all %>%
    mutate(match_variant = ifelse(match_variant == "M0",
                                  "REF-SUB", match_variant)) %>%
    filter(match_variant %in% row_order) %>%
    mutate(match_variant = factor(match_variant,
                                  levels = row_order,
                                  labels = row_labels[row_order]))

  plot_sum <- plot_df %>%
    group_by(match_variant) %>%
    summarise(
      coef_med    = median(coef, na.rm = TRUE),
      coef_min    = min(coef, na.rm = TRUE),
      coef_max    = max(coef, na.rm = TRUE),
      ci_low_med  = median(ci_low, na.rm = TRUE),
      ci_high_med = median(ci_high, na.rm = TRUE),
      .groups     = "drop"
    )

  # mark reference / variant rows for fill aesthetics
  ref_rows <- c(row_labels["REF-PROD"], row_labels["REF-SUB"])
  plot_sum$row_kind <- ifelse(plot_sum$match_variant %in% ref_rows,
                              "reference", "variant")

  prod_coef <- results_all$coef[results_all$match_variant == "REF-PROD"][1]

  p <- ggplot(plot_sum,
              aes(x = match_variant, y = coef_med)) +
    geom_hline(yintercept = 0,
               linetype = "solid", color = "black", linewidth = 0.4) +
    geom_hline(yintercept = prod_coef,
               linetype = "dashed", color = "grey40", linewidth = 0.3) +
    geom_errorbar(aes(ymin = ci_low_med, ymax = ci_high_med),
                  width = 0.15, color = "grey30", linewidth = 0.4) +
    geom_point(aes(shape = row_kind),
               color = "black", size = 2.4, fill = "black") +
    scale_shape_manual(values = c(reference = 18, variant = 16),
                       guide = "none") +
    labs(x = NULL,
         y = "ATT on threat_composite",
         title = NULL) +
    theme_classic(base_size = 10) +
    theme(
      panel.border = element_blank(),
      axis.line.y  = element_line(color = "black", linewidth = 0.3),
      axis.line.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.x  = element_text(angle = 35, hjust = 1),
      plot.margin  = margin(6, 10, 6, 6)
    )

  ggsave(filename = file.path(paths$figures_dir,
                              "fig.randomrob.coef.est.jpg"),
         plot = p, width = 9, height = 5, dpi = 300)
  cat(sprintf("  saved: %sfig.randomrob.coef.est.jpg\n", paths$figures_dir))
}


############################
### 6. LaTeX table #########
############################
{
  cat("\n=== building LaTeX table ===\n")

  star_fn <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.01) "\\sym{***}"
    else if (p < 0.05) "\\sym{**}"
    else if (p < 0.1)  "\\sym{*}"
    else ""
  }

  tab_df <- results_all %>%
    mutate(match_variant = ifelse(match_variant == "M0",
                                  "REF-SUB", match_variant)) %>%
    filter(match_variant %in% row_order) %>%
    group_by(match_variant) %>%
    summarise(
      coef_med = median(coef, na.rm = TRUE),
      se_med   = median(se,   na.rm = TRUE),
      pval_med = median(pval, na.rm = TRUE),
      n_seeds  = sum(!is.na(seed)),
      .groups  = "drop"
    )

  row_tex_label <- c(
    "REF-PROD"     = "Main text (full sample)",
    "MATCH-SUB-20" = "20\\% matched subsample",
    "MATCH-SUB-10" = "10\\% matched subsample",
    "REF-SUB"      = "1\\% unmatched, baseline matching",
    "M1"           = "1\\% unmatched, parsimonious covariates",
    "M2"           = "1\\% unmatched, expanded covariates (with climate quadratics)",
    "M3"           = "1\\% unmatched, 1:2 matching ratio",
    "M4"           = "1\\% unmatched, 1:5 matching ratio",
    "M5"           = "1\\% unmatched, without replacement",
    "M6"           = "1\\% unmatched, propensity score",
    "M7"           = "1\\% unmatched, propensity score (with 0.25 SD caliper)",
    "M8"           = "1\\% unmatched, coarsened exact matching"
  )

  # two TeX rows per variant: coef row (with label, stars) and SE row (in parens)
  tex_rows <- character()
  for (mv in row_order) {
    sub <- tab_df[tab_df$match_variant == mv, ]
    if (nrow(sub) == 0) {
      tex_rows <- c(tex_rows,
                    sprintf("%s & {--} \\\\", row_tex_label[mv]),
                    "             & {} \\\\")
      next
    }
    coef_cell <- sprintf("%.3f%s", sub$coef_med, star_fn(sub$pval_med))
    se_cell   <- sprintf("(%.3f)", sub$se_med)
    tex_rows <- c(tex_rows,
                  sprintf("%s & %s \\\\", row_tex_label[mv], coef_cell),
                  sprintf("             & %s \\\\", se_cell))
  }
  # strip trailing \\ from last row (siunitx S-column convention)
  tex_rows[length(tex_rows)] <- sub(" *\\\\\\\\$", "",
                                    tex_rows[length(tex_rows)])

  out_tex <- "pub/tables/tab.randomrob.tex"
  writeLines(tex_rows, out_tex)
  cat(sprintf("  saved: %s\n", out_tex))
}

cat("\n=== randomrob.analysis complete ===\n")
