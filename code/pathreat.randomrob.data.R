############################################
### pathreat.randomrob.data.R ##############
### Random subsample + matching-variant grid
############################################

library(Matching)
library(MatchIt)
library(dplyr)
library(tidyr)
library(fst)
library(parallel)
library(pbapply)

source("code/pathreat.analysis.config.R")


############################
### config #################
############################
{
  # random subsample settings
  seeds          <- 1001:1005
  subsample_frac <- 0.01   # ~1% of unmatched pixels
  k_min          <- 50     # per-cell floor: keep all rows if n(cell) < k_min

  # matching-variant grid
  # method: "nn" = Matching::Match (nearest-neighbor); "matchit" = MatchIt::matchit
  variants_tbl <- data.frame(
    id       = c("M0", "M1", "M2", "M3", "M4", "M5",
                 "M6", "M7", "M8"),
    label    = c("Base: Mah 1:1 repl",
                 "Parsimonious covs",
                 "Expanded covs (+quad)",
                 "Mah 1:2 repl",
                 "Mah 1:5 repl",
                 "Mah 1:1 no repl",
                 "PSM 1:1 repl",
                 "PSM 1:1 + caliper",
                 "CEM"),
    method   = c("nn", "nn", "nn", "nn", "nn", "nn",
                 "matchit", "matchit", "matchit"),
    distance = c("mahalanobis", "mahalanobis", "mahalanobis",
                 "mahalanobis", "mahalanobis", "mahalanobis",
                 "glm", "glm", "cem"),
    ratio    = c(1, 1, 1, 2, 5, 1, 1, 1, NA_integer_),
    replace  = c(TRUE, TRUE, TRUE, TRUE, TRUE, FALSE,
                 TRUE, TRUE, NA),
    caliper  = c(NA, NA, NA, NA, NA, NA, NA, 0.25, NA),
    cov_set  = c("full", "parsimonious", "expanded",
                 "full", "full", "full",
                 "full", "full", "full"),
    stringsAsFactors = FALSE
  )

  # covariate sets (resolve to column names)
  cov_sets <- list(
    full         = mlist,
    parsimonious = c("elevation", "slope", "access", "forest2000_parea"),
    # note: expanded adds quadratics constructed at runtime
    expanded     = mlist
  )

  # output directory
  randomrob_dir <- "data/store/randomrob"
  if (!dir.exists(randomrob_dir)) {
    dir.create(randomrob_dir, recursive = TRUE)
  }

  cat("=== randomrob.data config ===\n")
  cat(sprintf("  seeds:          %s\n", paste(seeds, collapse = ", ")))
  cat(sprintf("  subsample_frac: %.3f\n", subsample_frac))
  cat(sprintf("  k_min:          %d\n", k_min))
  cat(sprintf("  variants:       %d (%s)\n",
              nrow(variants_tbl), paste(variants_tbl$id, collapse = ", ")))
  cat(sprintf("  output dir:     %s\n", randomrob_dir))
}


############################
### load + filter data #####
############################
{
  cat("\n=== Loading unmatched data ===\n")
  match_cols <- c("pixel_id", "country_rast", "country", "biome_raster",
                  "treat", "buffer_5km", "pa_year_designated", "admin1_id",
                  "threat_composite", mlist)
  d0 <- read_fst(paths$data_unmatched, columns = match_cols)
  cat(sprintf("  loaded: %s rows\n", format(nrow(d0), big.mark = ",")))

  # same pre-match filters as matching.R
  mlist_check <- c("biome_raster", "country_rast", mlist)
  i <- which(rowSums(is.na(d0[mlist_check])) > 0)
  if (length(i) > 0) d0 <- d0[-i, ]

  i <- which(d0$pa_year_designated > 2020); if (length(i) > 0) d0 <- d0[-i, ]
  i <- which(d0$pa_year_designated <= 2000); if (length(i) > 0) d0 <- d0[-i, ]
  i <- which(d0$treat == 0 & d0$buffer_5km == 1); if (length(i) > 0) d0 <- d0[-i, ]
  i <- which(is.na(d0$admin1_id)); if (length(i) > 0) d0 <- d0[-i, ]

  # bigregion (admin1 for big 4 countries, else 0)
  big_countries <- c("Brazil", "Australia", "Canada", "Russia")
  is_big <- d0$country %in% big_countries
  d0$bigregion <- 0L
  d0$bigregion[is_big] <- d0$admin1_id[is_big]
  rm(is_big)

  # drop threat_composite NAs (can't estimate without outcome)
  i <- which(is.na(d0$threat_composite))
  if (length(i) > 0) d0 <- d0[-i, ]

  cat(sprintf("  after filters: %s rows (treat=1: %s, treat=0: %s)\n",
              format(nrow(d0), big.mark = ","),
              format(sum(d0$treat == 1), big.mark = ","),
              format(sum(d0$treat == 0), big.mark = ",")))

  # add quadratic climate columns for the "expanded" covariate set
  d0$annual_total_precipitation_sq <- d0$annual_total_precipitation^2
  d0$temperature_sq                <- d0$temperature^2
  cov_sets$expanded <- c(mlist,
                         "annual_total_precipitation_sq",
                         "temperature_sq")
}


############################
### helpers ################
############################

# stratified subsample: within each (country, biome, bigregion, treat) cell,
# sample max(k_min, ceil(n*frac)) rows (or keep all if n < k_min)
stratified_subsample <- function(d, frac, k_min, seed) {
  set.seed(seed)
  cell_key <- paste(d$country_rast, d$biome_raster,
                    d$bigregion, d$treat, sep = ".")
  idx_by_cell <- split(seq_len(nrow(d)), cell_key)
  sampled <- lapply(idx_by_cell, function(ix) {
    n <- length(ix)
    n_keep <- min(n, max(k_min, ceiling(n * frac)))
    if (n_keep >= n) ix else sample(ix, n_keep)
  })
  d[unlist(sampled, use.names = FALSE), , drop = FALSE]
}


# run Matching::Match on a single stratum; returns matched pairs in
# (pixel_id, treat, matched_pair_id, weight) form; NULL on failure/empty
match_stratum_nn <- function(stratum_rows, covs, variant) {
  m.c <- as.matrix(stratum_rows[, covs, drop = FALSE])
  m.t <- stratum_rows$treat
  pids <- stratum_rows$pixel_id

  # drop zero-variance covariates (stratum-specific safety net)
  keep_col <- apply(m.c, 2, function(v) sd(v, na.rm = TRUE) > 0)
  if (any(!keep_col)) m.c <- m.c[, keep_col, drop = FALSE]
  if (ncol(m.c) == 0) return(NULL)

  mr <- tryCatch(
    Match(Tr = m.t, X = m.c,
          replace = variant$replace,
          M       = as.integer(variant$ratio),
          Weight  = 2,
          ties    = FALSE),
    error = function(e) NULL
  )
  if (is.null(mr) || length(mr$index.treated) == 0) return(NULL)

  idx_t <- mr$index.treated
  idx_c <- mr$index.control
  w     <- mr$weights
  np    <- length(idx_t)

  data.frame(
    pixel_id        = c(pids[idx_t], pids[idx_c]),
    treat           = c(rep(1L, np),  rep(0L, np)),
    matched_pair_id = rep(seq_len(np), times = 2),
    weight          = c(w, w),
    row_in_stratum  = c(idx_t, idx_c),
    stringsAsFactors = FALSE
  )
}


# run MatchIt::matchit on a single stratum.
# For method="nearest" use get_matches() — it duplicates controls matched to
# multiple treated and guarantees a subclass column. For method="cem" use
# match.data() (get_matches is documented as unsuitable for CEM).
match_stratum_matchit <- function(stratum_rows, covs, variant) {
  stratum_rows <- stratum_rows[, c("pixel_id", "treat", covs), drop = FALSE]
  stratum_rows <- stratum_rows[complete.cases(stratum_rows), ]
  if (sum(stratum_rows$treat == 1) == 0 || sum(stratum_rows$treat == 0) == 0) {
    return(NULL)
  }

  fml <- as.formula(paste("treat ~", paste(covs, collapse = " + ")))

  args <- list(formula = fml, data = stratum_rows)
  if (variant$distance == "glm") {
    args$method   <- "nearest"
    args$distance <- "glm"
    args$ratio    <- as.integer(variant$ratio)
    args$replace  <- variant$replace
    if (!is.na(variant$caliper)) args$caliper <- variant$caliper
  } else if (variant$distance == "cem") {
    args$method   <- "cem"
  } else {
    stop(sprintf("Unsupported MatchIt distance: %s", variant$distance))
  }

  m <- tryCatch(suppressWarnings(do.call(MatchIt::matchit, args)),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)

  if (variant$distance == "cem") {
    md <- tryCatch(MatchIt::match.data(m, data = stratum_rows,
                                       weights = "weight",
                                       subclass = "matched_pair_id",
                                       drop.unmatched = TRUE),
                   error = function(e) NULL)
    if (is.null(md) || nrow(md) == 0 || is.null(md$matched_pair_id)) {
      return(NULL)
    }
    return(data.frame(
      pixel_id        = md$pixel_id,
      treat           = md$treat,
      matched_pair_id = as.integer(md$matched_pair_id),
      weight          = md$weight,
      row_in_stratum  = NA_integer_,
      stringsAsFactors = FALSE
    ))
  } else {
    mg <- tryCatch(MatchIt::get_matches(m, data = stratum_rows,
                                        weights = "weight",
                                        subclass = "matched_pair_id",
                                        id = "row_id_match"),
                   error = function(e) NULL)
    if (is.null(mg) || nrow(mg) == 0 || is.null(mg$matched_pair_id)) {
      return(NULL)
    }
    return(data.frame(
      pixel_id        = mg$pixel_id,
      treat           = mg$treat,
      matched_pair_id = as.integer(as.character(mg$matched_pair_id)),
      weight          = mg$weight,
      row_in_stratum  = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }
}


# balance diagnostics: weighted |SMD| per covariate on the matched sample
# `matched_df` must already contain the covariate columns (merged upstream).
compute_balance <- function(matched_df, covs) {
  wmean <- function(x, w) sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
  wvar  <- function(x, w) {
    m <- wmean(x, w)
    sum(w * (x - m)^2, na.rm = TRUE) / sum(w, na.rm = TRUE)
  }
  covs_present <- intersect(covs, names(matched_df))
  res <- lapply(covs_present, function(v) {
    xt <- matched_df[[v]][matched_df$treat == 1]
    wt <- matched_df$weight[matched_df$treat == 1]
    xc <- matched_df[[v]][matched_df$treat == 0]
    wc <- matched_df$weight[matched_df$treat == 0]
    mu_t <- wmean(xt, wt); mu_c <- wmean(xc, wc)
    sd_pooled <- sqrt(0.5 * (wvar(xt, wt) + wvar(xc, wc)))
    smd <- if (is.finite(sd_pooled) && sd_pooled > 0) {
      (mu_t - mu_c) / sd_pooled
    } else NA_real_
    data.frame(covariate = v, smd = smd, abs_smd = abs(smd))
  })
  do.call(rbind, res)
}


############################
### driver loop ############
############################

for (seed in seeds) {
  t_seed <- Sys.time()
  cat(sprintf("\n=== Seed %d ===\n", seed))

  # draw subsample (once per seed, reused across variants)
  d <- stratified_subsample(d0, subsample_frac, k_min, seed)
  cat(sprintf("  subsample: %s rows (treat=1: %s, treat=0: %s)\n",
              format(nrow(d), big.mark = ","),
              format(sum(d$treat == 1), big.mark = ","),
              format(sum(d$treat == 0), big.mark = ",")))

  # save subsample for reproducibility
  saveRDS(d, sprintf("%s/subsample.seed%d.Rds", randomrob_dir, seed),
          compress = FALSE)

  # split row indices by stratum once per seed
  strata_key <- paste(d$country_rast, d$biome_raster, d$bigregion, sep = ".")
  strata_idx <- split(seq_len(nrow(d)), strata_key)
  # only keep strata that have both treated and controls
  strata_idx <- strata_idx[vapply(strata_idx, function(r) {
    any(d$treat[r] == 1) && any(d$treat[r] == 0)
  }, logical(1))]
  cat(sprintf("  valid strata: %d\n", length(strata_idx)))

  for (v in seq_len(nrow(variants_tbl))) {
    variant <- variants_tbl[v, ]
    t_var <- Sys.time()
    covs <- cov_sets[[variant$cov_set]]

    out_file <- sprintf("%s/match.seed%d.variant-%s.Rds",
                        randomrob_dir, seed, variant$id)
    bal_file <- sprintf("%s/balance.seed%d.variant-%s.Rds",
                        randomrob_dir, seed, variant$id)

    cat(sprintf("  [%s] %s — covs=%s, ratio=%s, replace=%s, caliper=%s\n",
                variant$id, variant$label, variant$cov_set,
                variant$ratio, variant$replace, variant$caliper))

    # loop strata
    match_fun <- if (variant$method == "nn") match_stratum_nn else match_stratum_matchit

    strata_results <- lapply(seq_along(strata_idx), function(s) {
      rows <- strata_idx[[s]]
      stratum <- d[rows, , drop = FALSE]
      nt <- sum(stratum$treat == 1); nc <- sum(stratum$treat == 0)
      if (nt == 0 || nc == 0) return(NULL)

      res <- match_fun(stratum, covs, variant)
      if (is.null(res)) return(NULL)

      res$country_rast <- stratum$country_rast[1]
      res$biome_raster <- stratum$biome_raster[1]
      res$bigregion    <- stratum$bigregion[1]
      res$stratum_key  <- names(strata_idx)[s]
      # globally unique pair id: stratum_seq + local id
      res$matched_pair_id <- paste0(res$stratum_key, "_",
                                    res$matched_pair_id)
      res
    })
    n_strata_matched <- sum(!vapply(strata_results, is.null, logical(1)))
    strata_results <- strata_results[!vapply(strata_results, is.null, logical(1))]
    if (length(strata_results) == 0) {
      cat(sprintf("    variant %s produced no matches — skipping\n",
                  variant$id))
      next
    }

    matched <- do.call(rbind, strata_results)
    matched$matched_pair_id <- as.integer(factor(matched$matched_pair_id))
    rm(strata_results); gc(verbose = FALSE)

    # merge threat_composite and covariates back for estimation
    keep_cols <- c("pixel_id", "threat_composite", mlist,
                   "annual_total_precipitation_sq", "temperature_sq")
    keep_cols <- intersect(keep_cols, names(d))
    matched <- merge(matched, d[, keep_cols], by = "pixel_id", all.x = TRUE)

    # compute delta per pair (mean within pair × treat arm, to handle 1:k and CEM)
    matched <- matched %>%
      group_by(matched_pair_id) %>%
      mutate(
        delta = {
          mt <- mean(threat_composite[treat == 1], na.rm = TRUE)
          mc <- mean(threat_composite[treat == 0], na.rm = TRUE)
          mt - mc
        },
        tc_threat_composite = mean(threat_composite[treat == 0], na.rm = TRUE)
      ) %>%
      ungroup() %>%
      data.frame()

    # balance diagnostics
    bal <- compute_balance(matched, covs)
    bal$seed    <- seed
    bal$variant <- variant$id

    saveRDS(matched, out_file, compress = FALSE)
    saveRDS(bal, bal_file, compress = FALSE)

    cat(sprintf("    → pairs: %s | strata matched: %d | mean |SMD|: %.3f | %.1fs\n",
                format(length(unique(matched$matched_pair_id)), big.mark = ","),
                n_strata_matched,
                mean(bal$abs_smd, na.rm = TRUE),
                as.numeric(difftime(Sys.time(), t_var, units = "secs"))))
  }

  cat(sprintf("  seed %d done in %.1f min\n",
              seed,
              as.numeric(difftime(Sys.time(), t_seed, units = "mins"))))
}

cat("\n=== randomrob.data complete ===\n")
