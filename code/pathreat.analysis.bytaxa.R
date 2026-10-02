##################################################
### pathreat.analysis.bytaxa.R ###################
### Central taxon-median producer ################
##################################################
#
# Computes every taxon-level species summary and bootstrap once. Figure and
# headline-statistics scripts consume the saved bundle rather than resampling
# species independently.
#
# Inputs:
#   results/pathreat.byspecies.est.Rds
#   data/store/pathreat.data.species.Rds
#
# Output:
#   results/pathreat.bytaxa.est.Rds
#

library(dplyr)

source("code/pathreat.analysis.config.R")

############################
### helpers ################
############################

weighted_median <- function(x, w) {
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cumw <- cumsum(w) / sum(w)
  x[which(cumw >= 0.5)[1]]
}

two_sided_sign_p <- function(draws) {
  min(1, 2 * min(mean(draws >= 0), mean(draws <= 0)))
}

############################
### load + validate ########
############################
{
cat("=== Loading canonical species estimates ===\n")

sp_path <- file.path(paths$results_dir, "pathreat.byspecies.est.Rds")
if (!file.exists(sp_path)) {
  stop("Missing ", sp_path, ". Run code/pathreat.analysis.byspecies.R first.")
}
sp_est <- readRDS(sp_path)

required_est <- c(
  "species_id",
  "taxon",
  "delta_s",
  "se_delta_s",
  "p_delta_s",
  "pct_effect_s",
  "control_mean_s",
  "status"
)
missing_est <- setdiff(required_est, names(sp_est))
if (length(missing_est) > 0) {
  stop("Canonical species estimates are missing: ",
       paste(missing_est, collapse = ", "))
}
if (anyDuplicated(sp_est$species_id)) {
  stop("Canonical species estimates contain duplicated species_id values.")
}

sp_meta_path <- "data/store/pathreat.data.species.Rds"
sp_meta <- readRDS(sp_meta_path)
weight_cols <- c("tbl_median", "ed_median", "edge_median")
required_meta <- c("species_id", "taxon", weight_cols)
missing_meta <- setdiff(required_meta, names(sp_meta))
if (length(missing_meta) > 0) {
  stop("Species metadata are missing: ", paste(missing_meta, collapse = ", "))
}
if (anyDuplicated(sp_meta[, c("species_id", "taxon")])) {
  stop("Species metadata contain duplicated species_id/taxon rows.")
}

sp_est <- sp_est %>%
  left_join(
    sp_meta[, required_meta],
    by = c("species_id", "taxon")
  )

sp_est$eligible_estimated <- with(
  sp_est,
  status == "estimated" &
    is.finite(delta_s) &
    is.finite(se_delta_s) &
    is.finite(p_delta_s) &
    is.finite(pct_effect_s) &
    is.finite(control_mean_s)
)

if (anyNA(sp_est$status) || anyNA(sp_est$eligible_estimated)) {
  stop("Canonical species estimates contain missing status or eligibility values.")
}

canonical_numeric_fields <- c(
  "delta_s",
  "se_delta_s",
  "p_delta_s",
  "pct_effect_s",
  "control_mean_s"
)
canonical_numeric_finite <- Reduce(
  `&`,
  lapply(canonical_numeric_fields, function(field) is.finite(sp_est[[field]]))
)
canonical_estimated <- sp_est$status == "estimated"

if (any(canonical_estimated & !canonical_numeric_finite)) {
  stop("Canonical estimated-species rows contain nonfinite model values.")
}

canonical_eligible <- canonical_estimated & canonical_numeric_finite
if (!identical(sp_est$eligible_estimated, canonical_eligible)) {
  stop("Taxon eligibility is inconsistent with canonical species estimates.")
}

n_eligible_estimated <- sum(canonical_eligible)
if (n_eligible_estimated == 0L) {
  stop("No regression-estimated species have finite canonical values.")
}

cat(sprintf(
  "  Species: %s total; %s regression-estimated and eligible\n",
  format(nrow(sp_est), big.mark = ","),
  format(n_eligible_estimated, big.mark = ",")
))
}

############################
### settings + counts ######
############################
{
seed <- 42L
n_boot <- 2000L
effect_taxa <- c("amphibian", "bird", "mammal", "reptile")
level_taxa <- c("mammal", "bird", "amphibian", "reptile")

settings <- list(
  estimand = paste(
    "Species-range treatment effect from matched pairs whose protected pixel",
    "lies inside the species range, retaining both pair members and using",
    "adaptive country/biome fixed effects"
  ),
  seed = seed,
  bootstrap_draws = n_boot,
  eligibility_rule = paste(
    "status == 'estimated' and finite delta_s, se_delta_s, p_delta_s,",
    "pct_effect_s, and control_mean_s; single-pair records are excluded"
  ),
  weighting_definitions = c(
    equal = "Equal species weight within taxon",
    tbl = "Weighted median using positive finite terminal-branch-length median",
    ed = "Weighted median using positive finite evolutionary-distinctiveness median",
    edge = "Weighted median using positive finite EDGE2 median"
  ),
  percentage_definition = "100 * delta_s / control_mean_s",
  effect_taxon_order = effect_taxa,
  level_taxon_order = level_taxa
)

counts_list <- lapply(effect_taxa, function(tx) {
  d <- sp_est[sp_est$taxon == tx, , drop = FALSE]
  eligible <- d$eligible_estimated

  data.frame(
    taxon = tx,
    n_species_total = nrow(d),
    n_species_estimated = sum(eligible),
    n_species_equal = sum(eligible),
    n_species_tbl = sum(eligible & is.finite(d$tbl_median) & d$tbl_median > 0),
    n_species_ed = sum(eligible & is.finite(d$ed_median) & d$ed_median > 0),
    n_species_edge = sum(eligible & is.finite(d$edge_median) & d$edge_median > 0),
    stringsAsFactors = FALSE
  )
})
counts <- do.call(rbind, counts_list)
rownames(counts) <- NULL

if (anyDuplicated(counts$taxon) || !setequal(counts$taxon, effect_taxa)) {
  stop("Taxon counts do not contain exactly one row for each configured taxon.")
}

eligible_taxa <- sp_est$taxon[sp_est$eligible_estimated]
if (anyNA(eligible_taxa) || !all(eligible_taxa %in% effect_taxa)) {
  stop("Eligible canonical species contain missing or unsupported taxa.")
}

canonical_taxon_counts <- table(factor(eligible_taxa, levels = effect_taxa))
observed_taxon_counts <- setNames(
  counts$n_species_estimated,
  counts$taxon
)[effect_taxa]

if (!identical(
      as.integer(observed_taxon_counts),
      as.integer(canonical_taxon_counts)
    )) {
  stop("Per-taxon estimated-species counts disagree with canonical species rows.")
}
if (sum(counts$n_species_estimated) != n_eligible_estimated) {
  stop("Summed taxon counts disagree with the canonical eligible-species count.")
}
if (any(counts$n_species_equal != counts$n_species_estimated)) {
  stop("Equal-weight taxon counts disagree with estimated-species counts.")
}
}

##########################################
### raw + percentage taxon effects #######
##########################################
{
cat("\n=== Taxon effect medians and percentile intervals ===\n")

weight_specs <- list(
  equal = NULL,
  tbl = "tbl_median",
  ed = "ed_median",
  edge = "edge_median"
)

# Preserve the production figure's exact RNG sequence and loop order.
set.seed(seed)
effect_rows <- list()

for (tx in effect_taxa) {
  d_base <- sp_est[
    sp_est$taxon == tx & sp_est$eligible_estimated,
    ,
    drop = FALSE
  ]

  for (weighting in names(weight_specs)) {
    weight_col <- weight_specs[[weighting]]
    d <- d_base

    if (!is.null(weight_col)) {
      d <- d[
        is.finite(d[[weight_col]]) & d[[weight_col]] > 0,
        ,
        drop = FALSE
      ]
    }
    if (nrow(d) < 2) {
      stop("Fewer than two eligible species for ", tx, " / ", weighting, ".")
    }

    raw_values <- d$delta_s
    pct_values <- d$pct_effect_s
    weights <- if (is.null(weight_col)) NULL else d[[weight_col]]

    bootstrap_draws <- replicate(n_boot, {
      idx <- sample.int(nrow(d), replace = TRUE)
      if (is.null(weights)) {
        c(
          raw = median(raw_values[idx]),
          pct = median(pct_values[idx])
        )
      } else {
        c(
          raw = weighted_median(raw_values[idx], weights[idx]),
          pct = weighted_median(pct_values[idx], weights[idx])
        )
      }
    })

    raw_draws <- bootstrap_draws["raw", ]
    pct_draws <- bootstrap_draws["pct", ]
    raw_ci <- quantile(raw_draws, c(0.025, 0.975))
    pct_ci <- quantile(pct_draws, c(0.025, 0.975))

    raw_median <- if (is.null(weights)) {
      median(raw_values)
    } else {
      weighted_median(raw_values, weights)
    }
    pct_median <- if (is.null(weights)) {
      median(pct_values)
    } else {
      weighted_median(pct_values, weights)
    }

    effect_rows[[length(effect_rows) + 1L]] <- data.frame(
      taxon = tx,
      weighting = weighting,
      n_species = nrow(d),
      raw_median = raw_median,
      raw_ci_low = unname(raw_ci[1]),
      raw_ci_high = unname(raw_ci[2]),
      raw_prob_negative = mean(raw_draws < 0),
      raw_prob_positive = mean(raw_draws > 0),
      raw_sign_p = two_sided_sign_p(raw_draws),
      pct_median = pct_median,
      pct_ci_low = unname(pct_ci[1]),
      pct_ci_high = unname(pct_ci[2]),
      pct_prob_negative = mean(pct_draws < 0),
      pct_prob_positive = mean(pct_draws > 0),
      pct_sign_p = two_sided_sign_p(pct_draws),
      stringsAsFactors = FALSE
    )
  }
}

effects <- do.call(rbind, effect_rows)
rownames(effects) <- NULL
}

##########################################
### marginal control/protected levels ####
##########################################
{
cat("\n=== Equal-species marginal level medians ===\n")

# Preserve the production level bars' exact RNG sequence and taxon order.
set.seed(seed)
level_rows <- list()

for (tx in level_taxa) {
  d <- sp_est[
    sp_est$taxon == tx & sp_est$eligible_estimated,
    ,
    drop = FALSE
  ]

  control_values <- d$control_mean_s
  protected_values <- d$control_mean_s + d$delta_s
  raw_values <- d$delta_s

  bootstrap_draws <- replicate(n_boot, {
    idx <- sample.int(nrow(d), replace = TRUE)
    c(
      control = median(control_values[idx]),
      protected = median(protected_values[idx]),
      raw_effect = median(raw_values[idx])
    )
  })

  control_draws <- bootstrap_draws["control", ]
  protected_draws <- bootstrap_draws["protected", ]
  raw_draws <- bootstrap_draws["raw_effect", ]
  control_ci <- quantile(control_draws, c(0.025, 0.975))
  protected_ci <- quantile(protected_draws, c(0.025, 0.975))
  raw_ci <- quantile(raw_draws, c(0.025, 0.975))
  raw_sign_p <- max(two_sided_sign_p(raw_draws), 1 / n_boot)

  level_rows[[length(level_rows) + 1L]] <- data.frame(
    taxon = tx,
    n_species = nrow(d),
    control_median = median(control_values),
    control_ci_low = unname(control_ci[1]),
    control_ci_high = unname(control_ci[2]),
    control_approx_se = unname(diff(control_ci) / (2 * 1.96)),
    protected_median = median(protected_values),
    protected_ci_low = unname(protected_ci[1]),
    protected_ci_high = unname(protected_ci[2]),
    protected_approx_se = unname(diff(protected_ci) / (2 * 1.96)),
    raw_effect_median = median(raw_values),
    raw_effect_ci_low = unname(raw_ci[1]),
    raw_effect_ci_high = unname(raw_ci[2]),
    raw_effect_prob_negative = mean(raw_draws < 0),
    raw_effect_prob_positive = mean(raw_draws > 0),
    raw_effect_sign_p = raw_sign_p,
    stringsAsFactors = FALSE
  )
}

levels <- do.call(rbind, level_rows)
rownames(levels) <- NULL
}

############################
### save bundle ############
############################
{
taxon_bundle <- list(
  settings = settings,
  counts = counts,
  effects = effects,
  levels = levels
)

outpath <- file.path(paths$results_dir, "pathreat.bytaxa.est.Rds")
saveRDS(taxon_bundle, outpath)

cat("\nCounts by taxon:\n")
print(counts, row.names = FALSE)
cat("\nEqual-weight percentage medians:\n")
print(
  effects[
    effects$weighting == "equal",
    c("taxon", "n_species", "pct_median", "pct_ci_low", "pct_ci_high")
  ],
  row.names = FALSE
)
cat("\nSaved: ", outpath, "\n", sep = "")
}
