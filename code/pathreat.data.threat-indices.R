library(fst)

source("code/pathreat.analysis.config.R")



# read from merged data (all threats + GADM-derived country_rast)
d <- read_fst("data/store/pathreat.data.merge.fst",
              columns = c("pixel_id", "country_rast", threat_list_raw))


####################################
### resolve presence thresholds ####
####################################
{
cat("\n=== Resolving presence thresholds ===\n")
# dep_thresholds defined in config (variable, value, quantile_prob, dir)
# resolve quantile-based thresholds from unmatched sample
for (v in dep_thresholds$variable) {
  # v<-dep_thresholds$variable[1]
  qp <- dep_thresholds[v, "quantile_prob"]
  if (!is.na(qp)) {
    dep_thresholds[v, "value"] <- quantile(d[[v]], qp, na.rm = TRUE)
  }
}
cat("Resolved presence thresholds:\n")
print(dep_thresholds)
# saveRDS(dep_thresholds, "data/store/pathreat.data.1_dep_thresholds.Rds")
}



##########################################
### threat composite index ###############
### (global P1/P99 normalization) ########
##########################################
{
cat("\n=== Computing threat composite index (global P1/P99) ===\n")
# threat_list_raw and threat_categories defined in config

# step 1: compute global P1/P99 normalization params
norm_params_list <- list()
for (v in threat_list_raw) {
  cat("  Normalizing:", v, "\n")
  qs <- quantile(d[[v]], probs = c(0.01, 0.99), na.rm = TRUE)
  v_min <- qs[1]; v_max <- qs[2]; v_range <- v_max - v_min
  if (v_range == 0) {
    # sparse variable (P1=P99=0): use presence indicator as normalized value
    v_min <- NA; v_max <- NA; v_range <- NA
    thresh <- dep_thresholds[v, "value"]
    dir    <- dep_thresholds[v, "dir"]
    x <- d[[v]]
    d[[paste0("norm_", v)]] <- as.numeric(
      if (dir == "<") (!is.na(x) & x < thresh) else (!is.na(x) & x > thresh)
    )
    cat("    -> sparse (P1=P99=0), using presence indicator\n")
  } else {
    # standard winsorized min-max normalization, clip to [0, 1]
    d[[paste0("norm_", v)]] <- pmin(pmax((d[[v]] - v_min) / v_range, 0), 1)
    if (v == "swu") d[[paste0("norm_", v)]] <- 1 - d[[paste0("norm_", v)]]
  }

  norm_params_list[[v]] <- data.frame(
    variable  = v,
    v_min     = v_min,
    v_max     = v_max,
    v_range   = v_range,
    mean_norm = mean(d[[paste0("norm_", v)]], na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

norm_params <- do.call(rbind, norm_params_list)
rownames(norm_params) <- norm_params$variable
# saveRDS(norm_params, "data/store/pathreat.data.1_norm_params.Rds")
cat("Saved normalization params:\n")
print(norm_params)

# step 2: country-specific presence-proportion weights
# compute country-level presence proportions (small: ~158 rows × 14 cols)
cat("  Computing country-specific presence weights...\n")
countries <- sort(unique(d$country_rast))
cw_mat <- matrix(0, nrow = length(countries), ncol = length(threat_list_raw),
                 dimnames = list(as.character(countries), threat_list_raw))
for (v in threat_list_raw) {
  thresh <- dep_thresholds[v, "value"]
  dir    <- dep_thresholds[v, "dir"]
  x <- d[[v]]
  pres <- if (dir == "<") (!is.na(x) & x < thresh) else (!is.na(x) & x > thresh)
  # tapply to get per-country means
  cm <- tapply(pres, d$country_rast, mean, na.rm = TRUE)
  cw_mat[as.character(names(cm)), v] <- cm
}
rm(pres, cm, x); gc()

# convert raw prevalence fractions into weights nested within IUCN categories
within_cat_weight_mat <- cw_mat
cat_weight_mat <- sapply(threat_categories, function(vars) {
  rowSums(cw_mat[, vars, drop = FALSE], na.rm = TRUE)
})
for (cat_name in names(threat_categories)) {
  cat_vars <- threat_categories[[cat_name]]
  cat_weight <- cat_weight_mat[, cat_name]
  for (v in cat_vars) {
    within_cat_weight_mat[, v] <- ifelse(cat_weight > 0, cw_mat[, v] / cat_weight, 0)
  }
}

within_cat_sums <- sapply(threat_categories, function(vars) {
  rowSums(within_cat_weight_mat[, vars, drop = FALSE], na.rm = TRUE)
})
active_cat <- cat_weight_mat > 0
if (any(abs(within_cat_sums[active_cat] - 1) > 1e-8)) {
  stop("Within-category threat weights do not sum to one for all active country-categories")
}
cat("  Within-category weights normalized for", sum(active_cat), "active country-categories\n")

# integer index: maps each row to its country's row in cw_mat
cidx <- match(d$country_rast, countries)
cat("  Country weights computed for", nrow(cw_mat), "countries\n")

# step 3: category indices (weighted average within category, country-specific weights)
# uses vector-based loop to avoid allocating full 131M-row matrices
for (cat_name in names(threat_categories)) {
  cat_vars <- threat_categories[[cat_name]]
  norm_cols <- paste0("norm_", cat_vars)

  if (length(cat_vars) == 1) {
    # single-threat category: the index is the normalized indicator itself.
    # No positive-weight gate: a zero country weight is handled in step 4,
    # where the category simply receives no weight.
    d[[cat_name]] <- d[[norm_cols]]
  } else {
    w_sum <- numeric(nrow(d))
    v_wsum <- numeric(nrow(d))
    for (j in seq_along(cat_vars)) {
      w_j <- within_cat_weight_mat[cidx, cat_vars[j]]
      v_j <- d[[norm_cols[j]]]
      na_mask <- is.na(v_j)
      w_j[na_mask] <- 0
      v_j[na_mask] <- 0
      w_sum <- w_sum + w_j
      v_wsum <- v_wsum + v_j * w_j
    }
    d[[cat_name]] <- ifelse(w_sum > 0, v_wsum / w_sum, NA)
    rm(w_sum, v_wsum, w_j, v_j, na_mask)
  }
}
gc()

# step 4: overall composite = country-weighted mean of the 7 category indices.
# The category weight is the sum of its member threats' presence proportions in
# the pixel's country, so categories carrying more prevalent threats weigh more:
#   threat_composite = sum_C W_Cc * I_C / sum_C W_Cc,  W_Cc = sum_{k in C} wbar_kc
# This is the index behind the published results. An equal-weight average over
# active categories was used between commits e261d72 and this one; it produces a
# different index (corr = 0.87, global ATT -22.9% instead of -31.5%) and was
# never propagated into merge.unmatched.fst. Do not reintroduce it without
# re-running the full downstream pipeline. See
# ideas/20260826_q75_threshold.md for the diagnosis.
# cat_weight_mat (n_countries x 7) was already built above for the
# within-category normalization
cat_names <- names(threat_categories)
W_sum <- numeric(nrow(d))
tc_wsum <- numeric(nrow(d))
for (j in seq_along(cat_names)) {
  w_j <- cat_weight_mat[cidx, j]
  tc_j <- d[[cat_names[j]]]
  na_mask <- is.na(tc_j)
  w_j[na_mask] <- 0
  tc_j[na_mask] <- 0
  W_sum <- W_sum + w_j
  tc_wsum <- tc_wsum + tc_j * w_j
}
d$threat_composite <- ifelse(W_sum > 0, tc_wsum / W_sum, NA)
rm(W_sum, tc_wsum, w_j, tc_j, na_mask)

cat("threat_composite summary:\n")
print(summary(d$threat_composite))
cat("Category index summaries:\n")
for (cat_name in names(threat_categories)) {
  cat(" ", cat_name, "- range:", round(min(d[[cat_name]], na.rm = TRUE), 4),
      "to", round(max(d[[cat_name]], na.rm = TRUE), 4), "\n")
}

# clean up norm_ columns (keep tc_ category columns and threat_composite)
norm_cols_all <- paste0("norm_", threat_list_raw)
d <- d[, !(names(d) %in% norm_cols_all)]
rm(norm_params, norm_params_list, cw_mat, within_cat_weight_mat, within_cat_sums,
   active_cat, cidx, cat_weight_mat)
gc()
}

##########################################
### any_threat (binary) ##################
##########################################
{
cat("\n=== Computing any_threat (binary presence/absence) ===\n")
# dep_thresholds already resolved above

# compute presence matrix
presence_mat <- sapply(threat_list_raw, function(v) {
  thresh <- dep_thresholds[v, "value"]
  dir    <- dep_thresholds[v, "dir"]
  x <- d[[v]]
  if (dir == "<") as.integer(!is.na(x) & x < thresh) else as.integer(!is.na(x) & x > thresh)
})
d$any_threat <- as.integer(rowSums(presence_mat, na.rm = TRUE) > 0)
cat("any_threat:", sum(d$any_threat, na.rm = TRUE), "of", nrow(d), "obs",
    sprintf("(%.1f%%)\n", 100 * mean(d$any_threat, na.rm = TRUE)))

rm(presence_mat)
gc()
}


##############
### saving ###
##############
# only save pixel_id + computed indices
out_cols <- c("pixel_id", names(threat_categories), "threat_composite", "any_threat")
write_fst(d[, out_cols], "data/store/pathreat.data.threat-indices.fst")
