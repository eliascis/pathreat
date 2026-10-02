##############################################
### pathreat.analysis.robustness.est.R #####
### Robustness: sequential FE & controls #####
### + random (unmatched) sample ##############
##############################################

library(fixest)
library(texreg)
library(fst)

source("code/pathreat.analysis.config.R")

######################
### load data ########
######################
{
cat("Loading matched data...\n")
if (!"d.mbase" %in% ls()) {
  d.mbase <- load_matched_data()
}
cat(sprintf("  Loaded %d matched observations\n", nrow(d.mbase)))
}

star_fn <- function(p) {
  if (p < 0.01) "\\sym{***}"
  else if (p < 0.05) "\\sym{**}"
  else if (p < 0.1) "\\sym{*}"
  else ""
}

############################################
### 1. Unmatched sample: no FE #############
############################################
{
cat("\n=== Block 1: Unmatched sample, no FE ===\n")

# load unmatched data from main merge file (has threat_composite)
needed_cols <- c("treat", "country_rast", "biome_raster", "threat_composite")

cat("  Loading unmatched data (selected columns)...\n")
d_unmatched <- read_fst(paths$data_unmatched, columns = needed_cols)
cat(sprintf("  Loaded %d unmatched observations\n", nrow(d_unmatched)))

cat(sprintf("  Unmatched observations: %d\n", nrow(d_unmatched)))

# threat_composite already computed in data.1

ctrl_mean_unmatched <- mean(d_unmatched$threat_composite[d_unmatched$treat == 0], na.rm = TRUE)

# Model 1: Unmatched, no FE, no controls
cat("  Model 1: Unmatched, no FE\n")
rob1 <- feols(threat_composite ~ treat, data = d_unmatched, cluster = "country_rast")

# Model 1b: Unmatched, Country + Biome FE
cat("  Model 1b: Unmatched, Country + Biome FE\n")
rob1b <- feols(threat_composite ~ treat | country_rast + biome_raster,
               data = d_unmatched, cluster = "country_rast")

rm(d_unmatched); gc()
}

############################################
### 2. Matched sample: threat_composite ##########
### Sequential addition of FE & controls ###
############################################
{
cat("\n=== Block 2: Robustness regressions on matched sample (threat_composite) ===\n")

d <- d.mbase

# apply the canonical complete-pair sample for threat_composite
pair_index <- build_matched_pair_index(d)
sample_mask <- matched_pair_sample_mask(d, "threat_composite", pair_index)
robust_columns <- unique(c(
  "matched_pair_id",
  "treat",
  "country_rast",
  "biome_raster",
  "threat_composite",
  "wdpaid",
  mlist
))
d_rob <- d[sample_mask, robust_columns, drop = FALSE]
n_rob_treated <- sum(d_rob$treat == 1)
n_rob_control <- sum(d_rob$treat == 0)
if (n_rob_treated != n_rob_control || nrow(d_rob) != 2L * n_rob_treated) {
  stop("Matched robustness sample is not pair-complete")
}

clustvar <- c("country_rast")

# Model 2: Matched, no FE, no controls
cat("  Model 2: Matched, no FE, no controls\n")
rob2 <- feols(threat_composite ~ treat, data = d_rob, cluster = clustvar)

# Model 3: Matched, Country FE
cat("  Model 3: Matched, Country FE\n")
rob3 <- feols(threat_composite ~ treat | country_rast, data = d_rob, cluster = clustvar)

# Model 4: Matched, Country + Biome FE
cat("  Model 4: Matched, Country + Biome FE\n")
rob4 <- feols(threat_composite ~ treat | country_rast + biome_raster, data = d_rob, cluster = clustvar)

# Model 5: Matched, Country + Biome FE + matching covariates
cat("  Model 5: Matched, Country + Biome FE + matching covariates\n")
controls_formula <- paste(mlist, collapse = " + ")
rob5 <- feols(as.formula(paste("threat_composite ~ treat +", controls_formula,
                                "| country_rast + biome_raster")),
              data = d_rob, cluster = clustvar)

ctrl_mean_matched <- mean(d_rob$threat_composite[d_rob$treat == 0])

# Create PA fixed-effect identifier: wdpaid for treated pixels,
# matched PA's wdpaid for control pixels (via matched_pair_id)
pa_lookup <- d_rob[d_rob$treat == 1, c("matched_pair_id", "wdpaid")]
d_rob$pa_fe_id <- d_rob$wdpaid
ctrl_idx <- d_rob$treat == 0
d_rob$pa_fe_id[ctrl_idx] <- pa_lookup$wdpaid[match(d_rob$matched_pair_id[ctrl_idx],
                                                     pa_lookup$matched_pair_id)]
cat(sprintf("  PA FE: %d unique PAs (controls assigned matched PA's wdpaid)\n",
            length(unique(d_rob$pa_fe_id[d_rob$treat == 1]))))
rm(pa_lookup, ctrl_idx)

# Model 6: Matched, Country + PA FE (no biome)
cat("  Model 6: Matched, Country + PA FE\n")
rob6 <- feols(threat_composite ~ treat | country_rast + pa_fe_id,
              data = d_rob, cluster = clustvar)

# Model 8: Matched, Country FE, Fractional logit (quasi-binomial)
cat("  Model 8: Matched, Country FE, Fractional logit\n")
rob8 <- feglm(threat_composite ~ treat | country_rast, data = d_rob,
              family = quasibinomial(link = "logit"), cluster = clustvar)

matched_model_nobs <- vapply(
  list(rob2, rob3, rob4, rob5, rob6, rob8),
  function(model) model$nobs,
  numeric(1)
)
if (any(matched_model_nobs != nrow(d_rob))) {
  stop("At least one matched robustness model changed the canonical sample")
}

rm(d, d_rob, pair_index, sample_mask); gc()
}

############################################
### 4. Conley spatial SEs ##################
### (10% subsample for feasibility) ########
### Note: the 20%/10% matched-sample      ##
### subsample placebos previously here     ##
### (Models 7, 8) moved to the randomrob  ##
### pipeline (Figure A.3 / Table A.1).     ##
############################################
{
cat("\n=== Block 4: Conley spatial SEs (10% subsample) ===\n")

library(sf)

d <- d.mbase

# apply the same canonical complete-pair outcome sample
pair_index <- build_matched_pair_index(d)
sample_mask <- matched_pair_sample_mask(d, "threat_composite", pair_index)
d <- d[
  sample_mask,
  c("matched_pair_id", "treat", "country_rast", "biome_raster",
    "threat_composite", "x", "y"),
  drop = FALSE
]
pair_index <- build_matched_pair_index(d)

# require usable coordinates on both pair members before pair sampling
coord_complete <- is.finite(d$x[pair_index$treated_rows]) &
  is.finite(d$y[pair_index$treated_rows]) &
  is.finite(d$x[pair_index$control_rows]) &
  is.finite(d$y[pair_index$control_rows])
eligible_conley_pairs <- pair_index$pair_ids[coord_complete]
if (length(eligible_conley_pairs) == 0) {
  stop("No coordinate-complete matched pairs are available for Conley SEs")
}

# draw 10% of coordinate-complete matched pairs
set.seed(123)
n_conley_pairs <- max(1L, round(length(eligible_conley_pairs) * 0.10))
sample_pair_ids <- sample(eligible_conley_pairs, n_conley_pairs)
d_conley <- d[d$matched_pair_id %in% sample_pair_ids, , drop = FALSE]
n_conley_treated <- sum(d_conley$treat == 1)
n_conley_control <- sum(d_conley$treat == 0)
if (n_conley_treated != n_conley_pairs ||
    n_conley_control != n_conley_pairs ||
    nrow(d_conley) != 2L * n_conley_pairs) {
  stop("Conley subsample is not pair-complete")
}
cat(sprintf("  Drawing 10%% subsample: %s of %s coordinate-complete pairs (%s rows)\n",
            formatC(n_conley_pairs, format = "d", big.mark = ","),
            formatC(length(eligible_conley_pairs), format = "d", big.mark = ","),
            formatC(nrow(d_conley), format = "d", big.mark = ",")))

# reproject Mollweide x/y to WGS84 lat/lon
cat("  Reprojecting coordinates to lat/lon...\n")
pts <- st_as_sf(d_conley, coords = c("x", "y"),
                crs = "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m")
pts <- st_transform(pts, crs = 4326)
xy_ll <- st_coordinates(pts)
d_conley$lon <- xy_ll[, 1]
d_conley$lat <- xy_ll[, 2]
rm(pts, xy_ll); gc()

ctrl_mean_conley <- mean(d_conley$threat_composite[d_conley$treat == 0])

# Model 7: 10% subsample, Country + Biome FE, Conley SEs (100 km cutoff, 10 km pixel aggregation)
cat("  Model 7: 10% subsample, Country + Biome FE, Conley SEs (100 km)\n")
rob7 <- feols(threat_composite ~ treat | country_rast + biome_raster,
              data = d_conley,
              vcov = conley(cutoff = 100, pixel = 10))
if (rob7$nobs != nrow(d_conley)) {
  stop("Conley model changed the pair-complete subsample")
}

rm(d, d_conley, pair_index, sample_mask, eligible_conley_pairs,
   sample_pair_ids, coord_complete)
gc()
}

############################################
### 5. Combined outputs ####################
############################################
{
cat("\n=== Block 5: Saving outputs ===\n")

rob_list <- list(rob1, rob1b, rob2, rob3, rob4, rob5, rob6, rob7, rob8)
rob_names <- c("Unmatched", "Unmatched + FE", "No FE", "Country FE",
               "Country + Biome FE", "Country + Biome FE + Controls",
               "Country + PA FE", "Conley SEs", "Frac. logit")

# --- extract results for Rds ---
ctrl_means_all <- c(rep(ctrl_mean_unmatched, 2), rep(ctrl_mean_matched, 5),
                    ctrl_mean_conley, ctrl_mean_matched)

rob_results <- do.call(rbind, lapply(seq_along(rob_list), function(i) {
  e <- rob_list[[i]]
  ct <- coeftable(e)["treat", , drop = FALSE]
  estimate <- regression_estimate(e)
  data.frame(
    variable = "threat_composite",
    estimate_type = paste0("robustness_", i),
    model_label = rob_names[i],
    coef = ct[1, 1],
    se = ct[1, 2],
    t_stat = ct[1, 3],
    pval = ct[1, 4],
    ci_low = estimate$ci_low,
    ci_high = estimate$ci_high,
    n_obs = e$nobs,
    control_mean = ctrl_means_all[i],
    row.names = NULL
  )
}))

# --- HTML output ---
omit_pattern <- paste(c("\\(Intercept\\)", mlist), collapse = "|")

h_rob <- htmlreg(
  rob_list,
  omit.coef = omit_pattern,
  custom.model.names = rob_names,
  custom.gof.rows = list(
    "Control Mean"        = round(ctrl_means_all, 3),
    "Sample"              = c(rep("Unmatched", 2), rep("Full matched", 5),
                              "10\\% subsample", "Full matched"),
    "Estimator"           = c(rep("OLS", 8), "Frac. logit"),
    "Country FE"          = c("No", "Yes", "No", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes"),
    "Biome FE"            = c("No", "Yes", "No", "No", "Yes", "Yes", "No", "Yes", "No"),
    "PA FE"               = c("No", "No", "No", "No", "No", "No", "Yes", "No", "No"),
    "Matching covariates" = c("No", "No", "No", "No", "No", "Yes", "No", "No", "No"),
    "SE type"             = c(rep("Country", 7), "Conley (100 km)", "Country")
  ),
  stars = c(0.01, 0.05, 0.1),
  digits = 3,
  table = FALSE
)

# append to newest.html
cat(file = paste0(paths$results_dir, "pathreat.est.newest.html"),
    c("<br><br><b>Robustness: threat_composite — Sequential addition of FE and controls</b>", h_rob),
    append = TRUE)

# write separate robustness HTML
cat(file = paste0(paths$results_dir, "pathreat.robustness.est.html"),
    c("<b>Robustness: threat_composite — Sequential addition of FE and controls</b>", h_rob),
    append = FALSE)

# --- LaTeX table: pub/tables/tab.rob.controls.tex ---
cat("  Generating LaTeX table for robustness regressions\n")

tex <- character()

# header rows
tex <- c(tex, "              & {(1)}        & {(2)}              & {(3)}    & {(4)}         & {(5)}              & {(6)}         & {(7)}       & {(8)}              & {(9)}             \\\\")
tex <- c(tex, "              & {Unmatched}  & {Unmatched + FE}   & {No FE}  & {Country FE}  & {Country + Biome}  & {+ Controls}  & {+ PA FE}  & {Conley SEs}       & {Frac. logit}     \\\\")
tex <- c(tex, "\\midrule")

# coefficient row
coefs <- sapply(rob_list, function(e) coeftable(e)["treat", 1])
ses   <- sapply(rob_list, function(e) coeftable(e)["treat", 2])
pvals <- sapply(rob_list, function(e) coeftable(e)["treat", 4])

coef_cells <- paste0(sprintf("%.3f", coefs), sapply(pvals, star_fn))
se_cells   <- paste0("(", sprintf("%.3f", ses), ")")

tex <- c(tex, paste("Protected area", paste0(" & ", coef_cells, collapse = ""), "\\\\"))
tex <- c(tex, paste("             ", paste0(" & ", se_cells, collapse = ""), "\\\\"))

# footer rows
tex <- c(tex, "\\midrule")
tex <- c(tex, paste("Control mean",
  paste0(" & ", sprintf("%.3f", ctrl_means_all), collapse = ""), "\\\\"))
tex <- c(tex, paste("Num.\\ obs.",
  paste0(" & {", formatC(sapply(rob_list, function(e) e$nobs), format = "d", big.mark = ","), "}", collapse = ""), "\\\\"))
tex <- c(tex, paste("Sample",
  paste0(" & {", c(rep("Unmatched", 2), rep("Full matched", 5),
                    "10\\% subsample", "Full matched"), "}", collapse = ""), "\\\\"))
tex <- c(tex, paste("Estimator",
  paste0(" & {", c(rep("OLS", 8), "Frac. logit"), "}", collapse = ""), "\\\\"))
tex <- c(tex, paste("Country FE",
  paste0(" & {", c("No", "Yes", "No", "Yes", "Yes", "Yes", "Yes", "Yes", "Yes"), "}", collapse = ""), "\\\\"))
tex <- c(tex, paste("Biome FE",
  paste0(" & {", c("No", "Yes", "No", "No", "Yes", "Yes", "No", "Yes", "No"), "}", collapse = ""), "\\\\"))
tex <- c(tex, paste("PA FE",
  paste0(" & {", c("No", "No", "No", "No", "No", "No", "Yes", "No", "No"), "}", collapse = ""), "\\\\"))
tex <- c(tex, paste("Matching covariates",
  paste0(" & {", c("No", "No", "No", "No", "No", "Yes", "No", "No", "No"), "}", collapse = ""), "\\\\"))
tex <- c(tex, paste("SE type",
  paste0(" & {", c(rep("Country", 7), "Conley (100 km)", "Country"), "}", collapse = "")))

# strip trailing \\ from last line (safety)
tex[length(tex)] <- sub(" *\\\\\\\\$", "", tex[length(tex)])

writeLines(tex, "pub/tables/tab.rob.controls.tex")
cat("  Saved: pub/tables/tab.rob.controls.tex\n")

# --- save Rds ---
saveRDS(rob_results, paste0(paths$results_dir, "pathreat.robustness.est.Rds"))
cat("  Saved: results/pathreat.robustness.est.Rds\n")

rm(rob1, rob1b, rob2, rob3, rob4, rob5, rob6, rob7, rob8, rob_list, h_rob, rob_results); gc()
}

cat("\n=== Robustness estimation complete ===\n")
cat("Output HTML: ", paths$results_dir, "pathreat.robustness.est.html\n")
cat("Output Rds:  ", paths$results_dir, "pathreat.robustness.est.Rds\n")
cat("Output LaTeX:", "pub/tables/tab.rob.controls.tex\n")
