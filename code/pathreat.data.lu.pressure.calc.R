##################################################
### pathreat.data.lu.pressure.calc.R #############
### Predict land-use pressure from covariates ####
##################################################
#
# Purpose: Create continuous land-use pressure indices at pixel level
#
# Output variables:
#   lu.pressure.est    - Model-based: predicted from threat outcomes using Poisson/OLS
#   lu.pressure.comp   - Composite: normalized mean of 5 covariates (pop, cropsuit, elev, slope, access)
#   lu.pressure.access - Simple: 1 - normalized accessibility (travel time)
#
# Note: All min-max normalization uses the full pixel distribution (treated + control)
#
# Method for lu.pressure.est:
#   For each threat outcome, fit Poisson/OLS on a RANDOM SAMPLE
#   of control pixels using fixest, then predict for all pixels
#
# Input:   data/store/pathreat.data.merge.fst
# Output:  data/store/pathreat.data.lu.pressure.fst
#          results/pathreat.newest.lu.pressure.est.html
#          data/store/pathreat.data.lu.pressure.models.Rds (intermediate)
#
# Note: Run AFTER pathreat.data.merge.R, BEFORE pathreat.data.merge.unmatched.R
#
# WORKFLOW:
#   1. Run STEPS 1-4 to fit models and generate HTML output
#   2. Review results in HTML file
#   3. Run STEP 5 to make predictions (memory intensive)
#   4. Run STEP 6 to save final output
#
library(dplyr)
library(fst)
library(texreg)
library(fixest)

# load config (provides mlist, deplist, paths)
source("code/pathreat.analysis.config.R")

# threat variables to model (14 raw threats from config — excludes any_threat
# and threat_composite which don't exist yet at this pipeline stage)
threat_vars <- threat_list_raw

# matching covariates (predictors)
covariates <- mlist

# SAMPLE PERCENT for model fitting (adjust based on memory)
SAMPLE_PERCENT <- 1  # percent of controls to sample

# CONTROL FLAGS - set to TRUE/FALSE to run specific steps
RUN_ESTIMATION <- TRUE   # Steps 1-4: Load data, fit models, output HTML
RUN_PREDICTION <- TRUE   # Steps 5-6: Make predictions and save output

#########################################
### STEP 1: Load data ####################
#########################################
if (RUN_ESTIMATION) {
cat("STEP 1: Loading data (only necessary columns)...\n")

# define columns to keep
keep_cols <- c(
  "pixel_id", "x", "y",
  "country_rast", "biome_raster",
  "treat",
  covariates,
  threat_vars
)

# load from merged data (has all variables including treat)
d <- read_fst("data/store/pathreat.data.merge.fst", columns = keep_cols)
cat("Loaded merge data:", nrow(d), "rows x", ncol(d), "columns\n")
gc()
}

##############################################
### STEP 2: Prepare data for modeling ########
##############################################
if (RUN_ESTIMATION) {
cat("\nSTEP 2: Preparing data for modeling...\n")

# check which covariates and threats are available
available_covs <- covariates[covariates %in% names(d)]
available_threats <- threat_vars[threat_vars %in% names(d)]
cat("Available covariates:", length(available_covs), "/", length(covariates), "\n")
cat("Available threats:", length(available_threats), "/", length(threat_vars), "\n")

# identify control pixels (treat == 0)
controls_idx <- which(d$treat == 0)
cat("Control pixels:", length(controls_idx), "\n")

# create complete cases indicator for covariates (include country_rast)
covariate_complete <- complete.cases(d[, c(available_covs, "country_rast")])
cat("Pixels with complete covariates:", sum(covariate_complete), "\n")

# controls with complete covariates (for model fitting)
controls_complete_idx <- which(d$treat == 0 & covariate_complete)
cat("Controls with complete covariates:", length(controls_complete_idx), "\n")

# RANDOM SAMPLE of controls for model fitting
set.seed(42)  # reproducibility
sample_size <- round(length(controls_complete_idx) * SAMPLE_PERCENT / 100)
sample_idx <- sample(controls_complete_idx, sample_size)
cat("Random sample for modeling:", format(length(sample_idx), big.mark = ","),
    "controls (", SAMPLE_PERCENT, "% )\n")

# create sample data frame for modeling
d_sample <- d[sample_idx, c(available_covs, "country_rast", available_threats)]
cat("Sample data frame created:", nrow(d_sample), "rows x", ncol(d_sample), "cols\n")

# identify countries in sample (needed for prediction)
countries_in_sample <- unique(d_sample$country_rast)
cat("Countries in sample:", length(countries_in_sample), "\n")

gc()
}

##############################################
### STEP 3: Fit global models (fixest) #######
##############################################
if (RUN_ESTIMATION) {
cat("\nSTEP 3: Fitting GLOBAL models (fixest, no country FE)...\n")

# model formulas for fixest (no fixed effects)
formula_feols <- as.formula(paste("threat_y ~", paste(available_covs, collapse = " + ")))
formula_fepois <- as.formula(paste("threat_y ~", paste(available_covs, collapse = " + ")))
cat("Model formula:", deparse(formula_feols), "\n\n")

# store fitted models
model_list <- list()

# store model diagnostics
model_diagnostics <- data.frame(
  threat = available_threats,
  n_obs = NA,
  n_nonzero = NA,
  r2 = NA,
  stringsAsFactors = FALSE
)

for (i in seq_along(available_threats)) {
  threat <- available_threats[i]
  cat("Modeling:", threat, "(", i, "/", length(available_threats), ")...\n")

  # add response variable to sample data
  d_sample$threat_y <- d_sample[[threat]]

  # check for negative values (swu can be negative)
  has_negative <- any(d_sample$threat_y < 0, na.rm = TRUE)

  if (has_negative) {
    cat("  - Has negative values, using feols (linear)\n")
    model <- tryCatch({
      feols(formula_feols, data = d_sample)
    }, error = function(e) {
      cat("  - feols failed:", e$message, "\n")
      NULL
    })
  } else {
    cat("  - Non-negative, using fepois (Poisson)\n")
    model <- tryCatch({
      fepois(formula_fepois, data = d_sample)
    }, error = function(e) {
      cat("  - fepois failed, trying feols:", e$message, "\n")
      tryCatch({
        feols(formula_feols, data = d_sample)
      }, error = function(e2) {
        cat("  - feols also failed:", e2$message, "\n")
        NULL
      })
    })
  }

  if (!is.null(model)) {
    # store model
    model_list[[threat]] <- model

    # diagnostics
    model_diagnostics$n_obs[i] <- model$nobs
    model_diagnostics$n_nonzero[i] <- sum(d_sample$threat_y > 0, na.rm = TRUE)
    model_diagnostics$r2[i] <- tryCatch(r2(model, type = "r2"), error = function(e) NA)

    cat("  - N obs:", model_diagnostics$n_obs[i],
        ", Non-zero:", model_diagnostics$n_nonzero[i],
        ", R2:", round(model_diagnostics$r2[i], 4), "\n")
  }

  d_sample$threat_y <- NULL
  gc()
}

cat("\nModel diagnostics summary:\n")
print(model_diagnostics)

# save models and metadata for later use
model_output <- list(
  models = model_list,
  diagnostics = model_diagnostics,
  countries_in_sample = countries_in_sample,
  available_covs = available_covs,
  available_threats = available_threats,
  sample_percent = SAMPLE_PERCENT
)
saveRDS(model_output, "data/store/pathreat.data.lu.pressure.models.Rds")
cat("\nModels saved to: data/store/pathreat.data.lu.pressure.models.Rds\n")

}

##############################################
### STEP 4: Output model results (texreg) ####
##############################################
if (RUN_ESTIMATION) {
cat("\nSTEP 4: Generating HTML regression tables...\n")

# generate HTML table using texreg
h <- htmlreg(
  model_list,
  custom.model.names = names(model_list),
  stars = c(0.01, 0.05, 0.1),
  digits = 4,
  table = FALSE,
  custom.note = paste0("Based on random sample of ", SAMPLE_PERCENT,
                       "% of control pixels.")
)

# header with method description
header_note <- paste0(
  "<h2>Land-Use Pressure Models</h2>",
  "<p><b>Method:</b> fixest fepois (Poisson) or feols (OLS for variables with negative values)</p>",
  "<p><b>Sample:</b> Random sample of ", SAMPLE_PERCENT, "% of control pixels</p>",
  "<p><b>Covariates:</b> ", paste(available_covs, collapse = ", "), "</p>",
  "<p><b>Fixed effects:</b> None</p>"
)

# save to HTML file
cat(file = paste0(paths$results_dir, "pathreat.newest.lu.pressure.est.html"),
    c(header_note, "<br>", h), append = FALSE)

cat("HTML table saved to:", paste0(paths$results_dir, "pathreat.newest.lu.pressure.est.html"), "\n")

# clean up sample data (after texreg uses it)
rm(d_sample)
gc()

cat("\n=== ESTIMATION COMPLETE ===\n")
cat("Review results, then set RUN_PREDICTION <- TRUE to continue.\n")
}

# helper function: normalize to [0,1] using min-max of full distribution
norm_minmax <- function(x) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  if (is.na(x_min) || is.na(x_max) || !is.finite(x_min) || !is.finite(x_max) || x_max <= x_min) {
    return(rep(NA_real_, length(x)))
  }
  (x - x_min) / (x_max - x_min)
}

##############################################
### STEP 5: Make predictions #################
##############################################
if (RUN_PREDICTION) {
cat("\nSTEP 5: Making predictions for all pixels...\n")

# load models if not in memory
if (!exists("model_list") || length(model_list) == 0) {
  cat("Loading saved models...\n")
  model_output <- readRDS("data/store/pathreat.data.lu.pressure.models.Rds")
  model_list <- model_output$models
  countries_in_sample <- model_output$countries_in_sample
  available_covs <- model_output$available_covs
  available_threats <- model_output$available_threats
  SAMPLE_PERCENT <- model_output$sample_percent
}

# load full data if not in memory
if (!exists("d") || !is.data.frame(d) || nrow(d) == 0) {
  cat("Loading full data from merge.fst...\n")
  keep_cols <- c(
    "pixel_id", "x", "y",
    "country_rast", "biome_raster",
    "treat",
    available_covs, available_threats
  )
  d <- read_fst("data/store/pathreat.data.merge.fst", columns = keep_cols)
  gc()
}

# always recreate indices (needed for prediction and simple pressure indices)
cat("Creating indices...\n")
controls_idx <- which(d$treat == 0)
covariate_complete <- complete.cases(d[, available_covs])
cat("Controls:", length(controls_idx), ", Complete cases:", sum(covariate_complete), "\n")

# pixels that can be predicted (no country restriction without FE)
can_predict <- covariate_complete
cat("Pixels that can be predicted:", sum(can_predict), "\n")

# store predictions incrementally
lu.pressure.est_sum <- rep(0, nrow(d))
lu.pressure.est_count <- rep(0L, nrow(d))

for (i in seq_along(available_threats)) {
  threat <- available_threats[i]

  if (is.null(model_list[[threat]])) {
    cat("Skipping", threat, "(no model)\n")
    next
  }

  cat("Predicting:", threat, "(", i, "/", length(available_threats), ")...\n")

  m <- model_list[[threat]]

  # extract coefficients (no fixed effects)
  coefs <- coef(m)

  # separate intercept from covariate coefficients
  intercept <- ifelse("(Intercept)" %in% names(coefs), coefs["(Intercept)"], 0)
  beta <- coefs[names(coefs) != "(Intercept)"]

  # compute linear predictor: intercept + X %*% beta
  X <- as.matrix(d[, names(beta), drop = FALSE])
  lin_pred <- intercept + X %*% beta

  # for Poisson models (fepois), predictions are exp(linear predictor)
  # for OLS models (feols), predictions are just linear predictor
  is_poisson <- inherits(m, "fixest") && !is.null(m$family) && m$family$family == "poisson"
  if (is.null(m$family)) {
    is_poisson <- "fml_all" %in% names(m) && grepl("fepois", deparse(m$call)[1])
  }

  pred_threat <- rep(NA_real_, nrow(d))
  if (is_poisson) {
    pred_threat[can_predict] <- exp(lin_pred[can_predict])
  } else {
    pred_threat[can_predict] <- lin_pred[can_predict]
  }

  rm(X, lin_pred)

  # normalize predictions (0-1 using min-max of full distribution)
  p_min <- min(pred_threat, na.rm = TRUE)
  p_max <- max(pred_threat, na.rm = TRUE)

  if (!is.na(p_min) && !is.na(p_max) && is.finite(p_min) && is.finite(p_max) && p_max > p_min) {
    norm_pred <- (pred_threat - p_min) / (p_max - p_min)

    # add to running sum
    valid_idx <- !is.na(norm_pred)
    lu.pressure.est_sum[valid_idx] <- lu.pressure.est_sum[valid_idx] + norm_pred[valid_idx]
    lu.pressure.est_count[valid_idx] <- lu.pressure.est_count[valid_idx] + 1L

    rm(norm_pred)
  }

  rm(pred_threat)
  gc()
}

# compute final lu.pressure.est as mean of normalized predictions
d$lu.pressure.est <- ifelse(lu.pressure.est_count > 0,
                               lu.pressure.est_sum / lu.pressure.est_count,
                               NA_real_)
rm(lu.pressure.est_sum, lu.pressure.est_count)
gc()

cat("\nlu.pressure.est summary:\n")
print(summary(d$lu.pressure.est))
}

##############################################
### STEP 5b: Simple pressure indices #########
##############################################
if (RUN_PREDICTION) {
cat("\nSTEP 5b: Computing simple pressure indices...\n")

# --- lu.pressure.access: accessibility-based (simple) ---
# Higher access value = longer travel time = more remote = LESS pressure
# So invert: lu.pressure.access = 1 - norm(access)
cat("Computing lu.pressure.access...\n")
d$lu.pressure.access <- 1 - norm_minmax(d$access)
cat("lu.pressure.access summary:\n")
print(summary(d$lu.pressure.access))

# --- lu.pressure.comp: composite index (5 covariates) ---
# Pressure-increasing: pop_2000, cropsuit2 (higher = more pressure)
# Pressure-decreasing: elevation, slope, access (higher = less pressure, so invert)
cat("\nComputing lu.pressure.comp...\n")

# normalize each covariate
norm_pop <- norm_minmax(d$pop_2000)
norm_cropsuit <- norm_minmax(d$cropsuit2)
norm_elev <- norm_minmax(d$elevation)
norm_slope <- norm_minmax(d$slope)
norm_access <- norm_minmax(d$access)

# composite: mean of (pressure-increasing) and (1 - pressure-decreasing)
d$lu.pressure.comp <- rowMeans(cbind(
  norm_pop,
  norm_cropsuit,
  1 - norm_elev,
  1 - norm_slope,
  1 - norm_access
), na.rm = TRUE)

# set to NA where all components are NA
all_na <- is.na(norm_pop) & is.na(norm_cropsuit) & is.na(norm_elev) &
          is.na(norm_slope) & is.na(norm_access)
d$lu.pressure.comp[all_na] <- NA_real_

rm(norm_pop, norm_cropsuit, norm_elev, norm_slope, norm_access, all_na)
gc()

cat("lu.pressure.comp summary:\n")
print(summary(d$lu.pressure.comp))
}

##############################################
### STEP 6: Save data output #################
##############################################
if (RUN_PREDICTION) {
cat("\nSTEP 6: Saving data output...\n")

# keep only essential columns for downstream analysis
output_cols <- c(
  "pixel_id", "x", "y",
  "treat", "country_rast", "biome_raster",
  "lu.pressure.est",
  "lu.pressure.comp",
  "lu.pressure.access"
)

# check which columns exist
output_cols <- output_cols[output_cols %in% names(d)]
d_out <- d[, output_cols]

cat("Output dimensions:", dim(d_out), "\n")
cat("Output columns:", paste(names(d_out), collapse = ", "), "\n")

# save data
write_fst(d_out, "data/store/pathreat.data.lu.pressure.fst", compress = 50)
cat("Saved to: data/store/pathreat.data.lu.pressure.fst\n")

cat("\nFile size (MB):", round(file.info("data/store/pathreat.data.lu.pressure.fst")$size / 10^6, 1), "\n")
}

cat("\n=== DONE ===\n")
cat("Outputs:\n")
cat("  - Models: data/store/pathreat.data.lu.pressure.models.Rds\n")
cat("  - HTML table:", paste0(paths$results_dir, "pathreat.newest.lu.pressure.est.html"), "\n")
cat("  - Data: data/store/pathreat.data.lu.pressure.fst\n")
