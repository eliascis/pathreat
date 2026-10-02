##############################################
### pathreat.analysis.sfa.est.R #############
### Stochastic Frontier Analysis ############
### PA-level efficiency gaps ################
###                                          #
### Main spec: cost frontier on the within- #
### PA threat level (tc_protected = delta + #
### tc_control) with PA size in the         #
### inefficiency variance (uhet ~ log_size).#
### See ideas/20260427_sfa-delta-vs-levels  #
### .md for the rationale.                  #
##############################################

library(dplyr)
library(ggplot2)
library(texreg)
library(patchwork)
library(sfaR)

source("code/pathreat.analysis.config.R")

pa_country_membership_path <-
  "data/store/pathreat.analysis.byPA-country-membership.Rds"

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

######################
### helpers ##########
######################

star_fn <- function(p) {
  if (p < 0.01) "\\sym{***}"
  else if (p < 0.05) "\\sym{**}"
  else if (p < 0.1) "\\sym{*}"
  else ""
}

# check if SFA model produced valid estimates (no hessian inversion failure)
sfa_valid <- function(m) {
  if (is.null(m)) return(FALSE)
  s <- tryCatch(summary(m)$mlRes, error = function(e) NULL)
  if (is.null(s)) return(FALSE)
  if (any(is.infinite(s[, "Std. Error"]))) return(FALSE)
  ll <- as.numeric(logLik(m))
  if (abs(ll) > 1e6) return(FALSE)
  TRUE
}

# safe efficiency extraction (sfaR bugs out for intercept-only models)
safe_efficiencies <- function(m) {
  if (is.null(m)) return(NULL)
  tryCatch(efficiencies(m), error = function(e) NULL)
}

safe_nobs <- function(m) {
  if (is.null(m)) return(NA_integer_)
  tryCatch(nobs(m), error = function(e) NA_integer_)
}

# robust sfacross fit: BFGS first, fall back to Nelder-Mead on NA-gradient
# failures (matching covariates with mixed scales sometimes trigger this).
fit_sfa <- function(formula, udist = "hnormal", uhet = NULL, data = d_pa) {
  call_args <- list(
    formula     = formula,
    data        = data,
    udist       = udist,
    S           = -1L,
    logDepVar   = FALSE,
    hessianType = 2L,
    itermax     = 4000,
    printInfo   = FALSE
  )
  if (!is.null(uhet)) call_args$uhet <- uhet

  m <- tryCatch(
    do.call(sfacross, c(call_args, list(method = "bfgs"))),
    error = function(e) { cat("    !! BFGS failed:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(m)) {
    cat("    Retrying with Nelder-Mead (no gradient required)...\n")
    m <- tryCatch(
      do.call(sfacross, c(call_args, list(method = "nm"))),
      error = function(e) { cat("    !! NM also failed:", conditionMessage(e), "\n"); NULL }
    )
  }
  m
}

conv_str <- function(m) if (is.null(m)) "FAILED" else as.character(m$convergence)

#######################################
### 1. Pair data + PA aggregation #####
#######################################
{
cat("\n=== Block 1: Pair data and PA-level aggregation ===\n")

# pressure variable selection (one of: lu.pressure.est, lu.pressure.comp,
# lu.pressure.access, lu.vulnerable)
# lu.pressure.comp and lu.pressure.access are collinear with matching covariates
# at PA level (access ↔ lu.pressure.access: r = -1; cropsuit2 ↔ lu.pressure.comp: r = 0.95)
pressure_var <- "lu.pressure.est"
frontier_covs <- c(mlist, pressure_var)

# use the canonical complete-pair sample for threat_composite
pair_index <- build_matched_pair_index(d.mbase)
sample_mask <- matched_pair_sample_mask(
  d.mbase,
  "threat_composite",
  pair_index
)
d <- d.mbase[sample_mask, ]
cat(sprintf("  Pair-complete outcome sample: %s observations\n",
            format(nrow(d), big.mark = ",")))

# one row per pair: treat == 1 has delta + all covariates
keep_cols <- c("matched_pair_id", "delta", frontier_covs, pa_cols)
d_pairs <- d[d$treat == 1, keep_cols]
d_pairs$log_size <- log(d_pairs$size + 1)

# control-side threat_composite for frontier covariate
d_ctrl <- d[d$treat == 0, c("matched_pair_id", "threat_composite")]
names(d_ctrl)[2] <- "tc_control"
d_pairs <- merge(d_pairs, d_ctrl, by = "matched_pair_id", all.x = TRUE)
rm(d_ctrl)

# treated-side threat composite (the within-PA level) — the new DV
d_treat <- d[d$treat == 1, c("matched_pair_id", "threat_composite")]
names(d_treat)[2] <- "tc_protected_pair"
d_pairs <- merge(d_pairs, d_treat, by = "matched_pair_id", all.x = TRUE)
rm(d_treat)

if (!file.exists(pa_country_membership_path)) {
  stop("Missing PA-country membership audit: ", pa_country_membership_path)
}
d_pa_country <- readRDS(pa_country_membership_path) %>%
  arrange(wdpaid, country_rast)
membership_required <- c("wdpaid", "country_rast", "n_pixels_country")
if (!is.data.frame(d_pa_country) ||
    !all(membership_required %in% names(d_pa_country)) ||
    anyNA(d_pa_country[, membership_required, drop = FALSE]) ||
    anyDuplicated(d_pa_country[c("wdpaid", "country_rast")]) ||
    any(d_pa_country$n_pixels_country <= 0)) {
  stop("Stored PA-country membership audit is malformed")
}

d_pa_country_current <- d_pairs %>%
  group_by(wdpaid, country_rast) %>%
  summarise(
    n_pixels_country = n(),
    .groups = "drop"
  ) %>%
  arrange(wdpaid, country_rast)
if (!identical(
      as.data.frame(d_pa_country),
      as.data.frame(d_pa_country_current)
    )) {
  stop("Stored PA-country membership audit is stale")
}

cat(sprintf("  Pairs: %s\n", format(nrow(d_pairs), big.mark = ",")))
cat(sprintf("  delta:        mean = %.4f, sd = %.4f\n",
            mean(d_pairs$delta, na.rm = TRUE), sd(d_pairs$delta, na.rm = TRUE)))
cat(sprintf("  tc_control:   mean = %.4f, NAs = %d\n",
            mean(d_pairs$tc_control, na.rm = TRUE), sum(is.na(d_pairs$tc_control))))
cat(sprintf("  tc_protected: mean = %.4f (= delta + tc_control = %.4f)\n",
            mean(d_pairs$tc_protected_pair, na.rm = TRUE),
            mean(d_pairs$delta + d_pairs$tc_control, na.rm = TRUE)))

rm(d); gc()

# PA-level aggregation
d_pa <- d_pairs %>%
  group_by(wdpaid) %>%
  summarise(
    observed_delta = mean(delta, na.rm = TRUE),
    observed_protected = mean(tc_protected_pair, na.rm = TRUE),
    observed_control = mean(tc_control, na.rm = TRUE),
    across(all_of(frontier_covs), ~ mean(.x, na.rm = TRUE)),
    log_size     = first(log_size),
    size         = first(size),
    iucn_class   = first(iucn_class),
    size_class   = first(size_class),
    n_countries  = n_distinct(country_rast),
    country_rast = min(country_rast),
    biome_raster = first(biome_raster),
    n_pixels     = n(),
    .groups = "drop"
  ) %>%
  data.frame()

pa_estimates <- readRDS("data/store/pathreat.analysis.byPA.est.Rds")
required_pa <- c("wdpaid", "coef", "adjusted_control", "adjusted_protected", "status")
if (!all(required_pa %in% names(pa_estimates)) ||
    anyDuplicated(pa_estimates$wdpaid) ||
    !setequal(d_pa$wdpaid, pa_estimates$wdpaid)) {
  stop("Canonical PA regressions or standardized levels are unavailable")
}
pa_estimates <- pa_estimates[match(d_pa$wdpaid, pa_estimates$wdpaid), ]
if (any(pa_estimates$n_pixels != d_pa$n_pixels)) stop("PA sample counts changed")
d_pa$delta <- pa_estimates$coef
d_pa$tc_control <- pa_estimates$adjusted_control
d_pa$tc_protected <- pa_estimates$adjusted_protected
d_pa$regression_status <- pa_estimates$status
if (any(!is.finite(c(d_pa$delta, d_pa$tc_control, d_pa$tc_protected))) ||
    max(abs(d_pa$delta - d_pa$observed_delta),
        abs(d_pa$tc_control - d_pa$observed_control),
        abs(d_pa$tc_protected - d_pa$observed_protected)) > 1e-8) {
  stop("PA regression inputs disagree with the complete-pair audit")
}

cat(sprintf("  PAs: %d\n", nrow(d_pa)))
cat(sprintf("  PA-level tc_protected: mean = %.4f, sd = %.4f, range = [%.4f, %.4f]\n",
            mean(d_pa$tc_protected), sd(d_pa$tc_protected),
            min(d_pa$tc_protected), max(d_pa$tc_protected)))

# drop PAs with NA in any frontier covariate
na_before <- nrow(d_pa)
d_pa <- d_pa[complete.cases(d_pa[, c("tc_protected", "tc_control", frontier_covs, "log_size")]), ]
if (nrow(d_pa) < na_before) {
  cat(sprintf("  Dropped %d PAs with NA in frontier covariates (%d remaining)\n",
              na_before - nrow(d_pa), nrow(d_pa)))
}

# quadratic term used by several specifications
d_pa$tc_control_sq <- d_pa$tc_control^2

rm(d_pairs); gc()
}


##############################
### 2. SFA estimation ########
##############################
{
cat("\n=== Block 2: SFA estimation (PA-level, levels DV: tc_protected) ===\n")

# Half-normal cost frontier (S = -1) with tc_protected as the dependent variable
# in every spec, all without an intercept. Models 1-3 progressively saturate the
# deterministic frontier; Model 4 adds the eight matching covariates (kept for
# diagnostics, NOT shown in Table C.8); Model 5 (MAIN) parameterizes the
# inefficiency variance via uhet ~ log_size. See ideas/20260427_sfa-delta-vs-
# levels.md for rationale.

# matching covariates (mlist) — used in model 4 only
match_rhs <- paste(mlist, collapse = " + ")

cat("  Model 1: tc_protected ~ 0 + tc_control\n")
sf_ctrl  <- fit_sfa(tc_protected ~ 0 + tc_control)
cat("    Convergence:", conv_str(sf_ctrl), "\n")

cat("  Model 2: + tc_control_sq\n")
sf_quad  <- fit_sfa(tc_protected ~ 0 + tc_control + tc_control_sq)
cat("    Convergence:", conv_str(sf_quad), "\n")

cat("  Model 3: + log_size\n")
sf_size  <- fit_sfa(tc_protected ~ 0 + tc_control + tc_control_sq + log_size)
cat("    Convergence:", conv_str(sf_size), "\n")

cat("  Model 4 [diagnostic, NOT in Table C.8]: + matching covariates\n")
sf_match <- fit_sfa(as.formula(paste(
  "tc_protected ~ 0 + tc_control + tc_control_sq + log_size +", match_rhs
)))
cat("    Convergence:", conv_str(sf_match), "\n")

cat("  Model 5 [MAIN]: tc_protected ~ 0 + tc_control + tc_control_sq, uhet ~ log_size\n")
sf_uhet  <- fit_sfa(
  formula = tc_protected ~ 0 + tc_control + tc_control_sq,
  uhet    = ~ log_size
)
cat("    Convergence:", conv_str(sf_uhet), "\n")

sf_main <- sf_uhet
}


##############################
### 3. Efficiency extraction #
##############################
{
cat("\n=== Block 3: Efficiency extraction (from sf_main = uhet) ===\n")

# list of all fitted models (Model 4 is kept for diagnostics but NOT shown
# in Table C.8 — see Block 5)
model_list <- list(
  ctrl  = sf_ctrl,   # Table C.8 col 1
  quad  = sf_quad,   # Table C.8 col 2
  size  = sf_size,   # Table C.8 col 3
  match = sf_match,  # diagnostic only, not in Table C.8
  uhet  = sf_uhet    # Table C.8 col 4 (MAIN)
)

# extract from main model (col 7: uhet)
eff <- efficiencies(sf_main)
d_pa$u_hat          <- eff$u
d_pa$u_lb           <- eff$uLB
d_pa$u_ub           <- eff$uUB
d_pa$frontier_level <- fitted(sf_main)                    # NEW: level-space frontier (>= 0)
d_pa$frontier       <- d_pa$frontier_level - d_pa$tc_control  # back-compat: delta-space frontier (<= 0)
d_pa$te_jlms        <- exp(-d_pa$u_hat)

cat(sprintf("  Mean efficiency (te_jlms): %.3f\n", mean(d_pa$te_jlms)))
cat(sprintf("  Median efficiency:         %.3f\n", median(d_pa$te_jlms)))
cat(sprintf("  Mean gap (u_hat):          %.4f\n", mean(d_pa$u_hat)))
cat(sprintf("  Median gap:                %.4f\n", median(d_pa$u_hat)))
cat(sprintf("  Mean frontier_level:       %.4f (range [%.4f, %.4f])\n",
            mean(d_pa$frontier_level),
            min(d_pa$frontier_level), max(d_pa$frontier_level)))
cat(sprintf("  Mean frontier (delta-eq):  %.4f\n", mean(d_pa$frontier)))

# variance parameters of the main spec (averaged across PAs since uhet)
pars <- coef(sf_main, extraPar = TRUE)
ml   <- summary(sf_main)$mlRes
gamma_val <- unname(pars["gamma"])  # back-compat with downstream figure annotations
cat(sprintf("  sigma_u (avg): %.4f\n", sqrt(unname(pars["sigmauSq"]))))
cat(sprintf("  sigma_v (avg): %.4f\n", sqrt(unname(pars["sigmavSq"]))))
cat(sprintf("  gamma   (avg): %.3f\n", gamma_val))
cat(sprintf("  Log-likelihood: %.1f\n", as.numeric(logLik(sf_main))))

if ("Zu_log_size" %in% rownames(ml)) {
  zu_int   <- ml["Zu_(Intercept)", "Coefficient"]
  zu_slope <- ml["Zu_log_size",    "Coefficient"]
  ls_q <- quantile(d_pa$log_size, c(0.05, 0.50, 0.95), na.rm = TRUE)
  sigma_u_at <- function(ls) sqrt(exp(zu_int + zu_slope * ls))
  cat(sprintf("  uhet: Zu_log_size = %.4f (z = %.1f)\n",
              zu_slope, ml["Zu_log_size", "z value"]))
  cat(sprintf("    sigma_u at log_size P05 (%.2f, PA = %5.0f km^2): %.4f\n",
              ls_q[1], exp(ls_q[1]) - 1, sigma_u_at(ls_q[1])))
  cat(sprintf("    sigma_u at log_size P50 (%.2f, PA = %5.0f km^2): %.4f\n",
              ls_q[2], exp(ls_q[2]) - 1, sigma_u_at(ls_q[2])))
  cat(sprintf("    sigma_u at log_size P95 (%.2f, PA = %5.0f km^2): %.4f\n",
              ls_q[3], exp(ls_q[3]) - 1, sigma_u_at(ls_q[3])))
}
}


#####################################
### 4. Country-level aggregation ####
#####################################
{
cat("\n=== Block 4: Country-level aggregation ===\n")

d_country_input <- d_pa_country %>%
  inner_join(
    d_pa %>%
      select(-country_rast),
    by = "wdpaid"
  ) %>%
  mutate(
    country_pa_area = size * n_pixels_country / n_pixels
  )

d_country_pa_check <- d_country_input %>%
  group_by(wdpaid) %>%
  summarise(
    n_pixels_country = sum(n_pixels_country),
    n_countries_membership = n(),
    .groups = "drop"
  )
d_country_pa_index <- match(d_pa$wdpaid, d_country_pa_check$wdpaid)
if (anyNA(d_country_pa_index) ||
    nrow(d_country_pa_check) != nrow(d_pa) ||
    any(d_pa$n_pixels !=
        d_country_pa_check$n_pixels_country[d_country_pa_index]) ||
    any(d_pa$n_countries !=
        d_country_pa_check$n_countries_membership[d_country_pa_index])) {
  stop("Country aggregation does not reconstruct the PA-level sample")
}

d_country <- d_country_input %>%
  group_by(country_rast) %>%
  summarise(
    delta_mean = weighted.mean(
      delta,
      n_pixels_country,
      na.rm = TRUE
    ),
    tc_protected_mean = weighted.mean(
      tc_protected,
      n_pixels_country,
      na.rm = TRUE
    ),
    frontier_mean = weighted.mean(
      frontier,
      n_pixels_country,
      na.rm = TRUE
    ),
    frontier_level_mean = weighted.mean(
      frontier_level,
      n_pixels_country,
      na.rm = TRUE
    ),
    gap_mean = weighted.mean(u_hat, n_pixels_country, na.rm = TRUE),
    te_mean = weighted.mean(te_jlms, n_pixels_country, na.rm = TRUE),
    tc_control_mean = weighted.mean(
      tc_control,
      n_pixels_country,
      na.rm = TRUE
    ),
    n_pas = n_distinct(wdpaid),
    n_pixels = sum(n_pixels_country),
    total_pa_area = sum(country_pa_area, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  data.frame()

if (sum(d_country$n_pixels) != sum(d_pa$n_pixels) ||
    abs(sum(d_country$total_pa_area) - sum(d_pa$size)) > 1e-6) {
  stop("Country aggregation does not preserve global pixel and PA-area totals")
}

# merge country names (derived from matched data)
country_lookup <- unique(d.mbase[, c("country_rast", "country")])
names(country_lookup) <- c("country_id", "country_name")
d_country <- merge(d_country, country_lookup,
                   by.x = "country_rast", by.y = "country_id", all.x = TRUE)

cat(sprintf("  Countries: %d\n", nrow(d_country)))
cat(sprintf("  Mean country-level efficiency: %.3f\n", mean(d_country$te_mean)))
}


##############################
### 5. LaTeX table ###########
##############################
{
cat("\n=== Block 5: LaTeX table ===\n")

# models for table (matching-covariates fit is kept in model_list but not
# shown in the manuscript table — only the four-column build-up + main spec)
tab_models <- list(sf_ctrl, sf_quad, sf_size, sf_uhet)
tab_names  <- c("Control", "+ Quadratic", "+ Size", "Main: uhet")
n_cols <- length(tab_models)

# covariate labels for display
cov_display <- c(
  "(Intercept)"                = "Intercept",
  "elevation"                  = "Elevation",
  "pop_2000"                   = "Population density",
  "annual_total_precipitation" = "Precipitation",
  "slope"                      = "Slope",
  "temperature"                = "Temperature",
  "access"                     = "Travel time",
  "cropsuit2"                  = "Crop suitability",
  "forest2000_parea"           = "Forest cover (2000)",
  "lu.pressure.est"            = "LU pressure (est.)",
  "lu.vulnerable"              = "Land vulnerability",
  "log_size"                   = "log(PA size + 1)",
  "tc_control"                 = "Control threat level",
  "tc_control_sq"              = "Control threat level$^2$",
  "Zu_(Intercept)"             = "$\\sigma_u$ eq.: intercept",
  "Zu_log_size"                = "$\\sigma_u$ eq.: log(PA size + 1)"
)

# frontier coefficient names (union across the four shown models, in display order)
all_frontier_names <- c("tc_control", "tc_control_sq", "log_size",
                        "Zu_(Intercept)", "Zu_log_size")

tex <- character()

# header
col_headers <- paste0(" & {(", seq_len(n_cols), ")}", collapse = "")
tex <- c(tex, paste0("            ", col_headers, " \\\\"))
name_headers <- paste0(" & {", tab_names, "}", collapse = "")
tex <- c(tex, paste0("            ", name_headers, " \\\\"))
tex <- c(tex, "\\midrule")
tex <- c(tex, "\\addlinespace")
tex <- c(tex, paste0("\\textit{Frontier coefficients}", paste(rep(" &", n_cols), collapse = ""), " \\\\"))

# frontier coefficients with SEs
for (v in all_frontier_names) {
  label <- ifelse(v %in% names(cov_display), cov_display[v], v)

  cells_coef <- character(n_cols)
  cells_se   <- character(n_cols)

  for (j in seq_len(n_cols)) {
    if (!sfa_valid(tab_models[[j]])) {
      # model failed (e.g. hessian inversion) — show dash if variable present
      ct <- tryCatch(summary(tab_models[[j]])$mlRes, error = function(e) NULL)
      if (!is.null(ct) && v %in% rownames(ct)) {
        cells_coef[j] <- "{---}"
        cells_se[j]   <- ""
      } else {
        cells_coef[j] <- ""
        cells_se[j]   <- ""
      }
      next
    }
    ct <- tryCatch(summary(tab_models[[j]])$mlRes, error = function(e) NULL)
    if (!is.null(ct) && v %in% rownames(ct)) {
      coef_val <- ct[v, "Coefficient"]
      se_val   <- ct[v, "Std. Error"]
      p_val    <- ct[v, "Pr(>|z|)"]
      cells_coef[j] <- paste0(sprintf("%.4f", coef_val), star_fn(p_val))
      cells_se[j]   <- paste0("(", sprintf("%.4f", se_val), ")")
    } else {
      cells_coef[j] <- ""
      cells_se[j]   <- ""
    }
  }

  tex <- c(tex, paste0("{", label, "}", paste0(" & ", cells_coef, collapse = ""), " \\\\"))
  tex <- c(tex, paste0("          ", paste0(" & ", cells_se, collapse = ""), " \\\\"))
}

# variance parameters and diagnostics
tex <- c(tex, "\\addlinespace")
tex <- c(tex, "\\midrule")
tex <- c(tex, "\\addlinespace")
tex <- c(tex, paste0("\\textit{Variance parameters}", paste(rep(" &", n_cols), collapse = ""), " \\\\"))

# sigma_u
su_cells <- sapply(tab_models, function(m) {
  if (!sfa_valid(m)) return("{---}")
  p <- tryCatch(coef(m, extraPar = TRUE), error = function(e) NULL)
  if (!is.null(p) && "sigmauSq" %in% names(p)) sprintf("%.4f", sqrt(p["sigmauSq"])) else ""
})
tex <- c(tex, paste0("{$\\sigma_u$}", paste0(" & ", su_cells, collapse = ""), " \\\\"))

# sigma_v
sv_cells <- sapply(tab_models, function(m) {
  if (!sfa_valid(m)) return("{---}")
  p <- tryCatch(coef(m, extraPar = TRUE), error = function(e) NULL)
  if (!is.null(p) && "sigmavSq" %in% names(p)) sprintf("%.4f", sqrt(p["sigmavSq"])) else ""
})
tex <- c(tex, paste0("{$\\sigma_v$}", paste0(" & ", sv_cells, collapse = ""), " \\\\"))

# gamma
gamma_cells <- sapply(tab_models, function(m) {
  if (!sfa_valid(m)) return("{---}")
  p <- tryCatch(coef(m, extraPar = TRUE), error = function(e) NULL)
  if (!is.null(p) && "gamma" %in% names(p)) sprintf("%.3f", p["gamma"]) else ""
})
tex <- c(tex, paste0("{$\\gamma$}", paste0(" & ", gamma_cells, collapse = ""), " \\\\"))

# diagnostics
tex <- c(tex, "\\addlinespace")

# log-likelihood
ll_cells <- sapply(tab_models, function(m) {
  if (!sfa_valid(m)) return("{---}")
  sprintf("%.1f", as.numeric(logLik(m)))
})
tex <- c(tex, paste0("{Log-likelihood}", paste0(" & ", ll_cells, collapse = ""), " \\\\"))

# AIC
aic_cells <- sapply(tab_models, function(m) {
  if (!sfa_valid(m)) return("{---}")
  aic_val <- tryCatch(AIC(m), error = function(e) NA_real_)
  if (!is.finite(aic_val)) "{---}" else sprintf("%.1f", aic_val)
})
tex <- c(tex, paste0("{AIC}", paste0(" & ", aic_cells, collapse = ""), " \\\\"))

# N
n_cells <- sapply(tab_models, function(m) {
  nv <- safe_nobs(m)
  if (is.na(nv)) "{---}" else paste0("{", formatC(nv, format = "d", big.mark = ","), "}")
})
tex <- c(tex, paste0("{N (PAs)}", paste0(" & ", n_cells, collapse = ""), " \\\\"))

# mean efficiency
te_cells <- character(n_cols)
for (j in seq_len(n_cols)) {
  if (!sfa_valid(tab_models[[j]])) { te_cells[j] <- "{---}"; next }
  eff_j <- safe_efficiencies(tab_models[[j]])
  if (!is.null(eff_j)) {
    te_cells[j] <- sprintf("%.3f", mean(exp(-eff_j$u)))
  } else {
    te_cells[j] <- ""
  }
}
tex <- c(tex, paste0("{Mean efficiency}", paste0(" & ", te_cells, collapse = ""), " \\\\"))

# median gap
gap_cells <- character(n_cols)
for (j in seq_len(n_cols)) {
  if (!sfa_valid(tab_models[[j]])) { gap_cells[j] <- "{---}"; next }
  eff_j <- safe_efficiencies(tab_models[[j]])
  if (!is.null(eff_j)) {
    gap_cells[j] <- sprintf("%.4f", median(eff_j$u))
  } else {
    gap_cells[j] <- ""
  }
}
tex <- c(tex, paste0("{Median gap}", paste0(" & ", gap_cells, collapse = "")))

# strip trailing \\ from last line
tex[length(tex)] <- sub(" *\\\\\\\\$", "", tex[length(tex)])

writeLines(tex, "pub/tables/tab.sfa.frontier.tex")
cat("  Saved: pub/tables/tab.sfa.frontier.tex\n")

# --- HTML table via htmlreg for pathreat.est.newest.html ---
cat("  Generating HTML table (htmlreg) for SFA results\n")

# build texreg objects from sfacross models (no native texreg support)
sfa_to_texreg <- function(m, label = "") {
  if (is.null(m)) return(NULL)
  valid <- sfa_valid(m)
  s <- tryCatch(summary(m)$mlRes, error = function(e) NULL)
  if (is.null(s)) return(NULL)
  coef_names <- rownames(s)
  display_names <- vapply(coef_names, function(v)
    if (v %in% names(cov_display)) cov_display[v] else v, character(1))
  p <- tryCatch(coef(m, extraPar = TRUE), error = function(e) NULL)
  eff_j <- safe_efficiencies(m)
  gof_vals <- c(
    if (valid) as.numeric(logLik(m)) else NA,
    if (valid) tryCatch(unname(AIC(m)), error = function(e) NA_real_) else NA,
    safe_nobs(m),
    if (!is.null(p) && valid) unname(sqrt(p["sigmauSq"])) else NA,
    if (!is.null(p) && valid) unname(sqrt(p["sigmavSq"])) else NA,
    if (!is.null(p) && valid) unname(p["gamma"]) else NA,
    if (!is.null(eff_j)) mean(exp(-eff_j$u)) else NA,
    if (!is.null(eff_j)) median(eff_j$u) else NA
  )
  # for invalid models, set coefficients to NA so htmlreg shows blanks
  coefs  <- if (valid) unname(s[, "Coefficient"]) else rep(NA_real_, nrow(s))
  ses    <- if (valid) unname(s[, "Std. Error"])   else rep(NA_real_, nrow(s))
  pvals  <- if (valid) unname(s[, "Pr(>|z|)"])    else rep(NA_real_, nrow(s))
  createTexreg(
    coef.names  = unname(display_names),
    coef        = coefs,
    se          = ses,
    pvalues     = pvals,
    gof.names   = c("Log-likelihood", "AIC", "N (PAs)",
                     "sigma_u", "sigma_v", "gamma",
                     "Mean efficiency", "Median gap"),
    gof         = gof_vals,
    gof.decimal = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)
  )
}

tr_list <- lapply(tab_models, sfa_to_texreg)
tr_list <- Filter(Negate(is.null), tr_list)

h_sfa <- htmlreg(
  tr_list,
  custom.model.names = tab_names,
  stars = c(0.01, 0.05, 0.1),
  digits = 4,
  table = FALSE
)

cat(file = paste0(paths$results_dir, "pathreat.est.newest.html"),
    c("<br><br><b>Stochastic Frontier Analysis: PA-level efficiency gaps (threat_composite)</b>", h_sfa),
    append = TRUE)
cat("  Appended to: results/pathreat.est.newest.html\n")
}


##############################
### 6. Save results ##########
##############################
{
cat("\n=== Block 6: Saving results ===\n")

# PA-level results
saveRDS(d_pa, paste0(paths$results_dir, "pathreat.sfa.est.Rds"))
cat("  Saved: results/pathreat.sfa.est.Rds\n")

# country-level results
saveRDS(d_country, paste0(paths$results_dir, "pathreat.sfa.country.est.Rds"))
cat("  Saved: results/pathreat.sfa.country.est.Rds\n")

# model objects
saveRDS(model_list, paste0(paths$results_dir, "pathreat.sfa.models.est.Rds"))
cat("  Saved: results/pathreat.sfa.models.est.Rds\n")

# HTML summary
h_main <- summary(sf_main)
summary_lines <- trimws(capture.output(print(h_main)), which = "right")
writeLines(summary_lines, paste0(paths$results_dir, "pathreat.sfa.summary.est.txt"))
cat("  Saved: results/pathreat.sfa.summary.est.txt\n")

rm(d.mbase); gc()
}

##############################
### 7. Figures ###############
##############################
{
cat("\n=== Block 7: Figures ===\n")

# IUCN class labels
iucn_labels <- c(
  "Strict protection" = "Strict (Ia\u2013IV)",
  "Less strict"       = "Multi-use (V\u2013VI)",
  "Not Reported"      = "Not Reported"
)
d_pa$iucn_label <- factor(iucn_labels[d_pa$iucn_class],
                          levels = iucn_labels)

iucn_colors <- c(
  "Strict (Ia\u2013IV)"    = "dodgerblue3",
  "Multi-use (V\u2013VI)"  = "darkorange2",
  "Not Reported"            = "grey55"
)


### Global means (pixel-weighted, for reference lines) ###
global_delta    <- weighted.mean(d_pa$delta, d_pa$n_pixels)
global_frontier <- weighted.mean(d_pa$frontier, d_pa$n_pixels)

### Figure A: Actual vs. Frontier scatter ###
cat("  Figure A: Actual vs. frontier scatter\n")

ax_lo <- -0.22
ax_hi <- 0.12

p_scatter <- ggplot(d_pa, aes(x = frontier, y = delta)) +
  # shade regions
  annotate("polygon",
           x = c(ax_lo, ax_hi, ax_hi, ax_lo),
           y = c(ax_lo, ax_hi, ax_hi + 0.5, ax_lo + 0.5),
           fill = "coral", alpha = 0.04) +
  annotate("polygon",
           x = c(ax_lo, ax_hi, ax_hi, ax_lo),
           y = c(ax_lo, ax_hi, ax_hi - 0.5, ax_lo - 0.5),
           fill = "steelblue", alpha = 0.04) +
  # global reference lines
  geom_vline(xintercept = global_frontier, linetype = "dashed",
             color = "firebrick", linewidth = 0.4, alpha = 0.7) +
  geom_hline(yintercept = global_delta, linetype = "dashed",
             color = "steelblue", linewidth = 0.4, alpha = 0.7) +
  annotate("text", x = global_frontier, y = ax_lo + 0.008,
           label = sprintf("Global frontier: %.3f", global_frontier),
           hjust = -0.05, size = 2.5, color = "firebrick4") +
  annotate("text", x = ax_lo + 0.008, y = global_delta,
           label = sprintf("Global actual: %.3f", global_delta),
           vjust = -0.5, hjust = 0, size = 2.5, color = "steelblue4") +
  # 45-degree line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "grey30", linewidth = 0.5) +
  geom_point(aes(fill = iucn_label, size = log(n_pixels)),
             shape = 21, alpha = 0.35, stroke = 0.1, color = "grey40") +
  scale_fill_manual(values = iucn_colors, name = NULL) +
  scale_size_continuous(range = c(0.4, 3.5), guide = "none") +
  # region labels
  annotate("text", x = ax_lo + 0.01, y = ax_hi - 0.005,
           label = "Underperformance\n(actual > frontier)",
           hjust = 0, vjust = 1, size = 2.8, fontface = "italic",
           color = "coral4") +
  annotate("text", x = ax_hi - 0.01, y = ax_lo + 0.005,
           label = "Overperformance\n(actual < frontier)",
           hjust = 1, vjust = 0, size = 2.8, fontface = "italic",
           color = "steelblue4") +
  # stats
  annotate("text", x = ax_hi - 0.005, y = ax_hi - 0.005,
           label = sprintf("Mean eff. = %.3f\n\u03B3 = %.3f\nN = %s",
                           mean(d_pa$te_jlms), gamma_val,
                           format(nrow(d_pa), big.mark = ",")),
           hjust = 1, vjust = 1, size = 2.8, color = "grey20") +
  labs(x = expression("Frontier " * hat(delta)[frontier] *
                      " (maximum achievable threat reduction)"),
       y = expression("Observed " * delta *
                      " (actual treatment effect)")) +
  coord_fixed(ratio = 1, xlim = c(ax_lo, ax_hi), ylim = c(ax_lo, ax_hi)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3, alpha = 0.7)))

ggsave(paste0(paths$figures_dir, "fig.sfa.actual_vs_frontier.jpg"),
       plot = p_scatter, width = 16, height = 16, units = "cm", dpi = 300)
cat("    Saved: pub/figures/fig.sfa.actual_vs_frontier.jpg\n")


### Figure A2: Actual vs. Frontier scatter — country-level ###
cat("  Figure A2: Actual vs. frontier scatter (country-level)\n")

ax_lo_c <- min(c(d_country$delta_mean, d_country$frontier_mean), na.rm = TRUE) - 0.01
ax_hi_c <- max(c(d_country$delta_mean, d_country$frontier_mean), na.rm = TRUE) + 0.01

p_scatter_country <- ggplot(d_country, aes(x = frontier_mean, y = delta_mean)) +
  # zero line for actual effect
  geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +
  # shade regions
  annotate("polygon",
           x = c(ax_lo_c, ax_hi_c, ax_hi_c, ax_lo_c),
           y = c(ax_lo_c, ax_hi_c, ax_hi_c + 0.5, ax_lo_c + 0.5),
           fill = "coral", alpha = 0.04) +
  annotate("polygon",
           x = c(ax_lo_c, ax_hi_c, ax_hi_c, ax_lo_c),
           y = c(ax_lo_c, ax_hi_c, ax_hi_c - 0.5, ax_lo_c - 0.5),
           fill = "steelblue", alpha = 0.04) +
  # global reference lines
  geom_vline(xintercept = global_frontier, linetype = "dashed",
             color = "firebrick", linewidth = 0.4, alpha = 0.7) +
  geom_hline(yintercept = global_delta, linetype = "dashed",
             color = "steelblue", linewidth = 0.4, alpha = 0.7) +
  annotate("text", x = global_frontier, y = ax_lo_c + 0.005,
           label = sprintf("Global frontier: %.3f", global_frontier),
           hjust = -0.05, size = 2.5, color = "firebrick4") +
  annotate("text", x = ax_lo_c + 0.005, y = global_delta,
           label = sprintf("Global actual: %.3f", global_delta),
           vjust = -0.5, hjust = 0, size = 2.5, color = "steelblue4") +
  # 45-degree line
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "grey30", linewidth = 0.5) +
  geom_point(aes(size = n_pas), shape = 21,
             fill = "steelblue", alpha = 0.5, stroke = 0.3, color = "grey30") +
  scale_size_continuous(range = c(1.5, 8), name = "N (PAs)") +
  # country labels
  ggrepel::geom_text_repel(
    data = d_country[d_country$n_pas >= quantile(d_country$n_pas, 0.65), ],
    aes(label = country_name), size = 2.2, color = "grey30",
    max.overlaps = 20, segment.color = "grey70", segment.size = 0.2,
    min.segment.length = 0.1, box.padding = 0.25, seed = 42) +
  # region labels
  annotate("text", x = ax_lo_c + 0.005, y = ax_hi_c - 0.003,
           label = "Underperformance",
           hjust = 0, vjust = 1, size = 2.8, fontface = "italic",
           color = "coral4") +
  annotate("text", x = ax_hi_c - 0.005, y = ax_lo_c + 0.003,
           label = "Overperformance",
           hjust = 1, vjust = 0, size = 2.8, fontface = "italic",
           color = "steelblue4") +
  # stats
  annotate("text", x = ax_hi_c - 0.003, y = ax_hi_c - 0.003,
           label = sprintf("N = %d countries", nrow(d_country)),
           hjust = 1, vjust = 1, size = 2.8, color = "grey20") +
  labs(x = expression("Country mean frontier " * hat(delta)[frontier]),
       y = expression("Country mean observed " * delta)) +
  coord_fixed(ratio = 1, xlim = c(ax_lo_c, ax_hi_c), ylim = c(ax_lo_c, ax_hi_c)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  )

ggsave(paste0(paths$figures_dir, "fig.sfa.actual_vs_frontier.country.jpg"),
       plot = p_scatter_country, width = 16, height = 16, units = "cm", dpi = 300)
cat("    Saved: pub/figures/fig.sfa.actual_vs_frontier.country.jpg\n")


### Figure A3: Frontier curve vs. actual effects ###
cat("  Figure A3: Frontier curve (tc_control vs. delta)\n")

# predicted frontier line over range of tc_control (holding log_size at median)
# Note: sf_main fits the LEVEL frontier (tc_protected) without log_size; we
# convert to the delta-equivalent frontier (= level frontier - tc_control)
# so the y-axis interpretation (treatment effect delta) is preserved.
tc_seq <- seq(min(d_pa$tc_control), max(d_pa$tc_control), length.out = 200)
frontier_coefs <- coef(sf_main)
frontier_intercept <- ifelse("(Intercept)" %in% names(frontier_coefs),
                             frontier_coefs["(Intercept)"], 0)
get_coef <- function(nm) if (nm %in% names(frontier_coefs)) frontier_coefs[nm] else 0
size_med  <- median(d_pa$log_size)
fr_level  <- frontier_intercept +
  get_coef("tc_control")    * tc_seq +
  get_coef("tc_control_sq") * tc_seq^2 +
  get_coef("log_size")      * size_med
frontier_line   <- fr_level - tc_seq                               # delta-equivalent
d_frontier_line <- data.frame(tc_control = tc_seq, frontier = frontier_line)
frontier_label  <- sprintf("Frontier (median PA: %.0f km\u00B2)", exp(size_med) - 1)

# raw PA-level scatter
p_frontier_raw <- ggplot(d_pa, aes(x = tc_control, y = delta)) +
  geom_point(aes(fill = iucn_label, size = log(n_pixels)),
             shape = 21, alpha = 0.3, stroke = 0.1, color = "grey40") +
  geom_line(data = d_frontier_line, aes(x = tc_control, y = frontier),
            color = "firebrick", linewidth = 1) +
  scale_fill_manual(values = iucn_colors, name = NULL) +
  scale_size_continuous(range = c(0.4, 3.5), guide = "none") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_hline(yintercept = global_delta, linetype = "dashed",
             color = "steelblue", linewidth = 0.4, alpha = 0.7) +
  geom_hline(yintercept = global_frontier, linetype = "dashed",
             color = "firebrick", linewidth = 0.4, alpha = 0.7) +
  annotate("text", x = min(tc_seq), y = global_delta,
           label = sprintf("Global actual: %.3f", global_delta),
           vjust = -0.5, hjust = 0, size = 2.5, color = "steelblue4") +
  annotate("text", x = min(tc_seq), y = global_frontier,
           label = sprintf("Global frontier: %.3f", global_frontier),
           vjust = 1.5, hjust = 0, size = 2.5, color = "firebrick4") +
  annotate("text", x = max(tc_seq) * 0.95, y = min(frontier_line) * 0.9,
           label = frontier_label, color = "firebrick", fontface = "bold",
           size = 3, hjust = 1) +
  scale_x_continuous(limits = range(d_pa$tc_control), expand = expansion(mult = 0.02)) +
  labs(x = "Control threat level (counterfactual threat composite)",
       y = expression("Treatment effect " * delta * " (threat reduction)")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3, alpha = 0.7)))

ggsave(paste0(paths$figures_dir, "fig.sfa.frontier_curve.jpg"),
       plot = p_frontier_raw, width = 18, height = 14, units = "cm", dpi = 300)
cat("    Saved: pub/figures/fig.sfa.frontier_curve.jpg\n")

# binned version
cat("  Figure A3b: Frontier curve (quantile-binned)\n")

# quantile bins by IUCN class
n_bins <- 200
tc_breaks <- unique(quantile(d_pa$tc_control, probs = seq(0, 1, length.out = n_bins + 1)))
d_pa$tc_bin <- cut(d_pa$tc_control,
                   breaks = tc_breaks,
                   include.lowest = TRUE, labels = FALSE)

d_binned <- d_pa %>%
  group_by(tc_bin, iucn_label) %>%
  summarise(
    tc_mid     = mean(tc_control),
    delta_mean = mean(delta),
    n          = n(),
    .groups = "drop"
  )

p_frontier_curve <- ggplot() +
  geom_point(data = d_binned, aes(x = tc_mid, y = delta_mean,
                                   color = iucn_label),
             size = 2, alpha = 0.7) +
  geom_line(data = d_frontier_line, aes(x = tc_control, y = frontier),
            color = "firebrick", linewidth = 1) +
  scale_color_manual(values = iucn_colors, name = NULL) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_hline(yintercept = global_delta, linetype = "dashed",
             color = "steelblue", linewidth = 0.4, alpha = 0.7) +
  geom_hline(yintercept = global_frontier, linetype = "dashed",
             color = "firebrick", linewidth = 0.4, alpha = 0.7) +
  annotate("text", x = min(tc_seq), y = global_delta,
           label = sprintf("Global actual: %.3f", global_delta),
           vjust = -0.5, hjust = 0, size = 2.5, color = "steelblue4") +
  annotate("text", x = min(tc_seq), y = global_frontier,
           label = sprintf("Global frontier: %.3f", global_frontier),
           vjust = 1.5, hjust = 0, size = 2.5, color = "firebrick4") +
  annotate("text", x = max(tc_seq) * 0.95, y = min(frontier_line) * 0.9,
           label = frontier_label, color = "firebrick", fontface = "bold",
           size = 3, hjust = 1) +
  scale_x_continuous(limits = range(d_pa$tc_control), expand = expansion(mult = 0.02)) +
  labs(x = "Control threat level (counterfactual threat composite)",
       y = expression("Treatment effect " * delta * " (threat reduction)")) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.9)))

ggsave(paste0(paths$figures_dir, "fig.sfa.frontier_curve.binned.jpg"),
       plot = p_frontier_curve, width = 18, height = 14, units = "cm", dpi = 300)
cat("    Saved: pub/figures/fig.sfa.frontier_curve.binned.jpg\n")


### Figure A4: Normalized frontier curve with binned scatter ###
cat("  Figure A4: Frontier curve (normalized, % of control mean)\n")

# global control mean for normalization
g_global <- readRDS(paths$est_global)
g_tc <- g_global[g_global$variable == "threat_composite", ]
ctrl_mean <- g_tc$control_mean
global_effect_norm <- g_tc$coef / ctrl_mean * 100

# frontier line from model coefficients (evaluated at median log_size)
# tc_seq and frontier_coefs already computed above for Figure A3
d_frontier_norm <- data.frame(
  tc_control    = tc_seq,
  frontier_norm = frontier_line / ctrl_mean * 100
)

# binned scatter (150 equal-count bins)
d_pa$tc_bin_norm <- ntile(d_pa$tc_control, 150)
d_binned_norm <- d_pa %>%
  group_by(tc_bin_norm) %>%
  summarise(
    tc_mid     = mean(tc_control),
    delta_norm = mean(delta) / ctrl_mean * 100,
    size       = mean(size),
    .groups = "drop"
  )

p_frontier_norm <- ggplot(d_binned_norm, aes(x = tc_mid, y = delta_norm)) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 11, face = "bold")
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = global_effect_norm, linetype = "dashed",
             color = "steelblue", linewidth = 0.5) +
  # frontier curve (from SFA model, at median log_size)
  geom_line(data = d_frontier_norm, aes(x = tc_control, y = frontier_norm),
            color = "#B2182B", linewidth = 1.2) +
  # binned scatter
  geom_point(aes(size = size), alpha = 0.5, shape = 16, color = "grey30") +
  scale_size_continuous(name = expression("PA area (km"^2*")"), range = c(0.5, 4),
                        trans = "log10", breaks = c(10, 100, 1000, 10000)) +
  # loess smoother on binned data
  geom_smooth(method = "loess", span = 1.2, linewidth = 1,
              se = TRUE, alpha = 0.15, color = "#2166AC", fill = "#2166AC") +
  # labels
  annotate("text", x = -Inf, y = global_effect_norm, label = "Global average",
           hjust = -0.05, vjust = -0.5, size = 3, color = "steelblue4",
           fontface = "italic") +
  annotate("text", x = max(d_binned_norm$tc_mid) * 0.95,
           y = min(d_frontier_norm$frontier_norm[
             d_frontier_norm$tc_control <= max(d_binned_norm$tc_mid)]) * 0.9,
           label = frontier_label, color = "#B2182B", fontface = "bold",
           size = 3, hjust = 1) +
  annotate("text", x = max(d_binned_norm$tc_mid) * 0.55,
           y = min(d_binned_norm$delta_norm) * 0.6,
           label = sprintf("\u03B3 = %.3f\nMean eff. = %.3f",
                           gamma_val, mean(d_pa$te_jlms)),
           size = 3, color = "grey30", hjust = 0) +
  labs(
    x = "Control threat level (counterfactual threat composite)",
    y = "Effect on overall threat index (% of control mean)"
  )

ggsave(paste0(paths$figures_dir, "fig.sfa.frontier_curve.norm.jpg"),
       plot = p_frontier_norm, width = 18, height = 14, units = "cm", dpi = 300)
cat("    Saved: pub/figures/fig.sfa.frontier_curve.norm.jpg\n")

# unbinned version: 10% random sample of PAs
cat("  Figure A4b: Frontier curve (normalized, unbinned, 10% sample)\n")
set.seed(42)
d_pa_sample <- d_pa[sample(nrow(d_pa), size = round(nrow(d_pa) * 0.10)), ]
d_pa_sample$delta_norm <- d_pa_sample$delta / ctrl_mean * 100

p_frontier_norm_unbinned <- ggplot(d_pa_sample, aes(x = tc_control, y = delta_norm)) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 11, face = "bold")
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = global_effect_norm, linetype = "dashed",
             color = "steelblue", linewidth = 0.5) +
  geom_line(data = d_frontier_norm, aes(x = tc_control, y = frontier_norm),
            color = "#B2182B", linewidth = 1.2) +
  geom_point(aes(size = size), alpha = 0.25, shape = 16, color = "grey30") +
  scale_size_continuous(name = expression("PA area (km"^2*")"), range = c(0.3, 3),
                        trans = "log10", breaks = c(10, 100, 1000, 10000)) +
  annotate("text", x = -Inf, y = global_effect_norm, label = "Global average",
           hjust = -0.05, vjust = -0.5, size = 3, color = "steelblue4",
           fontface = "italic") +
  annotate("text", x = max(d_pa_sample$tc_control) * 0.95,
           y = min(d_frontier_norm$frontier_norm[
             d_frontier_norm$tc_control <= max(d_pa_sample$tc_control)]) * 0.9,
           label = frontier_label, color = "#B2182B", fontface = "bold",
           size = 3, hjust = 1) +
  labs(
    x = "Control threat level (counterfactual threat composite)",
    y = "Effect on overall threat index (% of control mean)"
  )

ggsave(paste0(paths$figures_dir, "fig.sfa.frontier_curve.norm.unbinned.jpg"),
       plot = p_frontier_norm_unbinned, width = 18, height = 14, units = "cm", dpi = 300)
cat("    Saved: pub/figures/fig.sfa.frontier_curve.norm.unbinned.jpg\n")


### Figure B: Efficiency distribution ###
cat("  Figure B: Efficiency distribution\n")

te_quants <- quantile(d_pa$te_jlms, probs = c(0.25, 0.5, 0.75))

p_dist <- ggplot(d_pa, aes(x = te_jlms)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 fill = "steelblue", color = "white", alpha = 0.7) +
  geom_density(color = "darkblue", linewidth = 0.8) +
  geom_vline(xintercept = te_quants, linetype = "dashed",
             color = "grey30", linewidth = 0.4) +
  annotate("text", x = te_quants, y = Inf,
           label = c(sprintf("P25: %.3f", te_quants[1]),
                     sprintf("Median: %.3f", te_quants[2]),
                     sprintf("P75: %.3f", te_quants[3])),
           vjust = 1.5, hjust = c(1.1, 1.1, -0.1), size = 2.8, color = "grey30") +
  labs(x = expression("PA efficiency score " * (e^{-hat(u)})),
       y = "Density") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave(paste0(paths$figures_dir, "fig.sfa.efficiency_distribution.jpg"),
       plot = p_dist, width = 16, height = 10, units = "cm", dpi = 300)
cat("    Saved: pub/figures/fig.sfa.efficiency_distribution.jpg\n")


### Figure C: Country lollipop (all countries) ###
cat("  Figure D: Country lollipop\n")

d_lollipop <- d_country %>%
  mutate(
    country_label = ifelse(is.na(country_name),
                           as.character(country_rast), country_name),
    delta_norm    = delta_mean / ctrl_mean * 100,
    frontier_norm = frontier_mean / ctrl_mean * 100
  ) %>%
  arrange(desc(delta_norm - frontier_norm))

d_lollipop$country_label <- factor(d_lollipop$country_label,
                                    levels = d_lollipop$country_label)

global_delta_norm    <- global_delta / ctrl_mean * 100
global_frontier_norm <- global_frontier / ctrl_mean * 100

p_lollipop <- ggplot(d_lollipop, aes(y = country_label)) +
  # global estimate lines
  geom_vline(xintercept = global_delta_norm, linewidth = 0.5, linetype = "dashed",
             color = "steelblue") +
  geom_vline(xintercept = global_frontier_norm, linewidth = 0.5, linetype = "dashed",
             color = "firebrick") +
  annotate("text", x = global_delta_norm, y = Inf,
           label = sprintf("Global actual: %.1f%%", global_delta_norm),
           vjust = -0.3, hjust = 0.5, size = 2.5, color = "steelblue4") +
  annotate("text", x = global_frontier_norm, y = Inf,
           label = sprintf("Global frontier: %.1f%%", global_frontier_norm),
           vjust = -1.5, hjust = 0.5, size = 2.5, color = "firebrick4") +
  # segments and points
  geom_segment(aes(x = delta_norm, xend = frontier_norm,
                   yend = country_label),
               color = "grey60", linewidth = 0.4) +
  geom_point(aes(x = delta_norm, color = "Actual effect", shape = "Actual effect"),
             size = 1.8) +
  geom_point(aes(x = frontier_norm, color = "Frontier", shape = "Frontier"),
             size = 1.8) +
  scale_color_manual(values = c("Actual effect" = "steelblue",
                                "Frontier" = "firebrick"),
                     name = NULL) +
  scale_shape_manual(values = c("Actual effect" = 16,
                                "Frontier" = 1),
                     name = NULL) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  labs(x = "Treatment effect (% of control mean)",
       y = NULL) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    axis.title.x = element_text(size = 10),
    axis.text.y = element_text(size = 5),
    axis.text.x = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.2),
    plot.margin = margin(20, 5, 5, 5, "pt")
  )

ggsave(paste0(paths$figures_dir, "fig.sfa.country_lollipop.jpg"),
       plot = p_lollipop, width = 18, height = 30, units = "cm", dpi = 300)
cat("    Saved: pub/figures/fig.sfa.country_lollipop.jpg\n")

### Figure E: Country scatter — effect vs. potential (normalized) ###
cat("  Figure E: Country scatter (effect vs. potential, normalized)\n")

d_scatter_norm <- d_country %>%
  mutate(
    delta_norm    = delta_mean / ctrl_mean * 100,
    frontier_norm = frontier_mean / ctrl_mean * 100
  )

p_country_scatter <- ggplot(d_scatter_norm, aes(x = frontier_norm, y = delta_norm)) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "black") +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed",
              color = "grey30", linewidth = 0.5) +
  geom_vline(xintercept = global_frontier_norm, linetype = "dashed",
             color = "firebrick", linewidth = 0.4, alpha = 0.7) +
  geom_hline(yintercept = global_delta_norm, linetype = "dashed",
             color = "steelblue", linewidth = 0.4, alpha = 0.7) +
  annotate("text", x = global_frontier_norm, y = min(d_scatter_norm$delta_norm),
           label = sprintf("Global potential: %.1f%%", global_frontier_norm),
           hjust = -0.05, vjust = 0, size = 2.5, color = "firebrick4") +
  annotate("text", x = min(d_scatter_norm$frontier_norm),
           y = global_delta_norm,
           label = sprintf("Global actual: %.1f%%", global_delta_norm),
           hjust = 0, vjust = -0.5, size = 2.5, color = "steelblue4") +
  geom_point(aes(size = total_pa_area), shape = 21,
             fill = "steelblue", alpha = 0.5, stroke = 0.3, color = "grey30") +
  scale_size_continuous(range = c(1.5, 8), name = expression("Total PA area (km"^2*")"),
                        trans = scales::trans_new("log2", log2, function(x) 2^x),
                        breaks = c(100, 1000, 10000, 100000),
                        labels = scales::comma) +
  ggrepel::geom_text_repel(
    data = d_scatter_norm[d_scatter_norm$total_pa_area >= quantile(d_scatter_norm$total_pa_area, 0.65), ],
    aes(label = country_name), size = 2.2, color = "grey30",
    max.overlaps = 20, segment.color = "grey70", segment.size = 0.2,
    min.segment.length = 0.1, box.padding = 0.25, seed = 42) +
  annotate("text", x = min(d_scatter_norm$frontier_norm) + 1,
           y = max(d_scatter_norm$delta_norm) - 1,
           label = "Underperformance",
           hjust = 0, vjust = 1, size = 2.8, fontface = "italic", color = "coral4") +
  annotate("text", x = max(d_scatter_norm$frontier_norm) - 1,
           y = min(d_scatter_norm$delta_norm) + (max(d_scatter_norm$delta_norm) - min(d_scatter_norm$delta_norm)) * 0.12,
           label = "Overperformance",
           hjust = 1, vjust = 0, size = 2.8, fontface = "italic", color = "steelblue4") +
  labs(x = "Potential effect (% of control mean)",
       y = "Effect on threat index (% of control mean)") +
  coord_cartesian(ylim = range(d_scatter_norm$delta_norm) * c(1.03, 1.03)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 9),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  )

ggsave(paste0(paths$figures_dir, "fig.sfa.country_scatter.norm.jpg"),
       plot = p_country_scatter, width = 16, height = 14, units = "cm", dpi = 300)
cat("    Saved: pub/figures/fig.sfa.country_scatter.norm.jpg\n")


cat("  All figures saved.\n")
}


cat("\n=== SFA estimation complete ===\n")
cat("PA-level results:  results/pathreat.sfa.est.Rds\n")
cat("Country-level:     results/pathreat.sfa.country.est.Rds\n")
cat("Model objects:     results/pathreat.sfa.models.est.Rds\n")
cat("LaTeX table:       pub/tables/tab.sfa.frontier.tex\n")
cat("Figures:           pub/figures/fig.sfa.*.jpg\n")
