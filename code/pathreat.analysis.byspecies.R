##################################################
### pathreat.analysis.byspecies.R ################
### Species-range fixed-effect estimates #########
##################################################
#
# Selects matched pairs by whether the protected pixel lies inside each
# species range. Both observations in every selected pair are retained.
# Species-specific regressions adapt the fixed effects to the variation in
# country and biome within the selected sample.
#
# Inputs:
#   data/store/pathreat.data.merge.matched.fst
#   data/store/pathreat.data.species.Rds
#   data/store/pathreat.data.species-pa-cover.Rds (coverage metadata only)
#   data/store/pathreat.data.rasterbase.tif
#   data/store/species_rasters/sp_{taxon}_{id}.tif
#
# Output:
#   results/pathreat.byspecies.est.Rds
#
# Optional environment controls:
#   PATHREAT_BYSPECIES_CACHE_DIR
#   PATHREAT_BYSPECIES_OUTPUT
#   PATHREAT_BYSPECIES_CORES
#

library(dplyr)
library(fixest)
library(fst)
library(terra)

source("code/pathreat.analysis.config.R")

env_or_default <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

tmp_dir <- env_or_default(
  "PATHREAT_BYSPECIES_CACHE_DIR",
  "data/tmp/byspecies_regression_v3_rows"
)
outpath <- env_or_default(
  "PATHREAT_BYSPECIES_OUTPUT",
  file.path(paths$results_dir, "pathreat.byspecies.est.Rds")
)
cores_value <- env_or_default("PATHREAT_BYSPECIES_CORES", "4")
if (!grepl("^[1-9][0-9]*$", cores_value)) {
  stop("PATHREAT_BYSPECIES_CORES must be a positive integer.")
}
no_cores <- as.integer(cores_value)
if (is.na(no_cores)) {
  stop("PATHREAT_BYSPECIES_CORES is too large to represent as an integer.")
}

estimator_contract_version <- "pathreat-byspecies-range-fe-v3"
species_metadata_path <- "data/store/pathreat.data.species.Rds"
species_coverage_path <- "data/store/pathreat.data.species-pa-cover.Rds"
rasterbase_path <- "data/store/pathreat.data.rasterbase.tif"
species_raster_dir <- "data/store/species_rasters"
cache_manifest_name <- "manifest.rds"

fingerprint_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required provenance input does not exist: ", path)
  }
  checksum <- unname(tools::md5sum(path))
  if (is.na(checksum)) {
    stop("Could not checksum required provenance input: ", path)
  }
  list(path = path, md5 = checksum)
}

save_rds_atomic <- function(object, path) {
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  temporary_path <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = dirname(path)
  )
  on.exit(unlink(temporary_path), add = TRUE)
  saveRDS(object, temporary_path)
  if (!file.rename(temporary_path, path)) {
    stop("Could not atomically move RDS into place: ", path)
  }
  invisible(path)
}

cat(sprintf(
  paste0(
    "=== Species estimator controls ===\n",
    "  Cache: %s\n",
    "  Output: %s\n",
    "  Cores: %d\n"
  ),
  tmp_dir,
  outpath,
  no_cores
))

############################
### metadata + coverage ####
############################
{
cat("=== Loading species metadata and coverage ===\n")

sp_meta <- readRDS(species_metadata_path)
required_meta <- c(
  "species_id",
  "species_name",
  "taxon",
  "category",
  "raster_file"
)
missing_meta <- setdiff(required_meta, names(sp_meta))
if (length(missing_meta) > 0) {
  stop("Species metadata is missing: ", paste(missing_meta, collapse = ", "))
}
if (anyDuplicated(sp_meta$species_id)) {
  stop("Species metadata contains duplicated species_id values.")
}

sp_cover <- readRDS(species_coverage_path)
required_cover <- c("species_id", "wdpaid", "overlap_area", "range_total")
missing_cover <- setdiff(required_cover, names(sp_cover))
if (length(missing_cover) > 0) {
  stop("Species coverage data are missing: ", paste(missing_cover, collapse = ", "))
}

sp_coverage <- sp_cover %>%
  group_by(species_id) %>%
  dplyr::summarize(
    range_total = first(range_total),
    protected_area = sum(overlap_area[wdpaid > 0], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(share_protected = protected_area / range_total) %>%
  select(species_id, range_total, share_protected)

cat(sprintf(
  "  Metadata: %s species; coverage: %s species\n",
  format(nrow(sp_meta), big.mark = ","),
  format(nrow(sp_coverage), big.mark = ",")
))

rm(sp_cover)
gc()
}

############################
### canonical row cache ####
############################
{
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

manifest_path <- file.path(tmp_dir, cache_manifest_name)
cache_entries <- list.files(tmp_dir, all.files = TRUE, no.. = TRUE)
if (!file.exists(manifest_path) && length(cache_entries) > 0L) {
  stop(
    "Refusing to reuse nonempty species cache without a manifest: ",
    tmp_dir
  )
}

raster_paths <- file.path(species_raster_dir, sp_meta$raster_file)
raster_info <- file.info(raster_paths)
raster_mtime_utc <- rep(NA_character_, nrow(raster_info))
has_mtime <- !is.na(raster_info$mtime)
raster_mtime_utc[has_mtime] <- format(
  raster_info$mtime[has_mtime],
  "%Y-%m-%dT%H:%M:%OS6Z",
  tz = "UTC"
)
raster_manifest <- data.frame(
  species_id = sp_meta$species_id,
  raster_file = as.character(sp_meta$raster_file),
  size_bytes = as.numeric(raster_info$size),
  mtime_utc = raster_mtime_utc,
  stringsAsFactors = FALSE,
  row.names = NULL
)

package_names <- c("dplyr", "fixest", "fst", "terra")
cache_manifest_current <- list(
  manifest_schema_version = 3L,
  estimator_source = fingerprint_file("code/pathreat.analysis.byspecies.R"),
  shared_source = fingerprint_file("code/pathreat.analysis.config.R"),
  estimator_contract_version = estimator_contract_version,
  matched_fst = fingerprint_file(paths$data_matched),
  species_metadata = fingerprint_file(species_metadata_path),
  rasterbase = fingerprint_file(rasterbase_path),
  species_rasters = raster_manifest,
  runtime = list(
    r_version = as.character(getRversion()),
    r_platform = R.version$platform,
    packages = setNames(
      vapply(
        package_names,
        function(package_name) as.character(utils::packageVersion(package_name)),
        character(1)
      ),
      package_names
    )
  )
)

if (file.exists(manifest_path)) {
  cache_manifest_stored <- tryCatch(
    readRDS(manifest_path),
    error = function(err) err
  )
  if (inherits(cache_manifest_stored, "error")) {
    stop(
      "Refusing to reuse species cache with an unreadable manifest: ",
      tmp_dir,
      " (",
      cache_manifest_stored$message,
      ")"
    )
  }

  manifest_components <- names(cache_manifest_current)
  incompatible_components <- manifest_components[
    !vapply(
      manifest_components,
      function(component) {
        identical(
          cache_manifest_stored[[component]],
          cache_manifest_current[[component]]
        )
      },
      logical(1)
    )
  ]
  if (length(incompatible_components) > 0L) {
    stop(
      "Refusing to reuse incompatible species cache ",
      tmp_dir,
      ". Manifest mismatch: ",
      paste(incompatible_components, collapse = ", ")
    )
  }
} else {
  save_rds_atomic(cache_manifest_current, manifest_path)
  cat("  Wrote cache manifest: ", manifest_path, "\n", sep = "")
}

done_files <- list.files(tmp_dir, pattern = "^sp_[0-9]+\\.rds$")
done_ids <- as.integer(sub("^sp_([0-9]+)\\.rds$", "\\1", done_files))
todo_idx <- which(!sp_meta$species_id %in% done_ids)

cat(sprintf(
  "  Species estimates: %s total, %s cached, %s remaining\n",
  format(nrow(sp_meta), big.mark = ","),
  format(length(done_ids), big.mark = ","),
  format(length(todo_idx), big.mark = ",")
))
}

########################################
### estimate missing cache rows only ###
########################################
if (length(todo_idx) > 0) {
  cat("\n=== Loading matched pixels for uncached species ===\n")

  d_match <- read_fst(
    paths$data_matched,
    columns = c(
      "matched_pair_id",
      "treat",
      "country_rast",
      "biome_raster",
      "x",
      "y",
      "threat_composite"
    )
  )
  cat(sprintf("  Matched observations: %s\n", format(nrow(d_match), big.mark = ",")))

  template <- rast("data/store/pathreat.data.rasterbase.tif")
  d_match$cell_id <- cellFromXY(template, cbind(d_match$x, d_match$y))
  d_match$x <- NULL
  d_match$y <- NULL

  pair_index <- build_matched_pair_index(d_match)
  outcome_mask <- matched_pair_sample_mask(d_match, "threat_composite", pair_index)
  treated_row_idx <- pair_index$treated_rows
  eligible_pairs <- outcome_mask[treated_row_idx] & !is.na(d_match$cell_id[treated_row_idx])

  estimate_species_delta_s <- function(i) {
    sid <- sp_meta$species_id[i]
    raster_file <- file.path(species_raster_dir, sp_meta$raster_file[i])
    tmp_file <- file.path(tmp_dir, paste0("sp_", sid, ".rds"))

    out <- data.frame(
      species_id = sid,
      delta_s = NA_real_,
      se_delta_s = NA_real_,
      p_delta_s = NA_real_,
      t_stat = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      adjusted_control = NA_real_,
      adjusted_protected = NA_real_,
      pct_effect_s = NA_real_,
      control_mean_s = NA_real_,
      n_pairs = 0L,
      n_obs = 0L,
      n_countries = 0L,
      n_biomes = 0L,
      fe_spec = NA_character_,
      se_type = NA_character_,
      status = NA_character_,
      fit_engine = NA_character_,
      inference_reason = NA_character_,
      row.names = NULL
    )

    if (!file.exists(raster_file)) {
      out$status <- "missing_raster"
      save_rds_atomic(out, tmp_file)
      return(NULL)
    }

    sp_rast <- rast(raster_file)
    range_values <- sp_rast[d_match$cell_id[treated_row_idx]][, 1]
    keep <- eligible_pairs & !is.na(range_values) & range_values == 1
    n_pairs <- sum(keep)
    if (n_pairs == 0L) {
      out$status <- "no_pairs_in_range"
      save_rds_atomic(out, tmp_file)
      return(NULL)
    }

    # The protected range membership selects a complete pair. Sorting the
    # original row indices preserves the previous regression row order and
    # every repeated donor appearance without re-grouping each species.
    selected_rows <- sort(c(treated_row_idx[keep], pair_index$control_rows[keep]))
    s <- d_match[selected_rows, , drop = FALSE]
    n_countries <- n_distinct(s$country_rast[s$treat == 1], na.rm = TRUE)
    n_biomes <- n_distinct(s$biome_raster[s$treat == 1], na.rm = TRUE)
    ctrl_mean <- mean(s$threat_composite[s$treat == 0])

    out$n_pairs <- as.integer(n_pairs)
    out$n_obs <- as.integer(nrow(s))
    out$n_countries <- as.integer(n_countries)
    out$n_biomes <- as.integer(n_biomes)
    out$control_mean_s <- ctrl_mean

    dd <- s[, c("threat_composite", "treat", "country_rast", "biome_raster")]
    fitted <- tryCatch(fit_pair_subgroup(dd), error = function(err) err)
    if (inherits(fitted, "error")) {
      out$status <- paste0("error: ", conditionMessage(fitted))
      save_rds_atomic(out, tmp_file)
      return(NULL)
    }
    out$delta_s <- fitted$coef
    out$se_delta_s <- fitted$se_delta
    out$p_delta_s <- fitted$p_value
    out$t_stat <- fitted$t_stat
    out$ci_low <- fitted$ci_low
    out$ci_high <- fitted$ci_high
    out$adjusted_control <- fitted$adjusted_control
    out$adjusted_protected <- fitted$adjusted_protected
    out$pct_effect_s <- ifelse(
      !is.finite(ctrl_mean) || ctrl_mean == 0,
      NA_real_, out$delta_s / ctrl_mean * 100
    )
    out$fe_spec <- fitted$fe_spec
    out$se_type <- if (n_pairs == 1L) "none" else fitted$se_type
    out$status <- fitted$status
    out$fit_engine <- fitted$fit_engine
    out$inference_reason <- fitted$inference_reason

    save_rds_atomic(out, tmp_file)
    return(NULL)
  }

  cat(sprintf("  Running %s uncached species on %d cores...\n",
              format(length(todo_idx), big.mark = ","), no_cores))
  parallel::mclapply(
    todo_idx,
    estimate_species_delta_s,
    mc.cores = no_cores,
    mc.preschedule = FALSE
  )
}

############################
### collect + validate #####
############################
{
cache_files <- list.files(
  tmp_dir,
  pattern = "^sp_[0-9]+\\.rds$",
  full.names = TRUE
)
cache_rows <- lapply(cache_files, readRDS)

if (length(cache_rows) == 0) {
  stop("No canonical per-species cache rows found in ", tmp_dir)
}
if (any(vapply(cache_rows, nrow, integer(1)) != 1L)) {
  stop("Every canonical species cache file must contain exactly one row.")
}

sp_est <- do.call(rbind, cache_rows)
rownames(sp_est) <- NULL

required_est <- c(
  "species_id",
  "delta_s",
  "se_delta_s",
  "p_delta_s",
  "t_stat",
  "ci_low",
  "ci_high",
  "adjusted_control",
  "adjusted_protected",
  "pct_effect_s",
  "control_mean_s",
  "n_pairs",
  "n_obs",
  "n_countries",
  "n_biomes",
  "fe_spec",
  "se_type",
  "status",
  "fit_engine",
  "inference_reason"
)
if (!identical(names(sp_est), required_est)) {
  stop("Canonical species cache schema does not match the required field order.")
}
if (anyDuplicated(sp_est$species_id)) {
  stop("Canonical species cache contains duplicated species_id values.")
}
if (!setequal(sp_est$species_id, sp_meta$species_id)) {
  stop("Canonical species cache IDs do not match species metadata IDs.")
}
pair_complete <- !is.na(sp_est$n_obs) &
  !is.na(sp_est$n_pairs) &
  sp_est$n_obs == 2L * sp_est$n_pairs
if (any(!pair_complete)) {
  stop("Canonical species cache contains non-pair-complete rows.")
}

# A converged complete-pair regression reproduces both observed group means
# after standardization. Reject invalid predictions before publishing the cache.
nonempty <- sp_est$n_pairs > 0L
if (any(!is.finite(sp_est$adjusted_control[nonempty])) ||
    any(!is.finite(sp_est$adjusted_protected[nonempty])) ||
    any(abs(sp_est$adjusted_control[nonempty] -
            sp_est$control_mean_s[nonempty]) > 1e-8) ||
    any(abs(sp_est$adjusted_protected[nonempty] -
            sp_est$adjusted_control[nonempty] - sp_est$delta_s[nonempty]) > 1e-8)) {
  stop("Canonical species regression prediction identities failed")
}
cat("Species regression engines:\n")
print(table(sp_est$fit_engine, useNA = "ifany"))

sp_result <- sp_meta %>%
  select(species_id, species_name, taxon, category) %>%
  left_join(sp_coverage, by = "species_id") %>%
  left_join(sp_est, by = "species_id") %>%
  arrange(taxon, species_name)

status_counts <- table(sp_result$status, useNA = "ifany")
cat("\nSpecies status counts:\n")
print(status_counts)

if (nrow(sp_result) != 7279L || n_distinct(sp_result$species_id) != 7279L) {
  stop("Expected 7,279 unique species rows.")
}

expected_status_counts <- c(
  estimated = 3310L,
  single_pair = 76L,
  no_pairs_in_range = 3893L
)
unexpected_status <- setdiff(
  unique(sp_result$status[!is.na(sp_result$status)]),
  names(expected_status_counts)
)
if (anyNA(sp_result$status) || length(unexpected_status) > 0L) {
  stop(
    "Unexpected species statuses: ",
    paste(c("NA"[anyNA(sp_result$status)], unexpected_status), collapse = ", ")
  )
}

observed_status_counts <- vapply(
  names(expected_status_counts),
  function(status_name) sum(sp_result$status == status_name),
  integer(1)
)
if (!identical(observed_status_counts, expected_status_counts)) {
  stop(
    "Species status counts do not match the established contract. Observed: ",
    paste(names(observed_status_counts), observed_status_counts, collapse = ", ")
  )
}

estimated_numeric_fields <- c(
  "delta_s",
  "se_delta_s",
  "p_delta_s",
  "t_stat",
  "ci_low",
  "ci_high",
  "adjusted_control",
  "adjusted_protected",
  "pct_effect_s",
  "control_mean_s"
)
estimated_rows <- sp_result$status == "estimated"
estimated_numeric_finite <- vapply(
  estimated_numeric_fields,
  function(field) all(is.finite(sp_result[[field]][estimated_rows])),
  logical(1)
)
if (!all(estimated_numeric_finite)) {
  stop(
    "Estimated species contain nonfinite model fields: ",
    paste(names(estimated_numeric_finite)[!estimated_numeric_finite],
          collapse = ", ")
  )
}

expected_fe_counts <- c(
  `country+biome` = 642L,
  country = 117L,
  biome = 1002L,
  none = 1549L
)
estimated_fe <- sp_result$fe_spec[sp_result$status == "estimated"]
unexpected_fe <- setdiff(
  unique(estimated_fe[!is.na(estimated_fe)]),
  names(expected_fe_counts)
)
if (anyNA(estimated_fe) || length(unexpected_fe) > 0L) {
  stop(
    "Unexpected fixed-effect specifications among estimated species: ",
    paste(c("NA"[anyNA(estimated_fe)], unexpected_fe), collapse = ", ")
  )
}

observed_fe_counts <- vapply(
  names(expected_fe_counts),
  function(fe_name) sum(estimated_fe == fe_name),
  integer(1)
)
if (!identical(observed_fe_counts, expected_fe_counts)) {
  stop(
    "Estimated-species FE counts do not match the established contract. Observed: ",
    paste(names(observed_fe_counts), observed_fe_counts, collapse = ", ")
  )
}

cat("Estimated-species fixed-effect counts:\n")
print(observed_fe_counts)
}

############################
### save ###################
############################
{
coverage_fingerprint <- fingerprint_file(species_coverage_path)
save_rds_atomic(sp_result, outpath)
cache_manifest_with_output <- cache_manifest_current
cache_manifest_with_output$output_provenance <- list(
  output = fingerprint_file(outpath),
  coverage_rds = coverage_fingerprint,
  created_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
)
save_rds_atomic(cache_manifest_with_output, manifest_path)
cat(sprintf("\nSaved: %s (%s rows)\n",
            outpath, format(nrow(sp_result), big.mark = ",")))
}
