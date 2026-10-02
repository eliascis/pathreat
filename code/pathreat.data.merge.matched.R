library(dplyr)
library(fst)
library(ggplot2)

source("code/pathreat.analysis.config.R")

# Load the matching producer's helper-only contract into an isolated
# environment so assembly validates the exact same version, specifications,
# and approved drop set without running matching.
matching_contract_env <- new.env(parent = .GlobalEnv)
matching_helper_setting <- Sys.getenv(
  "PATHREAT_MATCHING_HELPERS_ONLY",
  unset = NA_character_
)
Sys.setenv(PATHREAT_MATCHING_HELPERS_ONLY = "1")
sys.source(
  "code/pathreat.data.matching.R",
  envir = matching_contract_env
)
if (is.na(matching_helper_setting)) {
  Sys.unsetenv("PATHREAT_MATCHING_HELPERS_ONLY")
} else {
  Sys.setenv(PATHREAT_MATCHING_HELPERS_ONLY = matching_helper_setting)
}
approved_manifest_schema <-
  matching_contract_env$matching_manifest_schema_version
approved_contract_version <- matching_contract_env$matching_contract_version
approved_expected_strata <- matching_contract_env$matching_expected_strata
approved_drop_set <- matching_contract_env$matching_expected_drop_set
approved_matching_specification <-
  matching_contract_env$matching_specification()
rm(matching_helper_setting)


#################################
### paths and atomic writers ####
#################################

env_or_default <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) value else default
}

manifest_path <- env_or_default(
  "PATHREAT_MATCHING_MANIFEST",
  "data/store/pathreat.data.matching.manifest.Rds"
)
match_dir <- "data/store/newmatchresults"
matched_output_path <- paths$data_matched
balance_output_path <- file.path(
  paths$results_dir,
  "pathreat.covbalance.by_country.Rds"
)
summary_output_path <- file.path(
  paths$results_dir,
  "pathreat.matched-sample.summary.csv"
)
balance_figure_path <- file.path(
  paths$figures_dir,
  "fig.match.balance.by_country.jpg"
)
balance_threshold <- 0.10

atomic_rename <- function(temporary_path, output_path) {
  if (!file.rename(temporary_path, output_path)) {
    stop("Atomic rename failed: ", temporary_path, " -> ", output_path)
  }
  invisible(output_path)
}

atomic_save_rds <- function(object, output_path, compress = TRUE) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- tempfile(
    pattern = paste0(".", basename(output_path), "."),
    tmpdir = dirname(output_path),
    fileext = ".tmp"
  )
  on.exit(unlink(temporary_path), add = TRUE)
  saveRDS(object, temporary_path, compress = compress)
  atomic_rename(temporary_path, output_path)
}

atomic_write_csv <- function(object, output_path) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- tempfile(
    pattern = paste0(".", basename(output_path), "."),
    tmpdir = dirname(output_path),
    fileext = ".tmp"
  )
  on.exit(unlink(temporary_path), add = TRUE)
  write.csv(object, temporary_path, row.names = FALSE, na = "")
  atomic_rename(temporary_path, output_path)
}

atomic_write_fst <- function(object, output_path) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- tempfile(
    pattern = paste0(".", basename(output_path), "."),
    tmpdir = dirname(output_path),
    fileext = ".fst"
  )
  on.exit(unlink(temporary_path), add = TRUE)
  write_fst(object, temporary_path)
  check <- tryCatch(fst::fst(temporary_path), error = identity)
  if (inherits(check, "error") || nrow(check) != nrow(object)) {
    stop("Candidate matched FST failed its read-back validation.")
  }
  atomic_rename(temporary_path, output_path)
}

atomic_save_plot <- function(plot, output_path, width, height, dpi = 300) {
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  temporary_path <- tempfile(
    pattern = paste0(".", tools::file_path_sans_ext(basename(output_path)), "."),
    tmpdir = dirname(output_path),
    fileext = ".jpg"
  )
  on.exit(unlink(temporary_path), add = TRUE)
  ggsave(
    filename = temporary_path,
    plot = plot,
    units = "cm",
    width = width,
    height = height,
    dpi = dpi
  )
  if (!file.exists(temporary_path) || file.info(temporary_path)$size <= 0) {
    stop("Candidate balance figure is missing or empty.")
  }
  atomic_rename(temporary_path, output_path)
}


#################################
### generic validation ##########
#################################

require_fields <- function(object, required, label) {
  missing <- setdiff(required, names(object))
  if (length(missing) > 0L) {
    stop(label, " is missing fields: ", paste(missing, collapse = ", "))
  }
}

read_rds_checked <- function(path, label) {
  if (!file.exists(path)) stop("Missing ", label, ": ", path)
  object <- tryCatch(readRDS(path), error = identity)
  if (inherits(object, "error")) {
    stop("Unreadable ", label, " (", path, "): ", conditionMessage(object))
  }
  object
}

canonical_path <- function(path, must_work = TRUE) {
  normalizePath(path, winslash = "/", mustWork = must_work)
}

file_md5 <- function(path) {
  unname(tools::md5sum(path))
}


#################################
### validate matching manifest ##
#################################

cat("=== Validating matching manifest ===\n")
manifest <- read_rds_checked(manifest_path, "matching manifest")
require_fields(
  manifest,
  c(
    "schema_version", "contract_version", "created_at_utc", "input",
    "specification", "expected_strata", "expected_drop_set", "strata",
    "complete", "errors"
  ),
  "Matching manifest"
)
if (!identical(manifest$schema_version, approved_manifest_schema) ||
    !identical(manifest$contract_version, approved_contract_version)) {
  stop("Unsupported matching-manifest schema or contract version.")
}
if (!isTRUE(manifest$complete) || length(manifest$errors) > 0L) {
  stop("Matching manifest is incomplete: ", paste(manifest$errors, collapse = "; "))
}
if (!identical(as.integer(manifest$expected_strata), approved_expected_strata)) {
  stop("Matching manifest declares the wrong number of expected strata.")
}

require_fields(manifest$input, c("path", "fingerprint"), "Manifest input")
input_fingerprint <- manifest$input$fingerprint
require_fields(
  input_fingerprint,
  c(
    "algorithm", "value", "size_bytes", "mtime_unix", "n_rows", "columns"
  ),
  "Manifest input fingerprint"
)
if (!identical(input_fingerprint$algorithm, "md5") ||
    !grepl("^[0-9a-f]{32}$", input_fingerprint$value)) {
  stop("Manifest unmatched-input fingerprint is malformed.")
}
if (canonical_path(manifest$input$path) != canonical_path(paths$data_unmatched)) {
  stop("Manifest unmatched-input path does not match configured unmatched data.")
}

unmatched_info <- file.info(paths$data_unmatched)
unmatched_metadata <- fst::fst(paths$data_unmatched)
if (as.numeric(unmatched_info$size) != input_fingerprint$size_bytes ||
    abs(as.numeric(unmatched_info$mtime) - input_fingerprint$mtime_unix) > 1e-6 ||
    as.numeric(nrow(unmatched_metadata)) != input_fingerprint$n_rows ||
    !identical(colnames(unmatched_metadata), input_fingerprint$columns) ||
    !identical(file_md5(paths$data_unmatched), input_fingerprint$value)) {
  stop("Matching manifest is stale relative to the unmatched input.")
}
rm(unmatched_info, unmatched_metadata)

specification <- manifest$specification
require_fields(
  specification,
  c(
    "covariates", "near_constant_rule", "match", "filters", "strata",
    "big_countries", "pooling"
  ),
  "Manifest matching specification"
)
if (!identical(specification, approved_matching_specification)) {
  stop("Manifest matching specification differs from the approved v3 contract.")
}
if (!identical(specification$covariates, mlist) ||
    !identical(
      specification$strata,
      c("country_rast", "biome_raster", "bigregion")
    )) {
  stop("Manifest covariates or stratum definition differ from configuration.")
}
require_fields(
  specification$near_constant_rule,
  c("min_unique", "relative_sd_threshold", "formula"),
  "Manifest near-constant rule"
)
if (!identical(as.integer(specification$near_constant_rule$min_unique), 2L) ||
    !identical(
      as.numeric(specification$near_constant_rule$relative_sd_threshold),
      1e-6
    ) ||
    !identical(
      specification$near_constant_rule$formula,
      "sd(x) / max(1, abs(mean(x)))"
    )) {
  stop("Manifest near-constant rule differs from the approved contract.")
}
require_fields(
  specification$match,
  c("function", "M", "Weight", "replace", "ties", "caliper"),
  "Manifest Match specification"
)
if (!identical(specification$match[["function"]], "Matching::Match") ||
    as.integer(specification$match$M) != 1L ||
    as.integer(specification$match$Weight) != 2L ||
    !isTRUE(specification$match$replace) ||
    !identical(specification$match$ties, FALSE) ||
    !is.null(specification$match$caliper)) {
  stop("Manifest Match parameters differ from the approved contract.")
}

expected_drop_set <- manifest$expected_drop_set
if (!is.data.frame(expected_drop_set) ||
    !identical(names(expected_drop_set), c("stratum_id", "covariate")) ||
    nrow(expected_drop_set) == 0L ||
    anyNA(expected_drop_set) || anyDuplicated(expected_drop_set)) {
  stop("Manifest expected-drop set is malformed.")
}
drop_contract_key <- function(data) {
  sort(paste(data$stratum_id, data$covariate, sep = "::"))
}
if (!identical(
      drop_contract_key(expected_drop_set),
      drop_contract_key(approved_drop_set)
    ) ||
    nrow(expected_drop_set) != 67L ||
    length(unique(expected_drop_set$stratum_id)) != 62L) {
  stop("Manifest expected-drop set differs from the approved v3 contract.")
}

strata <- manifest$strata
if (!is.data.frame(strata)) stop("Manifest strata must be a data frame.")
stratum_fields <- c(
  "stratum_id", "country_rast", "country", "biome_raster", "bigregion",
  "source_treated", "source_controls_stratum", "source_controls_effective",
  "source_rows_effective", "pooling_scope", "output_file", "pair_count",
  "status", "error_text", "action", "validation_status", "output_md5",
  "reuse_reason", "output_size_bytes", "selected_covariates", "dropped_covariates",
  "dropped_diagnostics", "validation_errors"
)
require_fields(strata, stratum_fields, "Manifest strata")
if (nrow(strata) != approved_expected_strata ||
    anyDuplicated(strata$stratum_id) ||
    anyDuplicated(strata$output_file)) {
  stop("Manifest strata must contain 668 unique strata and output files.")
}
if (anyNA(strata$stratum_id) || anyNA(strata$output_file) ||
    any(strata$status != "success") ||
    any(!strata$action %in% c("reused", "rematched")) ||
    any(!strata$validation_status %in% c("valid_reused", "new_output_valid")) ||
    anyNA(strata$reuse_reason) || any(!nzchar(strata$reuse_reason)) ||
    anyNA(strata$pair_count) || any(strata$pair_count <= 0L) ||
    any(strata$pair_count != strata$source_treated)) {
  stop("Manifest contains an unsuccessful or internally inconsistent stratum.")
}
if (any((strata$action == "reused") !=
        (strata$validation_status == "valid_reused")) ||
    any((strata$action == "rematched") !=
        (strata$validation_status == "new_output_valid"))) {
  stop("Manifest action and validation-status fields disagree.")
}
affected_strata <- approved_drop_set$stratum_id
reused_affected <- strata$action == "reused" &
  strata$stratum_id %in% affected_strata
reused_unaffected <- strata$action == "reused" & !reused_affected
rematched <- strata$action == "rematched"
if (any(strata$reuse_reason[reused_unaffected] !=
        "valid_unaffected_output") ||
    any(strata$reuse_reason[reused_affected] !=
        "valid_affected_output_attested_by_prior_complete_manifest") ||
    any(!grepl(
      "^(affected_output_unattested:|invalid_or_missing_output:)",
      strata$reuse_reason[rematched]
    ))) {
  stop("Manifest reuse reasons disagree with v3 action semantics.")
}

identifier_fields <- c(
  "country_rast", "biome_raster", "bigregion", "source_treated",
  "source_controls_stratum", "source_controls_effective",
  "source_rows_effective", "pair_count"
)
invalid_identifier <- vapply(identifier_fields, function(field) {
  values <- strata[[field]]
  !is.numeric(values) || anyNA(values) || any(!is.finite(values)) ||
    any(values != floor(values))
}, logical(1))
if (any(invalid_identifier) ||
    any(strata$country_rast <= 0L) || any(strata$biome_raster <= 0L) ||
    any(strata$bigregion < 0L) || any(strata$source_treated <= 0L) ||
    any(strata$source_controls_stratum < 0L) ||
    any(strata$source_controls_effective <= 0L)) {
  stop("Manifest identifiers or source counts have invalid storage or values.")
}
if (any(strata$stratum_id != paste(
      strata$country_rast,
      strata$biome_raster,
      strata$bigregion,
      sep = "."
    )) ||
    anyNA(strata$country) || any(!nzchar(strata$country))) {
  stop("Manifest stratum IDs or country names are malformed.")
}
if (any(!strata$pooling_scope %in% c("none", "country", "country_biome")) ||
    any(strata$source_controls_effective < strata$source_controls_stratum) ||
    any(strata$pooling_scope == "none" &
          strata$source_controls_effective != strata$source_controls_stratum) ||
    any(strata$pooling_scope == "none" &
          strata$source_controls_stratum < 20L) ||
    any(strata$pooling_scope != "none" &
          strata$source_controls_stratum >= 20L) ||
    any(strata$pooling_scope == "country" & strata$bigregion != 0L) ||
    any(strata$pooling_scope == "country_biome" & strata$bigregion <= 0L)) {
  stop("Manifest pooling scopes disagree with source counts or bigregion IDs.")
}
source_count_fields <- c(
  "source_treated", "source_controls_stratum", "source_controls_effective",
  "source_rows_effective", "pair_count"
)
invalid_source_count <- vapply(source_count_fields, function(field) {
  values <- strata[[field]]
  !is.numeric(values) || anyNA(values) || any(!is.finite(values)) ||
    any(values != floor(values)) || any(values < 0)
}, logical(1))
if (any(invalid_source_count) ||
    any(strata$source_treated <= 0L) ||
    any(strata$source_controls_effective <= 0L) ||
    any(strata$source_controls_stratum > strata$source_controls_effective) ||
    any(strata$source_rows_effective !=
          strata$source_treated + strata$source_controls_effective)) {
  stop("Manifest source counts are malformed or internally inconsistent.")
}
if (anyNA(strata$output_md5) ||
    any(!grepl("^[0-9a-f]{32}$", strata$output_md5)) ||
    anyNA(strata$output_size_bytes) || any(strata$output_size_bytes <= 0)) {
  stop("Manifest contains malformed output fingerprints.")
}

actual_drop_keys <- unlist(lapply(seq_len(nrow(strata)), function(i) {
  dropped <- strata$dropped_covariates[[i]]
  if (length(dropped) == 0L) return(character())
  paste(strata$stratum_id[i], dropped, sep = "::")
}), use.names = FALSE)
expected_drop_keys <- paste(
  expected_drop_set$stratum_id,
  expected_drop_set$covariate,
  sep = "::"
)
if (!identical(sort(actual_drop_keys), sort(expected_drop_keys))) {
  stop(
    "Manifest recorded drop set differs from its expected-drop contract."
  )
}

diagnostic_fields <- c(
  "covariate", "n_unique", "min", "max", "mean", "sd", "relative_sd",
  "drop"
)
for (i in seq_len(nrow(strata))) {
  selected <- strata$selected_covariates[[i]]
  dropped <- strata$dropped_covariates[[i]]
  diagnostics <- strata$dropped_diagnostics[[i]]
  validation_errors <- strata$validation_errors[[i]]
  if (!is.character(selected) || !is.character(dropped) ||
      length(selected) == 0L || anyNA(selected) || anyNA(dropped) ||
      anyDuplicated(selected) || anyDuplicated(dropped) ||
      length(intersect(selected, dropped)) > 0L ||
      !setequal(c(selected, dropped), mlist)) {
    stop("Manifest selected/dropped covariates are malformed for ",
         strata$stratum_id[i], ".")
  }
  if (!is.character(validation_errors) || length(validation_errors) != 0L ||
      !is.na(strata$error_text[i])) {
    stop("Successful manifest stratum retains validation errors: ",
         strata$stratum_id[i], ".")
  }
  if (!is.data.frame(diagnostics) ||
      !all(diagnostic_fields %in% names(diagnostics)) ||
      nrow(diagnostics) != length(dropped) ||
      !setequal(diagnostics$covariate, dropped)) {
    stop("Manifest drop diagnostics are malformed for ",
         strata$stratum_id[i], ".")
  }
  if (nrow(diagnostics) > 0L &&
      (anyNA(diagnostics[, diagnostic_fields, drop = FALSE]) ||
       any(diagnostics$n_unique < 1L) ||
       any(!diagnostics$drop) ||
       any(!(diagnostics$n_unique < 2L |
             diagnostics$relative_sd <= 1e-6)))) {
    stop("Manifest drop diagnostics violate the near-constant rule for ",
         strata$stratum_id[i], ".")
  }
}


#################################
### reconstruct source contract #
#################################

# Do not accept a manifest and set of match files merely because they agree
# with one another. Reapply the matching producer's filters to the fingerprinted
# unmatched input, independently rebuild every expected effective stratum, and
# retain the resulting source state long enough to validate every matched row.
cat("=== Reconstructing effective matching strata from source ===\n")
contract_input <- matching_contract_env$load_matching_input(
  paths$data_unmatched
)
contract_source_strata <- matching_contract_env$build_expected_strata(
  contract_input
)
contract_state <- matching_contract_env$build_matching_state(contract_input)
rm(contract_input)
gc()

contract_prepared <- matching_contract_env$prepare_strata_contract(
  contract_source_strata,
  contract_state
)
contract_source_strata <- contract_prepared$source_strata
contract_manifest_strata <- contract_prepared$manifest_strata
rm(contract_prepared)

contract_stratum_index <- match(
  strata$stratum_id,
  contract_manifest_strata$stratum_id
)
if (anyNA(contract_stratum_index) ||
    anyDuplicated(contract_stratum_index) ||
    length(contract_stratum_index) != nrow(contract_manifest_strata)) {
  stop("Manifest strata differ from independently reconstructed strata.")
}
contract_manifest_ordered <-
  contract_manifest_strata[contract_stratum_index, , drop = FALSE]
contract_scalar_fields <- c(
  "stratum_id", "country_rast", "country", "biome_raster", "bigregion",
  "source_treated", "source_controls_stratum", "source_controls_effective",
  "source_rows_effective", "pooling_scope", "output_file"
)
for (field in contract_scalar_fields) {
  if (!identical(
    strata[[field]],
    contract_manifest_ordered[[field]]
  )) {
    stop(
      "Manifest field disagrees with reconstructed source contract: ",
      field,
      "."
    )
  }
}
for (i in seq_len(nrow(strata))) {
  if (!identical(
        strata$selected_covariates[[i]],
        contract_manifest_ordered$selected_covariates[[i]]
      ) ||
      !identical(
        strata$dropped_covariates[[i]],
        contract_manifest_ordered$dropped_covariates[[i]]
      )) {
    stop(
      "Manifest covariate selection disagrees with reconstructed source ",
      "contract for ",
      strata$stratum_id[i],
      "."
    )
  }
}
rm(contract_manifest_ordered)
cat("  Source reconstruction: 668 exact effective strata.\n")

expected_files <- vapply(
  strata$output_file,
  canonical_path,
  character(1),
  must_work = TRUE
)
actual_files_raw <- list.files(
  match_dir,
  pattern = "\\.Rds$",
  full.names = TRUE
)
actual_files <- vapply(
  actual_files_raw,
  canonical_path,
  character(1),
  must_work = TRUE
)
if (anyDuplicated(expected_files) || anyDuplicated(actual_files) ||
    length(expected_files) != length(actual_files) ||
    !setequal(expected_files, actual_files)) {
  missing_files <- setdiff(expected_files, actual_files)
  extra_files <- setdiff(actual_files, expected_files)
  stop(
    "Match-result file set disagrees with manifest. Missing: ",
    paste(basename(missing_files), collapse = ", "),
    "; extra: ", paste(basename(extra_files), collapse = ", ")
  )
}
cat("  Manifest: 668 complete strata; input fingerprint current.\n")


#################################
### validate exact match files ##
#################################

cat("=== Reading manifest-authorized match files ===\n")
match_columns <- c(
  "country_rast", "biome_raster", "row_id", "pair_id",
  "weight", "pixel_id", "treat", "nid"
)

validate_and_read_match <- function(i) {
  row <- strata[i, , drop = FALSE]
  source_i <- contract_stratum_index[i]
  source_g <- contract_source_strata[source_i, , drop = FALSE]
  source_effective <- matching_contract_env$effective_stratum_rows(
    source_g,
    contract_state
  )
  path <- row$output_file
  info <- file.info(path)
  if (as.numeric(info$size) != row$output_size_bytes ||
      !identical(file_md5(path), row$output_md5)) {
    stop("Stale match-result fingerprint for stratum ", row$stratum_id, ".")
  }

  result <- read_rds_checked(path, paste0("match result ", row$stratum_id))
  source_validation <- matching_contract_env$validate_match_result(
    result,
    source_g,
    contract_state,
    source_effective
  )
  if (!source_validation$valid) {
    stop(
      "Match result disagrees with reconstructed effective source for ",
      row$stratum_id,
      ": ",
      paste(source_validation$errors, collapse = "; ")
    )
  }
  if (!is.data.frame(result) || !identical(names(result), match_columns)) {
    stop("Malformed match-result schema for stratum ", row$stratum_id, ".")
  }
  if (any(!vapply(result[, match_columns, drop = FALSE], is.numeric, logical(1)))) {
    stop("Match-result columns have invalid storage types for stratum ",
         row$stratum_id, ".")
  }
  pair_count <- as.integer(row$pair_count)
  if (nrow(result) != 2L * pair_count ||
      anyNA(result$treat) || !all(result$treat %in% c(0, 1)) ||
      sum(result$treat == 1) != pair_count ||
      sum(result$treat == 0) != pair_count) {
    stop("Treatment/control parity failed for stratum ", row$stratum_id, ".")
  }
  if (anyNA(result$country_rast) ||
      any(result$country_rast != row$country_rast) ||
      anyNA(result$biome_raster) ||
      any(result$biome_raster != row$biome_raster)) {
    stop("Match contents disagree with manifest identifiers for ", row$stratum_id, ".")
  }
  if (anyNA(result$pixel_id) || any(!is.finite(result$pixel_id)) ||
      any(result$pixel_id <= 0) ||
      anyDuplicated(result$pixel_id[result$treat == 1])) {
    stop("Invalid or reused treated pixel ID in stratum ", row$stratum_id, ".")
  }
  if (anyNA(result$weight) || any(!is.finite(result$weight)) ||
      any(result$weight <= 0)) {
    stop("Invalid match weights in stratum ", row$stratum_id, ".")
  }
  integer_fields <- c(
    "country_rast", "biome_raster", "row_id", "pair_id", "pixel_id", "nid"
  )
  invalid_integer <- vapply(integer_fields, function(field) {
    values <- result[[field]]
    anyNA(values) || any(!is.finite(values)) ||
      any(values != floor(values)) || any(values <= 0)
  }, logical(1))
  if (any(invalid_integer) ||
      !identical(as.integer(result$nid), seq_len(nrow(result)))) {
    stop("Invalid row, pair, or sequence IDs in stratum ", row$stratum_id, ".")
  }

  treated_rows <- seq_len(pair_count)
  control_rows <- pair_count + seq_len(pair_count)
  if (!all(result$treat[treated_rows] == 1L) ||
      !all(result$treat[control_rows] == 0L) ||
      !all(result$row_id[treated_rows] == result$pair_id[control_rows]) ||
      !all(result$pair_id[treated_rows] == result$row_id[control_rows]) ||
      !all(result$weight[treated_rows] == result$weight[control_rows])) {
    stop("Pair cross-references failed for stratum ", row$stratum_id, ".")
  }

  result$bigregion <- as.integer(row$bigregion)
  result$stratum_id <- row$stratum_id
  result$.pooling_scope <- as.character(row$pooling_scope)
  result
}

manifest_read_order <- match(actual_files, expected_files)
if (anyNA(manifest_read_order) || anyDuplicated(manifest_read_order)) {
  stop("Could not map alphabetic match-file order to manifest rows.")
}
match_list <- lapply(manifest_read_order, validate_and_read_match)
d_matched_ids <- do.call(rbind, match_list)
rownames(d_matched_ids) <- NULL
rm(
  match_list,
  manifest_read_order,
  validate_and_read_match,
  contract_source_strata,
  contract_state,
  contract_manifest_strata,
  contract_stratum_index,
  matching_contract_env
)
gc()

if (sum(strata$pair_count) * 2 != nrow(d_matched_ids)) {
  stop("Combined match-result row count disagrees with manifest pair counts.")
}
indonesia_manifest_pairs <- sum(
  strata$pair_count[strata$country_rast == 102L]
)
indonesia_missing_stratum_pairs <- strata$pair_count[
  strata$stratum_id == "102.7.0"
]
if (length(indonesia_missing_stratum_pairs) != 1L ||
    indonesia_missing_stratum_pairs != 3645L ||
    indonesia_manifest_pairs != 3727L) {
  stop(
    "Indonesia manifest gate failed: stratum 102.7.0 must have 3,645 pairs ",
    "and Indonesia must total 3,727 pairs."
  )
}
cat(
  "  Valid match files: ", nrow(strata),
  "; Indonesia pairs: ", format(indonesia_manifest_pairs, big.mark = ","),
  ".\n",
  sep = ""
)


#################################
### construct global pair IDs ###
#################################

d <- d_matched_ids %>%
  rename(paired_row_id = pair_id)
d$pair_key <- ifelse(
  d$treat == 1L,
  paste(
    d$country_rast, d$biome_raster, d$bigregion,
    d$row_id, d$paired_row_id,
    sep = "_"
  ),
  paste(
    d$country_rast, d$biome_raster, d$bigregion,
    d$paired_row_id, d$row_id,
    sep = "_"
  )
)
d$matched_pair_id <- as.integer(factor(d$pair_key))
d$pair_key <- NULL

pair_check <- d %>%
  group_by(matched_pair_id) %>%
  summarize(
    n_obs = n(),
    n_treated = sum(treat == 1L),
    n_control = sum(treat == 0L),
    n_strata = n_distinct(stratum_id),
    .groups = "drop"
  )
if (any(pair_check$n_obs != 2L) ||
    any(pair_check$n_treated != 1L) ||
    any(pair_check$n_control != 1L) ||
    any(pair_check$n_strata != 1L) ||
    nrow(pair_check) != sum(strata$pair_count)) {
  stop("Global matched-pair construction is incomplete or duplicated.")
}
rm(pair_check, d_matched_ids)
gc()


#################################
### join unmatched attributes ###
#################################

cat("=== Loading unmatched attributes for matched pixels ===\n")
matched_pixel_ids <- unique(d$pixel_id)
unmatched_table <- fst::fst(paths$data_unmatched)
total_rows <- nrow(unmatched_table)
chunk_size <- 10000000L
chunks <- list()
for (start in seq.int(1L, total_rows, by = chunk_size)) {
  end <- min(start + chunk_size - 1L, total_rows)
  cat(
    "  Rows ", format(start, big.mark = ","), "-",
    format(end, big.mark = ","), "\n",
    sep = ""
  )
  chunk <- read_fst(paths$data_unmatched, from = start, to = end)
  chunk <- chunk[chunk$pixel_id %in% matched_pixel_ids, , drop = FALSE]
  if (nrow(chunk) > 0L) chunks[[length(chunks) + 1L]] <- chunk
  rm(chunk)
  gc()
}
unmatched_matched <- do.call(rbind, chunks)
rm(chunks, unmatched_table)
gc()

if (anyDuplicated(unmatched_matched$pixel_id) ||
    !setequal(unmatched_matched$pixel_id, matched_pixel_ids)) {
  stop("Matched pixel IDs are missing or duplicated in the unmatched input.")
}
require_fields(
  unmatched_matched,
  c(
    "treat", "country_rast", "country", "biome_raster", "admin1_id",
    "pixel_id"
  ),
  "Unmatched source rows"
)
source_numeric_fields <- c(
  "treat", "country_rast", "biome_raster", "admin1_id", "pixel_id"
)
if (any(!vapply(
  unmatched_matched[, source_numeric_fields, drop = FALSE],
  is.numeric,
  logical(1)
))) {
  stop("Unmatched source identifier columns have invalid storage types.")
}
unmatched_matched$.source_treat <- unmatched_matched$treat
unmatched_matched$.source_country_rast <- unmatched_matched$country_rast
unmatched_matched$.source_biome_raster <- unmatched_matched$biome_raster
unmatched_matched$.source_bigregion <- 0L
source_big_country <- unmatched_matched$country %in%
  specification$big_countries
if (any(source_big_country & is.na(unmatched_matched$admin1_id))) {
  stop("Big-country matched pixels contain missing source admin1 IDs.")
}
unmatched_matched$.source_bigregion[source_big_country] <-
  unmatched_matched$admin1_id[source_big_country]
unmatched_matched$treat <- NULL
unmatched_matched$biome_raster <- NULL
unmatched_matched$country_rast <- NULL
unmatched_matched$.unmatched_join_ok <- TRUE

expected_join_rows <- nrow(d)
d <- merge(
  d,
  unmatched_matched,
  by = "pixel_id",
  all.x = TRUE,
  sort = TRUE
)
if (nrow(d) != expected_join_rows || anyNA(d$.unmatched_join_ok)) {
  stop("Unmatched-attribute join changed row count or left unmatched pixels.")
}
d$.unmatched_join_ok <- NULL
rm(unmatched_matched, matched_pixel_ids, expected_join_rows)
gc()

required_joined <- c(
  "country", "wdpaid", "x", "y", "threat_composite", mlist
)
require_fields(d, required_joined, "Joined matched data")
if (anyNA(d$country_rast) || anyNA(d$biome_raster) ||
    anyNA(d$country) || anyNA(d$x) || anyNA(d$y) ||
    any(!is.finite(d$x)) || any(!is.finite(d$y)) ||
    anyNA(d$threat_composite) || any(!is.finite(d$threat_composite))) {
  stop("Joined matched data contain missing/non-finite identifiers or outcomes.")
}
if (anyNA(d$.source_treat) || anyNA(d$.source_country_rast) ||
    anyNA(d$.source_biome_raster) || anyNA(d$.source_bigregion) ||
    any(d$treat != d$.source_treat) ||
    any(d$country_rast != d$.source_country_rast)) {
  stop(
    "Match-file treatment or country identifiers disagree with source data."
  )
}

treated_rows <- d$treat == 1L
control_rows <- d$treat == 0L
no_pool_controls <- control_rows & d$.pooling_scope == "none"
country_pool_controls <- control_rows & d$.pooling_scope == "country"
country_biome_pool_controls <-
  control_rows & d$.pooling_scope == "country_biome"
if (any(!d$.pooling_scope %in% c("none", "country", "country_biome")) ||
    any(d$biome_raster[treated_rows] != d$.source_biome_raster[treated_rows]) ||
    any(d$bigregion[treated_rows] != d$.source_bigregion[treated_rows]) ||
    any(d$biome_raster[no_pool_controls] !=
          d$.source_biome_raster[no_pool_controls]) ||
    any(d$bigregion[no_pool_controls] !=
          d$.source_bigregion[no_pool_controls]) ||
    any(d$.source_bigregion[country_pool_controls] != 0L) ||
    any(d$biome_raster[country_biome_pool_controls] !=
          d$.source_biome_raster[country_biome_pool_controls])) {
  stop("Match-file rows violate the declared control-pooling scope.")
}
manifest_country_lookup <- unique(
  strata[, c("country_rast", "country"), drop = FALSE]
)
if (anyDuplicated(manifest_country_lookup$country_rast)) {
  stop("Manifest country IDs map to multiple country names.")
}
manifest_country_index <- match(
  d$country_rast,
  manifest_country_lookup$country_rast
)
if (anyNA(manifest_country_index) ||
    any(d$country != manifest_country_lookup$country[manifest_country_index])) {
  stop("Manifest country names disagree with source unmatched data.")
}
rm(manifest_country_lookup, manifest_country_index)
d$.source_treat <- NULL
d$.source_country_rast <- NULL
d$.source_biome_raster <- NULL
d$.source_bigregion <- NULL
d$.pooling_scope <- NULL


#################################
### pair outcomes and checks ####
#################################

cat("=== Computing and validating pair outcomes ===\n")
d <- d %>%
  group_by(matched_pair_id) %>%
  mutate(
    delta = threat_composite[treat == 1L] -
      threat_composite[treat == 0L],
    tc_threat_composite = threat_composite[treat == 0L]
  ) %>%
  ungroup() %>%
  data.frame()
if (anyNA(d$delta) || any(!is.finite(d$delta)) ||
    anyNA(d$tc_threat_composite) ||
    any(!is.finite(d$tc_threat_composite))) {
  stop("Matched-pair outcome construction produced missing/non-finite values.")
}

final_pair_check <- d %>%
  group_by(matched_pair_id) %>%
  summarize(
    n_obs = n(),
    n_treated = sum(treat == 1L),
    n_control = sum(treat == 0L),
    n_delta = n_distinct(delta),
    n_control_outcome = n_distinct(tc_threat_composite),
    .groups = "drop"
  )
if (any(final_pair_check$n_obs != 2L) ||
    any(final_pair_check$n_treated != 1L) ||
    any(final_pair_check$n_control != 1L) ||
    any(final_pair_check$n_delta != 1L) ||
    any(final_pair_check$n_control_outcome != 1L)) {
  stop("Pair completeness failed after the unmatched-data join.")
}
rm(final_pair_check)


#################################
### country covariate balance ###
#################################

cat("=== Evaluating country balance ===\n")
d_unmatched_balance <- read_fst(
  paths$data_unmatched,
  columns = c("treat", "country_rast", mlist)
)
pooled_sd_pre <- lapply(
  sort(unique(d_unmatched_balance$country_rast)),
  function(country_id) {
    country_data <- d_unmatched_balance[
      d_unmatched_balance$country_rast == country_id,
      ,
      drop = FALSE
    ]
    values <- vapply(mlist, function(covariate) {
      sd_treated <- sd(
        country_data[country_data$treat == 1L, covariate],
        na.rm = TRUE
      )
      sd_control <- sd(
        country_data[country_data$treat == 0L, covariate],
        na.rm = TRUE
      )
      sqrt((sd_treated^2 + sd_control^2) / 2)
    }, numeric(1))
    names(values) <- mlist
    values
  }
)
names(pooled_sd_pre) <- sort(unique(d_unmatched_balance$country_rast))
rm(d_unmatched_balance)
gc()

country_ids <- sort(unique(d$country_rast))
balance <- lapply(country_ids, function(country_id) {
  country_rows <- which(d$country_rast == country_id)
  country_treat <- d$treat[country_rows]
  denominators <- pooled_sd_pre[[as.character(country_id)]]
  smds <- vapply(mlist, function(covariate) {
    country_values <- d[[covariate]][country_rows]
    mean_treated <- mean(
      country_values[country_treat == 1L],
      na.rm = TRUE
    )
    mean_control <- mean(
      country_values[country_treat == 0L],
      na.rm = TRUE
    )
    denominator <- denominators[covariate]
    if (length(denominator) == 0L || is.na(denominator) || denominator == 0) {
      if (is.na(mean_treated) || is.na(mean_control)) return(NA_real_)
      if (mean_treated == mean_control) return(0)
      return(NA_real_)
    }
    (mean_treated - mean_control) / denominator
  }, numeric(1))
  data.frame(
    country_id = country_id,
    n_pairs = sum(country_treat == 1L),
    as.list(smds),
    avg_abs_smd = mean(abs(smds), na.rm = TRUE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})
balance <- do.call(rbind, balance)
country_lookup <- unique(d[, c("country_rast", "country")])
if (anyDuplicated(country_lookup$country_rast)) {
  country_name_counts <- country_lookup %>%
    count(country_rast, name = "n_names")
  if (any(country_name_counts$n_names != 1L)) {
    stop("Country raster IDs map to multiple country names.")
  }
  country_lookup <- country_lookup[!duplicated(country_lookup$country_rast), ]
}
names(country_lookup) <- c("country_id", "country_name")
balance <- merge(
  balance,
  country_lookup,
  by = "country_id",
  all.x = TRUE,
  sort = FALSE
) %>%
  relocate(country_name, .after = country_id) %>%
  arrange(desc(avg_abs_smd)) %>%
  data.frame()
if (anyNA(balance$country_name) || any(!is.finite(balance$avg_abs_smd))) {
  stop("Country-balance results contain missing names or non-finite SMDs.")
}
if (any(balance$n_pairs <= 0L) ||
    sum(balance$n_pairs) != sum(d$treat == 1L)) {
  stop("Country-balance pair counts disagree with the matched sample.")
}
rm(country_ids, country_lookup, pooled_sd_pre)
gc()

indonesia_balance <- balance[balance$country_id == 102L, , drop = FALSE]
if (nrow(indonesia_balance) != 1L ||
    indonesia_balance$n_pairs != 3727L) {
  stop("Indonesia balance row is missing or does not contain 3,727 pairs.")
}
indonesia_passes_balance <-
  indonesia_balance$avg_abs_smd <= balance_threshold

balanced_countries <- balance$country_id[
  balance$avg_abs_smd <= balance_threshold
]
d_post_balance <- d[
  d$country_rast %in% balanced_countries,
  ,
  drop = FALSE
]
if (nrow(d_post_balance) == 0L) stop("Country-balance filter removed all data.")

post_pair_check <- d_post_balance %>%
  count(matched_pair_id, treat, name = "n") %>%
  tidyr::pivot_wider(
    names_from = treat,
    values_from = n,
    values_fill = 0,
    names_prefix = "treat_"
  )
if (any(post_pair_check$treat_0 != 1L) ||
    any(post_pair_check$treat_1 != 1L)) {
  stop("Country-balance filter broke matched-pair completeness.")
}
rm(post_pair_check)

indonesia_post_pairs <- sum(
  d_post_balance$country_rast == 102L & d_post_balance$treat == 1L
)
if ((indonesia_passes_balance && indonesia_post_pairs != 3727L) ||
    (!indonesia_passes_balance && indonesia_post_pairs != 0L)) {
  stop("Indonesia retention disagrees with the approved country-balance rule.")
}


#################################
### matched-sample audit ########
#################################

count_cells_25km <- function(data, rows) {
  if (length(rows) == 0L) return(0L)
  as.integer(length(unique(paste(
    round(data$x[rows] / 25000) * 25000,
    round(data$y[rows] / 25000) * 25000,
    sep = "_"
  ))))
}

summarize_stage <- function(data, stage) {
  treated_rows <- which(data$treat == 1L)
  indonesia_rows <- which(data$country_rast == 102L)
  indonesia_treated_rows <- which(
    data$country_rast == 102L & data$treat == 1L
  )
  data.frame(
    stage = stage,
    observations = nrow(data),
    matched_pairs = n_distinct(data$matched_pair_id),
    countries = n_distinct(data$country_rast),
    protected_areas = n_distinct(
      data$wdpaid[treated_rows][!is.na(data$wdpaid[treated_rows])]
    ),
    treated_pixels = length(treated_rows),
    treated_cells_25km = count_cells_25km(data, treated_rows),
    indonesia_observations = length(indonesia_rows),
    indonesia_matched_pairs = n_distinct(
      data$matched_pair_id[indonesia_rows]
    ),
    indonesia_protected_areas = n_distinct(
      data$wdpaid[indonesia_treated_rows][
        !is.na(data$wdpaid[indonesia_treated_rows])
      ]
    ),
    indonesia_treated_pixels = length(indonesia_treated_rows),
    indonesia_cells_25km = count_cells_25km(
      data,
      indonesia_treated_rows
    ),
    indonesia_avg_abs_smd = indonesia_balance$avg_abs_smd,
    balance_threshold = balance_threshold,
    indonesia_passes_balance = indonesia_passes_balance,
    indonesia_in_stage = length(indonesia_treated_rows) > 0L,
    manifest_contract_version = manifest$contract_version,
    manifest_expected_strata = as.integer(manifest$expected_strata),
    manifest_input_md5 = input_fingerprint$value,
    stringsAsFactors = FALSE
  )
}

sample_summary <- rbind(
  summarize_stage(d, "pre_balance"),
  summarize_stage(d_post_balance, "post_balance")
)
if (sample_summary$indonesia_matched_pairs[1] != 3727L ||
    sample_summary$indonesia_observations[1] != 7454L ||
    any(sample_summary$observations != 2L * sample_summary$matched_pairs) ||
    any(sample_summary$treated_pixels != sample_summary$matched_pairs)) {
  stop("Matched-sample audit failed pair/count identities.")
}
rm(d)
gc()

cat(sprintf(
  paste0(
    "  Balance filter avg |SMD| <= %.2f: %d -> %d countries; ",
    "%s -> %s observations.\n"
  ),
  balance_threshold,
  sample_summary$countries[1],
  sample_summary$countries[2],
  format(sample_summary$observations[1], big.mark = ","),
  format(sample_summary$observations[2], big.mark = ",")
))
cat(sprintf(
  paste0(
    "  Indonesia avg |SMD| = %.6f; %s; post-balance pairs = %s.\n"
  ),
  indonesia_balance$avg_abs_smd,
  if (indonesia_passes_balance) "retained" else "excluded",
  format(sample_summary$indonesia_matched_pairs[2], big.mark = ",")
))


#################################
### stage figure and publish ####
#################################

balance_plot_data <- balance
balance_plot_data$country_name <- with(
  balance_plot_data,
  reorder(country_name, avg_abs_smd)
)
x_cap <- 1
balance_plot_data$avg_abs_smd_capped <- pmin(
  balance_plot_data$avg_abs_smd,
  x_cap
)
balance_plot_data$is_capped <- balance_plot_data$avg_abs_smd > x_cap
balance_plot <- ggplot(
  balance_plot_data,
  aes(x = avg_abs_smd_capped, y = country_name)
) +
  geom_vline(
    xintercept = balance_threshold,
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_vline(xintercept = 0, linetype = "solid", linewidth = 0.5) +
  geom_point(aes(shape = is_capped), colour = "#e66101", size = 2) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    guide = "none"
  ) +
  scale_x_continuous(limits = c(0, x_cap), breaks = seq(0, x_cap, 0.1)) +
  labs(x = "Average |Std. Mean Difference|", y = "") +
  theme(
    text = element_text(size = 10),
    axis.text.y = element_text(size = 6),
    panel.background = element_blank(),
    panel.grid.major.x = element_line(colour = "grey90"),
    panel.grid.major.y = element_line(colour = "grey95")
  )
figure_height <- max(20, nrow(balance_plot_data) * 0.4)

# No production artifact is replaced until every manifest, join, pair, delta,
# balance, and sample-summary gate above has passed and all four replacement
# artifacts have been written and read back successfully.
publish_outputs <- function() {
  output_paths <- c(
    matched = matched_output_path,
    balance = balance_output_path,
    summary = summary_output_path,
    figure = balance_figure_path
  )
  invisible(lapply(dirname(output_paths), dir.create,
                   recursive = TRUE, showWarnings = FALSE))
  staged_paths <- c(
    matched = tempfile(
      pattern = paste0(".", basename(matched_output_path), "."),
      tmpdir = dirname(matched_output_path),
      fileext = ".fst"
    ),
    balance = tempfile(
      pattern = paste0(".", basename(balance_output_path), "."),
      tmpdir = dirname(balance_output_path),
      fileext = ".rds"
    ),
    summary = tempfile(
      pattern = paste0(".", basename(summary_output_path), "."),
      tmpdir = dirname(summary_output_path),
      fileext = ".csv"
    ),
    figure = tempfile(
      pattern = paste0(
        ".",
        tools::file_path_sans_ext(basename(balance_figure_path)),
        "."
      ),
      tmpdir = dirname(balance_figure_path),
      fileext = ".jpg"
    )
  )
  backup_paths <- vapply(names(output_paths), function(name) {
    tempfile(
      pattern = paste0(".", basename(output_paths[name]), ".backup."),
      tmpdir = dirname(output_paths[name]),
      fileext = ".tmp"
    )
  }, character(1))
  backed_up <- character()
  promoted <- character()
  promotion_complete <- FALSE
  on.exit({
    if (!promotion_complete) {
      for (name in rev(promoted)) {
        if (file.exists(output_paths[name])) unlink(output_paths[name])
      }
      for (name in rev(backed_up)) {
        if (file.exists(backup_paths[name])) {
          if (!file.rename(backup_paths[name], output_paths[name])) {
            warning("Rollback failed for output: ", output_paths[name])
          }
        }
      }
    }
    unlink(staged_paths)
    if (promotion_complete) unlink(backup_paths)
  }, add = TRUE)

  write_fst(d_post_balance, staged_paths["matched"])
  matched_check <- tryCatch(
    fst::fst(staged_paths["matched"]),
    error = identity
  )
  if (inherits(matched_check, "error") ||
      nrow(matched_check) != nrow(d_post_balance) ||
      !identical(colnames(matched_check), names(d_post_balance))) {
    stop("Staged matched FST failed read-back validation.")
  }

  saveRDS(balance, staged_paths["balance"])
  balance_check <- tryCatch(
    readRDS(staged_paths["balance"]),
    error = identity
  )
  if (inherits(balance_check, "error") ||
      !identical(names(balance_check), names(balance)) ||
      nrow(balance_check) != nrow(balance)) {
    stop("Staged balance RDS failed read-back validation.")
  }

  write.csv(
    sample_summary,
    staged_paths["summary"],
    row.names = FALSE,
    na = ""
  )
  summary_check <- tryCatch(
    read.csv(staged_paths["summary"], stringsAsFactors = FALSE),
    error = identity
  )
  if (inherits(summary_check, "error") ||
      !identical(names(summary_check), names(sample_summary)) ||
      nrow(summary_check) != nrow(sample_summary)) {
    stop("Staged matched-sample CSV failed read-back validation.")
  }

  ggsave(
    filename = staged_paths["figure"],
    plot = balance_plot,
    units = "cm",
    width = 20,
    height = figure_height,
    dpi = 300
  )
  if (!file.exists(staged_paths["figure"]) ||
      file.info(staged_paths["figure"])$size <= 0L) {
    stop("Staged balance figure failed validation.")
  }

  # Promote ancillary audit artifacts first and the matched FST last so no
  # downstream reader can observe the new sample before its audits exist.
  promotion_order <- c("balance", "summary", "figure", "matched")
  for (name in promotion_order) {
    if (file.exists(output_paths[name])) {
      if (!file.rename(output_paths[name], backup_paths[name])) {
        stop("Could not stage existing output for rollback: ", output_paths[name])
      }
      backed_up <- c(backed_up, name)
    }
    atomic_rename(staged_paths[name], output_paths[name])
    promoted <- c(promoted, name)
  }
  promotion_complete <- TRUE
  unlink(backup_paths)
  invisible(output_paths)
}

publish_outputs()

cat("=== Matched-data assembly complete ===\n")
cat("  Matched FST: ", matched_output_path, "\n", sep = "")
cat("  Balance audit: ", balance_output_path, "\n", sep = "")
cat("  Sample audit: ", summary_output_path, "\n", sep = "")
