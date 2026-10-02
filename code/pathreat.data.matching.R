library(tidyr)
library(dplyr)
library(fst)
library(parallel)
library(pbapply)

## load config (provides mlist, paths, and other shared settings)
source("code/pathreat.analysis.config.R")


#################################
### matching contract ###########
#################################

matching_manifest_schema_version <- "1.0.0"
matching_contract_version <- "pathreat_matching_v3"
matching_manifest_path <- "data/store/pathreat.data.matching.manifest.Rds"
matching_output_dir <- "data/store/newmatchresults"
matching_expected_strata <- 668L
# Approved 2026-07-15 from the full 668-stratum effective-sample audit.
# This is embedded so production matching never depends on a data/tmp artifact.
matching_expected_drop_keys <- c(
  "100.6.0::cropsuit2",
  "100.6.0::forest2000_parea",
  "100.99.0::cropsuit2",
  "100.99.0::forest2000_parea",
  "102.7.0::pop_2000",
  "112.13.0::forest2000_parea",
  "112.8.0::forest2000_parea",
  "116.13.0::forest2000_parea",
  "123.13.0::forest2000_parea",
  "131.13.0::forest2000_parea",
  "131.9.0::forest2000_parea",
  "14.13.155::forest2000_parea",
  "14.13.156::forest2000_parea",
  "14.13.157::forest2000_parea",
  "14.13.158::forest2000_parea",
  "14.13.161::forest2000_parea",
  "145.13.0::forest2000_parea",
  "148.13.0::forest2000_parea",
  "155.13.0::forest2000_parea",
  "155.7.0::forest2000_parea",
  "155.9.0::forest2000_parea",
  "155.98.0::forest2000_parea",
  "156.10.0::forest2000_parea",
  "156.9.0::forest2000_parea",
  "16.8.0::forest2000_parea",
  "163.13.0::forest2000_parea",
  "163.8.0::forest2000_parea",
  "164.99.0::forest2000_parea",
  "166.13.0::forest2000_parea",
  "176.13.0::forest2000_parea",
  "180.11.2599::cropsuit2",
  "180.11.2617::cropsuit2",
  "180.11.2622::cropsuit2",
  "180.11.2632::cropsuit2",
  "180.11.2633::cropsuit2",
  "180.11.2647::cropsuit2",
  "180.11.2667::cropsuit2",
  "180.6.2627::cropsuit2",
  "192.13.0::forest2000_parea",
  "209.13.0::forest2000_parea",
  "209.7.0::forest2000_parea",
  "212.1.0::forest2000_parea",
  "213.11.0::cropsuit2",
  "225.13.0::forest2000_parea",
  "225.9.0::forest2000_parea",
  "232.13.0::forest2000_parea",
  "232.8.0::forest2000_parea",
  "245.13.0::forest2000_parea",
  "30.9.0::forest2000_parea",
  "37.7.0::forest2000_parea",
  "4.13.0::forest2000_parea",
  "4.9.0::forest2000_parea",
  "41.9.0::forest2000_parea",
  "41.9.0::pop_2000",
  "42.11.435::cropsuit2",
  "42.11.435::pop_2000",
  "42.11.438::cropsuit2",
  "42.11.441::cropsuit2",
  "42.6.438::cropsuit2",
  "42.6.438::pop_2000",
  "42.98.431::pop_2000",
  "46.13.0::forest2000_parea",
  "62.13.0::forest2000_parea",
  "66.13.0::forest2000_parea",
  "69.7.0::forest2000_parea",
  "7.13.0::forest2000_parea",
  "81.14.0::forest2000_parea"
)
matching_expected_drop_parts <- strsplit(
  matching_expected_drop_keys,
  "::",
  fixed = TRUE
)
matching_expected_drop_set <- data.frame(
  stratum_id = vapply(matching_expected_drop_parts, `[[`, character(1), 1L),
  covariate = vapply(matching_expected_drop_parts, `[[`, character(1), 2L),
  stringsAsFactors = FALSE
)
matching_affected_strata <- unique(matching_expected_drop_set$stratum_id)
rm(matching_expected_drop_parts)

matching_result_columns <- c(
  "country_rast", "biome_raster", "row_id", "pair_id",
  "weight", "pixel_id", "treat", "nid"
)

matching_filename <- function(country_rast, biome_raster, bigregion) {
  paste0(
    "matchresult.cty-", country_rast,
    ".biome-", biome_raster,
    ".bigregion-", bigregion,
    ".maha.1-1.caliper-100.Rds"
  )
}

matching_stratum_id <- function(country_rast, biome_raster, bigregion) {
  paste(country_rast, biome_raster, bigregion, sep = ".")
}

matching_timestamp <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
}


#################################
### atomic writes and hashes ####
#################################

atomic_save_rds <- function(object, path, compress = TRUE) {
  output_dir <- dirname(path)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  temporary_path <- tempfile(
    pattern = paste0(".", basename(path), "."),
    tmpdir = output_dir,
    fileext = ".tmp"
  )
  on.exit(unlink(temporary_path), add = TRUE)

  saveRDS(object, temporary_path, compress = compress)
  if (!file.rename(temporary_path, path)) {
    stop(sprintf("Atomic rename failed: %s -> %s", temporary_path, path))
  }
  invisible(path)
}

matching_file_md5 <- function(path) {
  unname(tools::md5sum(path))
}

matching_md5_is_valid <- function(value) {
  length(value) == 1L &&
    is.character(value) &&
    !is.na(value) &&
    grepl("^[0-9a-f]{32}$", value)
}

matching_input_fingerprint_is_valid <- function(fingerprint) {
  required <- c(
    "algorithm", "value", "size_bytes", "mtime_unix", "n_rows", "columns"
  )
  if (!is.list(fingerprint) ||
      any(!required %in% names(fingerprint))) {
    return(FALSE)
  }

  isTRUE(
    identical(fingerprint$algorithm, "md5") &&
      matching_md5_is_valid(fingerprint$value) &&
      length(fingerprint$size_bytes) == 1L &&
      is.numeric(fingerprint$size_bytes) &&
      is.finite(fingerprint$size_bytes) &&
      fingerprint$size_bytes > 0 &&
      length(fingerprint$mtime_unix) == 1L &&
      is.numeric(fingerprint$mtime_unix) &&
      is.finite(fingerprint$mtime_unix) &&
      length(fingerprint$n_rows) == 1L &&
      is.numeric(fingerprint$n_rows) &&
      is.finite(fingerprint$n_rows) &&
      fingerprint$n_rows >= 0 &&
      fingerprint$n_rows == floor(fingerprint$n_rows) &&
      is.character(fingerprint$columns) &&
      length(fingerprint$columns) > 0L &&
      !anyNA(fingerprint$columns)
  )
}

matching_input_fingerprint <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("Unmatched input is missing: %s", path))
  }

  file_metadata <- file.info(path)
  fst_metadata <- fst::fst(path)
  md5 <- matching_file_md5(path)
  fingerprint <- list(
    algorithm = "md5",
    value = md5,
    size_bytes = as.numeric(file_metadata$size),
    mtime_unix = as.numeric(file_metadata$mtime),
    n_rows = as.numeric(nrow(fst_metadata)),
    columns = colnames(fst_metadata)
  )
  if (!matching_input_fingerprint_is_valid(fingerprint)) {
    stop(sprintf("Could not compute a valid input fingerprint for %s", path))
  }

  fingerprint
}

matching_input_fingerprint_differences <- function(expected, observed) {
  fields <- c(
    "algorithm", "value", "size_bytes", "mtime_unix", "n_rows", "columns"
  )
  fields[!vapply(fields, function(field) {
    identical(expected[[field]], observed[[field]])
  }, logical(1))]
}

recheck_matching_input <- function(
    path,
    expected_fingerprint,
    fingerprint_fun = matching_input_fingerprint) {
  if (!matching_input_fingerprint_is_valid(expected_fingerprint)) {
    return(list(
      valid = FALSE,
      errors = "Original unmatched-input fingerprint is malformed",
      fingerprint = NULL
    ))
  }

  observed <- tryCatch(fingerprint_fun(path), error = identity)
  if (inherits(observed, "error")) {
    return(list(
      valid = FALSE,
      errors = paste0(
        "Could not re-fingerprint unmatched input after matching: ",
        conditionMessage(observed)
      ),
      fingerprint = NULL
    ))
  }
  if (!matching_input_fingerprint_is_valid(observed)) {
    return(list(
      valid = FALSE,
      errors = "Final unmatched-input fingerprint is malformed",
      fingerprint = observed
    ))
  }

  differences <- matching_input_fingerprint_differences(
    expected_fingerprint,
    observed
  )
  errors <- if (length(differences) == 0L) {
    character()
  } else {
    labels <- differences
    labels[labels == "value"] <- "MD5"
    sprintf(
      "Unmatched input fingerprint changed during matching: %s",
      paste(labels, collapse = ", ")
    )
  }

  list(
    valid = length(errors) == 0L,
    errors = errors,
    fingerprint = observed
  )
}


#################################
### covariate selection #########
#################################

covariate_diagnostics <- function(x, covariates = colnames(x)) {
  x <- as.data.frame(x)
  missing_covariates <- setdiff(covariates, names(x))
  if (length(missing_covariates) > 0) {
    stop(sprintf(
      "Covariates are missing from matching data: %s",
      paste(missing_covariates, collapse = ", ")
    ))
  }

  diagnostics <- lapply(covariates, function(covariate) {
    values <- x[[covariate]]
    if (!is.numeric(values)) {
      stop(sprintf("Matching covariate is not numeric: %s", covariate))
    }
    if (anyNA(values) || any(!is.finite(values))) {
      stop(sprintf("Matching covariate contains missing/non-finite values: %s", covariate))
    }

    n_unique <- length(unique(values))
    value_range <- range(values)
    value_mean <- mean(values)
    value_sd <- if (length(values) > 1L) sd(values) else NA_real_
    relative_sd <- value_sd / max(1, abs(value_mean))
    drop <- n_unique < 2L || is.na(relative_sd) || relative_sd <= 1e-6

    data.frame(
      covariate = covariate,
      n_unique = as.integer(n_unique),
      min = value_range[1],
      max = value_range[2],
      mean = value_mean,
      sd = value_sd,
      relative_sd = relative_sd,
      drop = drop,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, diagnostics)
}

select_matching_covariates <- function(x, covariates = colnames(x)) {
  diagnostics <- covariate_diagnostics(x, covariates)
  selected <- diagnostics$covariate[!diagnostics$drop]
  dropped <- diagnostics$covariate[diagnostics$drop]

  if (length(selected) == 0L) {
    stop("No usable matching covariates remain after near-constant screening")
  }

  list(
    selected = selected,
    dropped = dropped,
    diagnostics = diagnostics
  )
}

matching_drop_table <- function(strata) {
  pieces <- lapply(seq_len(nrow(strata)), function(i) {
    dropped <- strata$dropped_covariates[[i]]
    if (length(dropped) == 0L) return(NULL)
    data.frame(
      stratum_id = rep(strata$stratum_id[i], length(dropped)),
      covariate = dropped,
      stringsAsFactors = FALSE
    )
  })
  pieces <- Filter(Negate(is.null), pieces)
  if (length(pieces) == 0L) {
    return(data.frame(
      stratum_id = character(),
      covariate = character(),
      stringsAsFactors = FALSE
    ))
  }
  do.call(rbind, pieces)
}

matching_drop_set_is_expected <- function(actual, expected = matching_expected_drop_set) {
  actual_key <- sort(paste(actual$stratum_id, actual$covariate, sep = "::"))
  expected_key <- sort(paste(expected$stratum_id, expected$covariate, sep = "::"))
  identical(actual_key, expected_key)
}


#################################
### input and expected strata ###
#################################

load_matching_input <- function(path = paths$data_unmatched) {
  match_cols <- c(
    "pixel_id", "country_rast", "country", "biome_raster",
    "treat", "buffer_5km", "pa_year_designated", "admin1_id", mlist
  )
  d <- read_fst(path, columns = match_cols)

  n_before <- nrow(d)
  complete_matching_data <- complete.cases(d[, c("biome_raster", "country_rast", mlist)])
  d <- d[complete_matching_data, ]
  cat(sprintf(
    "Complete matching fields: dropped %s rows\n",
    format(n_before - nrow(d), big.mark = ",")
  ))

  drop_after_2020 <- which(d$pa_year_designated > 2020)
  if (length(drop_after_2020) > 0L) d <- d[-drop_after_2020, ]

  drop_before_2001 <- which(d$pa_year_designated <= 2000)
  if (length(drop_before_2001) > 0L) d <- d[-drop_before_2001, ]

  drop_buffer_controls <- which(d$treat == 0 & d$buffer_5km == 1)
  if (length(drop_buffer_controls) > 0L) d <- d[-drop_buffer_controls, ]

  drop_missing_admin1 <- which(is.na(d$admin1_id))
  cat(sprintf(
    "admin1_id: dropping %s pixels (%.1f%%) with NA\n",
    format(length(drop_missing_admin1), big.mark = ","),
    100 * length(drop_missing_admin1) / nrow(d)
  ))
  if (length(drop_missing_admin1) > 0L) d <- d[-drop_missing_admin1, ]

  if (anyNA(d$treat) || !all(d$treat %in% c(0, 1))) {
    stop("treat must be observed and coded 0/1 after matching-input filters")
  }
  if (anyDuplicated(d$pixel_id)) {
    stop("pixel_id must be unique in the filtered unmatched input")
  }

  big_countries <- c("Brazil", "Australia", "Canada", "Russia")
  is_big <- d$country %in% big_countries
  d$bigregion <- 0L
  d$bigregion[is_big] <- d$admin1_id[is_big]
  cat(sprintf(
    "bigregion: %d countries subdivided by admin1 (%d unique regions)\n",
    sum(big_countries %in% unique(d$country)),
    length(unique(d$bigregion[d$bigregion > 0]))
  ))

  d
}

build_expected_strata <- function(d) {
  strata_counts <- d %>%
    group_by(country_rast, country, biome_raster, bigregion, treat) %>%
    summarize(obs = n(), .groups = "drop") %>%
    pivot_wider(
      names_from = treat,
      values_from = obs,
      values_fill = 0,
      names_glue = "n_{treat}"
    ) %>%
    rename(
      st.no.treated = n_1,
      st.no.control = n_0
    ) %>%
    mutate(st.no.total = st.no.treated + st.no.control) %>%
    data.frame()

  country_counts <- d %>%
    group_by(country_rast, treat) %>%
    summarize(obs = n(), .groups = "drop") %>%
    pivot_wider(
      names_from = treat,
      values_from = obs,
      values_fill = 0,
      names_glue = "n_{treat}"
    ) %>%
    rename(
      c.no.treated = n_1,
      c.no.control = n_0
    ) %>%
    mutate(c.no.total = c.no.treated + c.no.control) %>%
    data.frame()

  strata <- merge(strata_counts, country_counts, by = "country_rast", all.x = TRUE)
  strata <- strata[strata$c.no.control >= 20 & strata$st.no.treated > 0, ]
  strata$too.few.controls <- as.integer(strata$st.no.control < 20)
  strata$stratum_id <- matching_stratum_id(
    strata$country_rast,
    strata$biome_raster,
    strata$bigregion
  )
  strata <- strata[order(strata$st.no.treated, strata$stratum_id), ]
  rownames(strata) <- NULL

  if (anyDuplicated(strata$stratum_id)) {
    stop("Expected matching strata contain duplicate identifiers")
  }
  if (nrow(strata) != matching_expected_strata) {
    stop(sprintf(
      "Expected %d matching strata under contract %s, found %d",
      matching_expected_strata,
      matching_contract_version,
      nrow(strata)
    ))
  }

  strata
}

build_matching_state <- function(d) {
  state <- list(
    d_covars = as.matrix(d[, mlist]),
    d_pixel_id = d$pixel_id,
    d_treat = d$treat,
    d_biome = d$biome_raster,
    d_rows_by_stratum = split(
      seq_len(nrow(d)),
      matching_stratum_id(d$country_rast, d$biome_raster, d$bigregion)
    ),
    d_rows_by_country = split(seq_len(nrow(d)), d$country_rast)
  )
  state
}

effective_stratum_rows <- function(g, state) {
  rows <- state$d_rows_by_stratum[[g$stratum_id]]
  if (is.null(rows) || length(rows) == 0L) {
    stop(sprintf("No source rows found for stratum %s", g$stratum_id))
  }

  pooling_scope <- "none"
  if (g$too.few.controls == 1L) {
    rows_country <- state$d_rows_by_country[[as.character(g$country_rast)]]
    if (g$bigregion > 0L) {
      rows_pool <- rows_country[state$d_biome[rows_country] == g$biome_raster]
      pooling_scope <- "country_biome"
    } else {
      rows_pool <- rows_country
      pooling_scope <- "country"
    }
    rows_treated <- rows[state$d_treat[rows] == 1]
    rows_control <- rows_pool[state$d_treat[rows_pool] == 0]
    rows <- c(rows_treated, rows_control)
  }

  list(
    rows = rows,
    pooling_scope = pooling_scope,
    source_treated = sum(state$d_treat[rows] == 1),
    source_controls = sum(state$d_treat[rows] == 0)
  )
}

prepare_strata_contract <- function(strata, state) {
  selected_covariates <- vector("list", nrow(strata))
  dropped_covariates <- vector("list", nrow(strata))
  dropped_diagnostics <- vector("list", nrow(strata))
  source_controls_effective <- integer(nrow(strata))
  source_rows_effective <- integer(nrow(strata))
  pooling_scope <- character(nrow(strata))

  for (i in seq_len(nrow(strata))) {
    g <- strata[i, ]
    effective <- effective_stratum_rows(g, state)
    selection <- select_matching_covariates(
      state$d_covars[effective$rows, , drop = FALSE],
      mlist
    )

    selected_covariates[[i]] <- selection$selected
    dropped_covariates[[i]] <- selection$dropped
    dropped_diagnostics[[i]] <- selection$diagnostics[selection$diagnostics$drop, ]
    source_controls_effective[i] <- effective$source_controls
    source_rows_effective[i] <- length(effective$rows)
    pooling_scope[i] <- effective$pooling_scope

    if (length(selection$dropped) > 0L) {
      for (j in which(selection$diagnostics$drop)) {
        diagnostic <- selection$diagnostics[j, ]
        cat(sprintf(
          paste0(
            "Dropped covariate | stratum=%s covariate=%s unique=%d ",
            "range=[%.17g, %.17g] sd=%.17g relative_sd=%.17g\n"
          ),
          g$stratum_id,
          diagnostic$covariate,
          diagnostic$n_unique,
          diagnostic$min,
          diagnostic$max,
          diagnostic$sd,
          diagnostic$relative_sd
        ))
      }
    }
  }

  manifest_strata <- data.frame(
    stratum_id = strata$stratum_id,
    country_rast = as.integer(strata$country_rast),
    country = as.character(strata$country),
    biome_raster = as.integer(strata$biome_raster),
    bigregion = as.integer(strata$bigregion),
    source_treated = as.integer(strata$st.no.treated),
    source_controls_stratum = as.integer(strata$st.no.control),
    source_controls_effective = source_controls_effective,
    source_rows_effective = source_rows_effective,
    pooling_scope = pooling_scope,
    output_file = file.path(
      matching_output_dir,
      matching_filename(
        strata$country_rast,
        strata$biome_raster,
        strata$bigregion
      )
    ),
    pair_count = rep(NA_integer_, nrow(strata)),
    status = rep("failure", nrow(strata)),
    error_text = rep(NA_character_, nrow(strata)),
    action = rep("failed", nrow(strata)),
    validation_status = rep("failed", nrow(strata)),
    reuse_reason = rep("not_evaluated", nrow(strata)),
    output_md5 = rep(NA_character_, nrow(strata)),
    output_size_bytes = rep(NA_real_, nrow(strata)),
    stringsAsFactors = FALSE
  )
  manifest_strata$selected_covariates <- I(selected_covariates)
  manifest_strata$dropped_covariates <- I(dropped_covariates)
  manifest_strata$dropped_diagnostics <- I(dropped_diagnostics)
  manifest_strata$validation_errors <- I(
    replicate(nrow(strata), character(), simplify = FALSE)
  )

  list(source_strata = strata, manifest_strata = manifest_strata)
}


#################################
### result validation ###########
#################################

validate_match_result <- function(result, g, state, effective = NULL) {
  errors <- character()
  if (!is.data.frame(result)) {
    return(list(valid = FALSE, errors = "not_data_frame", pair_count = NA_integer_))
  }
  if (!identical(names(result), matching_result_columns)) {
    return(list(valid = FALSE, errors = "wrong_schema", pair_count = NA_integer_))
  }
  if (is.null(effective)) effective <- effective_stratum_rows(g, state)

  if (nrow(result) == 0L) errors <- c(errors, "empty_result")
  if (nrow(result) %% 2L != 0L) errors <- c(errors, "odd_row_count")
  if (anyNA(result$country_rast) ||
      !all(result$country_rast == g$country_rast)) {
    errors <- c(errors, "country_identifier_mismatch")
  }
  if (anyNA(result$biome_raster) ||
      !all(result$biome_raster == g$biome_raster)) {
    errors <- c(errors, "biome_identifier_mismatch")
  }
  if (anyNA(result$treat) || !all(result$treat %in% c(0, 1))) {
    errors <- c(errors, "invalid_treatment")
  }

  n_treated <- sum(result$treat == 1, na.rm = TRUE)
  n_control <- sum(result$treat == 0, na.rm = TRUE)
  expected_pairs <- as.integer(g$st.no.treated)
  if (n_treated != n_control || n_treated * 2L != nrow(result)) {
    errors <- c(errors, "treatment_control_parity")
  }
  if (n_treated != expected_pairs) {
    errors <- c(errors, "pair_count_mismatch")
  }
  if (anyNA(result$weight) ||
      any(!is.finite(result$weight)) ||
      any(result$weight <= 0)) {
    errors <- c(errors, "invalid_weights")
  }
  if (anyNA(result$pixel_id) ||
      any(!is.finite(result$pixel_id)) ||
      any(result$pixel_id <= 0)) {
    errors <- c(errors, "invalid_pixel_id")
  }
  if (anyDuplicated(result$pixel_id[result$treat == 1])) {
    errors <- c(errors, "treated_pixel_reused")
  }
  if (!identical(result$nid, seq_len(nrow(result)))) {
    errors <- c(errors, "invalid_nid")
  }

  source_rows <- effective$rows
  source_pixel_id <- state$d_pixel_id[source_rows]
  source_treat <- state$d_treat[source_rows]
  source_position <- match(result$pixel_id, source_pixel_id)
  if (anyNA(source_position)) {
    errors <- c(errors, "pixel_not_in_effective_stratum")
  } else if (!all(source_treat[source_position] == result$treat)) {
    errors <- c(errors, "pixel_treatment_mismatch")
  }

  valid_index <- function(x) {
    !anyNA(x) && all(is.finite(x)) && all(x == floor(x)) &&
      all(x >= 1) && all(x <= length(source_rows))
  }
  row_id_valid <- valid_index(result$row_id)
  pair_id_valid <- valid_index(result$pair_id)
  if (!row_id_valid) errors <- c(errors, "invalid_row_id")
  if (!pair_id_valid) errors <- c(errors, "invalid_pair_id")
  if (row_id_valid &&
      !all(source_pixel_id[as.integer(result$row_id)] == result$pixel_id)) {
    errors <- c(errors, "row_id_pixel_mismatch")
  }
  if (row_id_valid &&
      !all(source_treat[as.integer(result$row_id)] == result$treat)) {
    errors <- c(errors, "row_id_treatment_mismatch")
  }
  if (pair_id_valid &&
      !all(source_treat[as.integer(result$pair_id)] == 1 - result$treat)) {
    errors <- c(errors, "pair_id_treatment_mismatch")
  }

  expected_treated <- source_pixel_id[source_treat == 1]
  observed_treated <- result$pixel_id[result$treat == 1]
  if (length(observed_treated) != length(expected_treated) ||
      !setequal(observed_treated, expected_treated)) {
    errors <- c(errors, "treated_source_coverage")
  }

  if (n_treated == n_control && n_treated > 0L) {
    treated_rows <- which(result$treat == 1)
    control_rows <- which(result$treat == 0)
    ordered_layout <- length(treated_rows) == n_treated &&
      identical(treated_rows, seq_len(n_treated)) &&
      identical(control_rows, n_treated + seq_len(n_control))
    if (!ordered_layout) {
      errors <- c(errors, "unexpected_pair_layout")
    } else if (!all(result$row_id[treated_rows] == result$pair_id[control_rows]) ||
               !all(result$pair_id[treated_rows] == result$row_id[control_rows]) ||
               !all(result$weight[treated_rows] == result$weight[control_rows])) {
      errors <- c(errors, "pair_cross_reference_mismatch")
    }
  }

  list(
    valid = length(errors) == 0L,
    errors = unique(errors),
    pair_count = as.integer(n_treated)
  )
}

validate_match_file <- function(path, g, state, effective = NULL) {
  if (!file.exists(path)) {
    return(list(
      valid = FALSE,
      errors = "missing_file",
      pair_count = NA_integer_,
      md5 = NA_character_,
      size_bytes = NA_real_
    ))
  }

  result <- tryCatch(readRDS(path), error = identity)
  if (inherits(result, "error")) {
    return(list(
      valid = FALSE,
      errors = paste0("unreadable_file: ", conditionMessage(result)),
      pair_count = NA_integer_,
      md5 = NA_character_,
      size_bytes = as.numeric(file.info(path)$size)
    ))
  }

  validation <- validate_match_result(result, g, state, effective)
  list(
    valid = validation$valid,
    errors = validation$errors,
    pair_count = validation$pair_count,
    md5 = if (validation$valid) matching_file_md5(path) else NA_character_,
    size_bytes = as.numeric(file.info(path)$size)
  )
}


#################################
### structured matching #########
#################################

run_matching_engine <- function(treatment, covariates,
                                match_fun = Matching::Match) {
  tryCatch(
    {
      result <- match_fun(
        Tr = treatment,
        X = covariates,
        replace = TRUE,
        M = 1,
        Weight = 2,
        ties = FALSE
      )
      if (is.null(result)) stop("Matching engine returned NULL")
      list(success = TRUE, result = result, error_text = NA_character_)
    },
    error = function(e) {
      list(
        success = FALSE,
        result = NULL,
        error_text = conditionMessage(e)
      )
    }
  )
}

match_one_stratum <- function(i, source_strata, manifest_strata, state,
                              match_fun = Matching::Match) {
  t0 <- Sys.time()
  g <- source_strata[i, ]
  selected <- manifest_strata$selected_covariates[[i]]
  output_file <- manifest_strata$output_file[i]

  tryCatch(
    {
      effective <- effective_stratum_rows(g, state)
      if (length(selected) == 0L) {
        stop("No usable matching covariates remain")
      }

      rows <- effective$rows
      covariates <- state$d_covars[rows, selected, drop = FALSE]
      treatment <- state$d_treat[rows]
      pixel_id <- state$d_pixel_id[rows]

      engine <- run_matching_engine(treatment, covariates, match_fun)
      if (!engine$success) stop(engine$error_text)

      match_result <- engine$result
      index_treated <- match_result$index.treated
      index_control <- match_result$index.control
      pair_count <- length(index_treated)
      if (length(index_control) != pair_count ||
          length(match_result$weights) != pair_count) {
        stop("Matching engine returned inconsistent pair vectors")
      }

      result <- data.frame(
        country_rast = rep.int(as.integer(g$country_rast), 2L * pair_count),
        biome_raster = rep.int(as.integer(g$biome_raster), 2L * pair_count),
        row_id = c(index_treated, index_control),
        pair_id = c(index_control, index_treated),
        weight = c(match_result$weights, match_result$weights),
        pixel_id = c(pixel_id[index_treated], pixel_id[index_control]),
        treat = c(treatment[index_treated], treatment[index_control]),
        nid = seq_len(2L * pair_count)
      )

      validation <- validate_match_result(result, g, state, effective)
      if (!validation$valid) {
        stop(sprintf(
          "New match result failed validation: %s",
          paste(validation$errors, collapse = "; ")
        ))
      }

      atomic_save_rds(result, output_file, compress = FALSE)
      saved_validation <- validate_match_file(output_file, g, state, effective)
      if (!saved_validation$valid) {
        stop(sprintf(
          "Saved match result failed validation: %s",
          paste(saved_validation$errors, collapse = "; ")
        ))
      }

      cat(sprintf(
        "Matched %s | T=%d C=%d pairs=%d | %.1f sec\n",
        g$stratum_id,
        effective$source_treated,
        effective$source_controls,
        saved_validation$pair_count,
        as.numeric(difftime(Sys.time(), t0, units = "secs"))
      ))

      list(
        index = i,
        success = TRUE,
        pair_count = saved_validation$pair_count,
        error_text = NA_character_,
        output_md5 = saved_validation$md5,
        output_size_bytes = saved_validation$size_bytes
      )
    },
    error = function(e) {
      error_text <- conditionMessage(e)
      cat(sprintf("Matching failed %s | %s\n", g$stratum_id, error_text))
      list(
        index = i,
        success = FALSE,
        pair_count = NA_integer_,
        error_text = error_text,
        output_md5 = NA_character_,
        output_size_bytes = NA_real_
      )
    }
  )
}

matching_worker_failure <- function(index, error_text) {
  error_text <- paste(as.character(error_text), collapse = " | ")
  if (length(error_text) != 1L ||
      is.na(error_text) ||
      !nzchar(error_text)) {
    error_text <- "Unspecified matching-worker failure"
  }

  list(
    index = as.integer(index),
    success = FALSE,
    pair_count = NA_integer_,
    error_text = error_text,
    output_md5 = NA_character_,
    output_size_bytes = NA_real_
  )
}

matching_worker_error_text <- function(result) {
  if (inherits(result, "condition")) {
    return(conditionMessage(result))
  }
  text <- paste(as.character(result), collapse = " | ")
  if (length(text) != 1L || is.na(text) || !nzchar(text)) {
    "unknown worker error"
  } else {
    text
  }
}

normalize_matching_worker_results <- function(results, rerun) {
  rerun <- as.integer(rerun)
  global_errors <- character()

  if (inherits(results, "condition") || inherits(results, "try-error")) {
    global_errors <- sprintf(
      "Parallel matching execution failed: %s",
      matching_worker_error_text(results)
    )
    results <- list()
  } else if (!is.list(results)) {
    global_errors <- sprintf(
      "Parallel matching returned a non-list result of class %s",
      paste(class(results), collapse = "/")
    )
    results <- list()
  }

  if (length(results) != length(rerun)) {
    global_errors <- c(
      global_errors,
      sprintf(
        "Parallel matching returned %d results for %d requested strata",
        length(results),
        length(rerun)
      )
    )
  }

  normalized <- lapply(seq_along(rerun), function(position) {
    expected_index <- rerun[position]
    if (position > length(results)) {
      return(matching_worker_failure(
        expected_index,
        sprintf(
          "Missing worker result for stratum %s",
          expected_index
        )
      ))
    }

    result <- results[[position]]
    if (inherits(result, "condition") || inherits(result, "try-error")) {
      return(matching_worker_failure(
        expected_index,
        sprintf(
          "Worker process failed for stratum %s: %s",
          expected_index,
          matching_worker_error_text(result)
        )
      ))
    }

    required <- c(
      "index", "success", "pair_count", "error_text",
      "output_md5", "output_size_bytes"
    )
    if (!is.list(result) || any(!required %in% names(result))) {
      return(matching_worker_failure(
        expected_index,
        sprintf(
          "Malformed worker result for stratum %s: wrong schema",
          expected_index
        )
      ))
    }

    valid_index <- length(result$index) == 1L &&
      is.numeric(result$index) &&
      is.finite(result$index) &&
      result$index == floor(result$index) &&
      as.integer(result$index) == expected_index
    valid_success <- length(result$success) == 1L &&
      is.logical(result$success) &&
      !is.na(result$success)
    if (!valid_index || !valid_success) {
      return(matching_worker_failure(
        expected_index,
        sprintf(
          "Malformed worker result for stratum %s: invalid index or status",
          expected_index
        )
      ))
    }

    if (!isTRUE(result$success)) {
      valid_error <- length(result$error_text) == 1L &&
        is.character(result$error_text) &&
        !is.na(result$error_text) &&
        nzchar(result$error_text)
      return(matching_worker_failure(
        expected_index,
        if (valid_error) {
          result$error_text
        } else {
          sprintf(
            "Malformed failed-worker result for stratum %s: missing error text",
            expected_index
          )
        }
      ))
    }

    valid_pair_count <- length(result$pair_count) == 1L &&
      is.numeric(result$pair_count) &&
      is.finite(result$pair_count) &&
      result$pair_count > 0 &&
      result$pair_count == floor(result$pair_count)
    valid_error_text <- length(result$error_text) == 1L &&
      is.character(result$error_text) &&
      is.na(result$error_text)
    valid_size <- length(result$output_size_bytes) == 1L &&
      is.numeric(result$output_size_bytes) &&
      is.finite(result$output_size_bytes) &&
      result$output_size_bytes > 0
    if (!valid_pair_count ||
        !valid_error_text ||
        !matching_md5_is_valid(result$output_md5) ||
        !valid_size) {
      return(matching_worker_failure(
        expected_index,
        sprintf(
          paste0(
            "Malformed successful-worker result for stratum %s: ",
            "invalid pair count, error state, or output fingerprint"
          ),
          expected_index
        )
      ))
    }

    list(
      index = expected_index,
      success = TRUE,
      pair_count = as.integer(result$pair_count),
      error_text = NA_character_,
      output_md5 = result$output_md5,
      output_size_bytes = as.numeric(result$output_size_bytes)
    )
  })

  list(
    results = normalized,
    errors = unique(global_errors)
  )
}


#################################
### manifest ####################
#################################

matching_specification <- function() {
  list(
    covariates = mlist,
    near_constant_rule = list(
      min_unique = 2L,
      relative_sd_threshold = 1e-6,
      formula = "sd(x) / max(1, abs(mean(x)))"
    ),
    match = list(
      "function" = "Matching::Match",
      M = 1L,
      Weight = 2L,
      replace = TRUE,
      ties = FALSE,
      caliper = NULL,
      filename_tag = "maha.1-1.caliper-100"
    ),
    filters = list(
      complete = c("biome_raster", "country_rast", mlist),
      pa_year_designated = "NA or 2001-2020",
      controls = "exclude buffer_5km == 1",
      admin1_id = "non-missing",
      country_controls_minimum = 20L
    ),
    strata = c("country_rast", "biome_raster", "bigregion"),
    big_countries = c("Brazil", "Australia", "Canada", "Russia"),
    pooling = list(
      trigger = "fewer than 20 controls in stratum",
      big_country = "controls from same country-biome",
      other_country = "controls from same country"
    )
  )
}

matching_manifest_complete <- function(strata) {
  required <- c(
    "stratum_id", "source_treated", "output_file", "pair_count",
    "status", "action", "validation_status", "reuse_reason",
    "output_md5", "output_size_bytes", "dropped_covariates"
  )
  if (!is.data.frame(strata) ||
      any(!required %in% names(strata))) {
    return(FALSE)
  }

  md5_valid <- vapply(
    as.character(strata$output_md5),
    matching_md5_is_valid,
    logical(1)
  )
  size_valid <- is.numeric(strata$output_size_bytes) &&
    length(strata$output_size_bytes) == nrow(strata) &&
    all(is.finite(strata$output_size_bytes)) &&
    all(strata$output_size_bytes > 0)
  pair_count_valid <- is.numeric(strata$pair_count) &&
    length(strata$pair_count) == nrow(strata) &&
    all(is.finite(strata$pair_count)) &&
    all(strata$pair_count == floor(strata$pair_count)) &&
    all(strata$pair_count == strata$source_treated)
  output_file_valid <- is.character(strata$output_file) &&
    !anyNA(strata$output_file) &&
    all(nzchar(strata$output_file)) &&
    !anyDuplicated(strata$output_file)

  isTRUE(
    nrow(strata) == matching_expected_strata &&
      !anyDuplicated(strata$stratum_id) &&
      output_file_valid &&
      all(strata$status == "success") &&
      all(strata$action %in% c("reused", "rematched")) &&
      all(strata$validation_status %in% c(
        "valid_reused",
        "new_output_valid"
      )) &&
      all(!is.na(strata$reuse_reason) & nzchar(strata$reuse_reason)) &&
      pair_count_valid &&
      all(md5_valid) &&
      size_valid &&
      matching_drop_set_is_expected(matching_drop_table(strata))
  )
}

build_matching_manifest <- function(input_fingerprint, strata, errors = character()) {
  complete <- matching_manifest_complete(strata) && length(errors) == 0L
  list(
    schema_version = matching_manifest_schema_version,
    contract_version = matching_contract_version,
    created_at_utc = matching_timestamp(),
    input = list(
      path = paths$data_unmatched,
      fingerprint = input_fingerprint
    ),
    specification = matching_specification(),
    expected_strata = matching_expected_strata,
    expected_drop_set = matching_expected_drop_set,
    strata = strata,
    complete = complete,
    errors = as.character(errors)
  )
}

matching_prior_manifest_failure <- function(reason) {
  list(valid = FALSE, reason = reason, manifest = NULL)
}

matching_covariate_contract_identical <- function(prior_strata,
                                                   current_strata) {
  required <- c(
    "stratum_id", "selected_covariates", "dropped_covariates"
  )
  if (!is.data.frame(prior_strata) ||
      !is.data.frame(current_strata) ||
      any(!required %in% names(prior_strata)) ||
      any(!required %in% names(current_strata)) ||
      anyDuplicated(prior_strata$stratum_id) ||
      anyDuplicated(current_strata$stratum_id)) {
    return(FALSE)
  }

  prior_index <- match(current_strata$stratum_id, prior_strata$stratum_id)
  if (anyNA(prior_index) || nrow(prior_strata) != nrow(current_strata)) {
    return(FALSE)
  }

  all(vapply(seq_len(nrow(current_strata)), function(i) {
    j <- prior_index[i]
    identical(
      prior_strata$selected_covariates[[j]],
      current_strata$selected_covariates[[i]]
    ) && identical(
      prior_strata$dropped_covariates[[j]],
      current_strata$dropped_covariates[[i]]
    )
  }, logical(1)))
}

matching_prior_manifest_proof <- function(path, input_fingerprint,
                                           current_strata) {
  if (!file.exists(path)) {
    return(matching_prior_manifest_failure("prior_manifest_missing"))
  }

  prior <- tryCatch(readRDS(path), error = identity)
  if (inherits(prior, "error")) {
    return(matching_prior_manifest_failure("prior_manifest_unreadable"))
  }

  required_top <- c(
    "schema_version", "contract_version", "input", "specification",
    "expected_strata", "expected_drop_set", "strata", "complete", "errors"
  )
  if (!is.list(prior) || any(!required_top %in% names(prior))) {
    return(matching_prior_manifest_failure("prior_manifest_schema_invalid"))
  }
  if (!identical(prior$schema_version, matching_manifest_schema_version) ||
      !identical(prior$contract_version, matching_contract_version)) {
    return(matching_prior_manifest_failure("prior_manifest_contract_mismatch"))
  }
  if (!isTRUE(prior$complete) || length(prior$errors) > 0L) {
    return(matching_prior_manifest_failure("prior_manifest_incomplete"))
  }
  if (!identical(prior$expected_strata, matching_expected_strata) ||
      !matching_drop_set_is_expected(prior$expected_drop_set)) {
    return(matching_prior_manifest_failure("prior_manifest_drop_contract_mismatch"))
  }
  if (!is.list(prior$input) ||
      !identical(prior$input$path, paths$data_unmatched) ||
      !identical(prior$input$fingerprint, input_fingerprint)) {
    return(matching_prior_manifest_failure("prior_manifest_input_mismatch"))
  }
  if (!identical(prior$specification, matching_specification())) {
    return(matching_prior_manifest_failure("prior_manifest_specification_mismatch"))
  }

  required_strata <- c(
    "stratum_id", "country_rast", "country", "biome_raster", "bigregion",
    "source_treated", "source_controls_stratum",
    "source_controls_effective", "source_rows_effective", "pooling_scope",
    "output_file", "pair_count", "status", "action", "validation_status",
    "reuse_reason", "output_md5", "output_size_bytes",
    "selected_covariates", "dropped_covariates"
  )
  if (!is.data.frame(prior$strata) ||
      any(!required_strata %in% names(prior$strata)) ||
      nrow(prior$strata) != matching_expected_strata ||
      anyDuplicated(prior$strata$stratum_id) ||
      !setequal(prior$strata$stratum_id, current_strata$stratum_id)) {
    return(matching_prior_manifest_failure("prior_manifest_strata_invalid"))
  }
  if (!matching_manifest_complete(prior$strata)) {
    return(matching_prior_manifest_failure("prior_manifest_completion_invalid"))
  }
  if (!matching_covariate_contract_identical(prior$strata, current_strata)) {
    return(matching_prior_manifest_failure("prior_manifest_covariate_contract_mismatch"))
  }

  list(
    valid = TRUE,
    reason = "complete_same_contract_manifest",
    manifest = prior
  )
}

matching_reuse_decision <- function(current_row, validation, prior_proof) {
  stratum_id <- as.character(current_row$stratum_id)
  if (!isTRUE(validation$valid)) {
    return(list(
      reuse = FALSE,
      reason = paste0(
        "invalid_or_missing_output:",
        paste(validation$errors, collapse = "|")
      )
    ))
  }

  if (!stratum_id %in% matching_affected_strata) {
    return(list(reuse = TRUE, reason = "valid_unaffected_output"))
  }

  if (!isTRUE(prior_proof$valid)) {
    return(list(
      reuse = FALSE,
      reason = paste0("affected_output_unattested:", prior_proof$reason)
    ))
  }

  prior_strata <- prior_proof$manifest$strata
  prior_index <- which(prior_strata$stratum_id == stratum_id)
  if (length(prior_index) != 1L) {
    return(list(reuse = FALSE, reason = "affected_output_unattested:stratum_missing"))
  }
  prior_row <- prior_strata[prior_index, , drop = FALSE]

  scalar_fields <- c(
    "country_rast", "country", "biome_raster", "bigregion",
    "source_treated", "source_controls_stratum",
    "source_controls_effective", "source_rows_effective", "pooling_scope",
    "output_file"
  )
  scalar_contract_matches <- all(vapply(scalar_fields, function(field) {
    identical(prior_row[[field]], current_row[[field]])
  }, logical(1)))
  covariate_contract_matches <- identical(
    prior_row$selected_covariates[[1]],
    current_row$selected_covariates[[1]]
  ) && identical(
    prior_row$dropped_covariates[[1]],
    current_row$dropped_covariates[[1]]
  )
  output_contract_matches <- identical(
    prior_row$pair_count,
    validation$pair_count
  ) && identical(
    prior_row$output_md5,
    validation$md5
  ) && identical(
    prior_row$output_size_bytes,
    validation$size_bytes
  )

  if (!scalar_contract_matches ||
      !covariate_contract_matches ||
      !output_contract_matches) {
    return(list(
      reuse = FALSE,
      reason = "affected_output_unattested:stratum_contract_or_fingerprint_mismatch"
    ))
  }

  list(
    reuse = TRUE,
    reason = "valid_affected_output_attested_by_prior_complete_manifest"
  )
}

mark_matching_stratum_failure <- function(strata, index, error_text) {
  existing_error <- strata$error_text[index]
  combined_error <- unique(c(
    if (!is.na(existing_error) && nzchar(existing_error)) {
      existing_error
    },
    error_text
  ))
  combined_validation_errors <- unique(c(
    strata$validation_errors[[index]],
    error_text
  ))

  strata$status[index] <- "failure"
  strata$error_text[index] <- paste(combined_error, collapse = "; ")
  strata$action[index] <- "failed"
  strata$validation_status[index] <- "failed"
  strata$validation_errors[[index]] <- combined_validation_errors
  strata
}

recheck_matching_outputs <- function(
    strata,
    output_dir = matching_output_dir,
    list_files_fun = list.files,
    file_info_fun = file.info,
    md5_fun = matching_file_md5) {
  errors <- character()

  add_stratum_error <- function(index, error_text) {
    stratum_error <- sprintf(
      "Stratum %s final output attestation failed: %s",
      strata$stratum_id[index],
      error_text
    )
    strata <<- mark_matching_stratum_failure(
      strata,
      index,
      error_text
    )
    errors <<- c(errors, stratum_error)
  }

  expected_files <- strata$output_file
  duplicated_files <- unique(expected_files[duplicated(expected_files)])
  if (length(duplicated_files) > 0L) {
    errors <- c(
      errors,
      sprintf(
        "Manifest contains duplicate output paths: %s",
        paste(duplicated_files, collapse = ", ")
      )
    )
    duplicated_indices <- which(expected_files %in% duplicated_files)
    for (index in duplicated_indices) {
      add_stratum_error(index, "duplicate manifest output path")
    }
  }

  if (!dir.exists(output_dir)) {
    actual_files <- character()
    errors <- c(
      errors,
      sprintf("Matching output directory is missing: %s", output_dir)
    )
  } else {
    actual_files <- tryCatch(
      list_files_fun(
        output_dir,
        pattern = "\\.Rds$",
        full.names = TRUE
      ),
      error = identity
    )
    if (inherits(actual_files, "error")) {
      errors <- c(
        errors,
        paste0(
          "Could not enumerate final match-result files: ",
          conditionMessage(actual_files)
        )
      )
      actual_files <- character()
    }
  }

  missing_files <- setdiff(expected_files, actual_files)
  extra_files <- setdiff(actual_files, expected_files)
  if (length(missing_files) > 0L) {
    errors <- c(
      errors,
      sprintf(
        "Final match-result file set is missing: %s",
        paste(basename(missing_files), collapse = ", ")
      )
    )
    for (index in which(expected_files %in% missing_files)) {
      add_stratum_error(index, "output file is missing")
    }
  }
  if (length(extra_files) > 0L) {
    errors <- c(
      errors,
      sprintf(
        "Final match-result file set has unexpected files: %s",
        paste(basename(extra_files), collapse = ", ")
      )
    )
  }

  for (index in seq_len(nrow(strata))) {
    expected_md5 <- strata$output_md5[index]
    expected_size <- strata$output_size_bytes[index]
    if (!matching_md5_is_valid(expected_md5)) {
      add_stratum_error(index, "manifest output MD5 is missing or invalid")
    }
    valid_expected_size <- length(expected_size) == 1L &&
      is.numeric(expected_size) &&
      is.finite(expected_size) &&
      expected_size > 0
    if (!valid_expected_size) {
      add_stratum_error(index, "manifest output size is missing or invalid")
    }

    path <- expected_files[index]
    if (!path %in% actual_files) next

    observed_info <- tryCatch(file_info_fun(path), error = identity)
    if (inherits(observed_info, "error") ||
        nrow(observed_info) != 1L ||
        is.na(observed_info$size) ||
        !is.finite(observed_info$size) ||
        observed_info$size <= 0 ||
        isTRUE(observed_info$isdir)) {
      error_text <- if (inherits(observed_info, "error")) {
        paste0("could not read output metadata: ", conditionMessage(observed_info))
      } else {
        "output metadata or size is invalid"
      }
      add_stratum_error(index, error_text)
      next
    }

    observed_md5 <- tryCatch(md5_fun(path), error = identity)
    if (inherits(observed_md5, "error") ||
        !matching_md5_is_valid(observed_md5)) {
      error_text <- if (inherits(observed_md5, "error")) {
        paste0("could not compute output MD5: ", conditionMessage(observed_md5))
      } else {
        "computed output MD5 is invalid"
      }
      add_stratum_error(index, error_text)
      next
    }

    if (valid_expected_size &&
        as.numeric(observed_info$size) != expected_size) {
      add_stratum_error(
        index,
        sprintf(
          "output size changed from %.0f to %.0f bytes",
          expected_size,
          as.numeric(observed_info$size)
        )
      )
    }
    if (matching_md5_is_valid(expected_md5) &&
        observed_md5 != expected_md5) {
      add_stratum_error(
        index,
        sprintf(
          "output MD5 changed from %s to %s",
          expected_md5,
          observed_md5
        )
      )
    }
  }

  list(
    strata = strata,
    errors = unique(errors)
  )
}

matching_failed_strata_errors <- function(strata) {
  failed <- which(is.na(strata$status) | strata$status != "success")
  if (length(failed) == 0L) return(character())

  vapply(failed, function(index) {
    error_text <- strata$error_text[index]
    if (is.na(error_text) || !nzchar(error_text)) {
      validation_errors <- strata$validation_errors[[index]]
      error_text <- if (length(validation_errors) > 0L) {
        paste(validation_errors, collapse = "; ")
      } else {
        sprintf(
          "status=%s action=%s validation=%s",
          strata$status[index],
          strata$action[index],
          strata$validation_status[index]
        )
      }
    }
    sprintf("Stratum %s failed: %s", strata$stratum_id[index], error_text)
  }, character(1))
}


#################################
### main ########################
#################################

main <- function() {
  t0_total <- Sys.time()
  if (!requireNamespace("Matching", quietly = TRUE)) {
    stop(
      "Package 'Matching' is required for production matching but is not available in this R library"
    )
  }
  cat(sprintf(
    "Matching contract %s | computing full unmatched-input MD5\n",
    matching_contract_version
  ))
  input_fingerprint <- matching_input_fingerprint(paths$data_unmatched)

  d <- load_matching_input(paths$data_unmatched)
  source_strata <- build_expected_strata(d)
  cat(sprintf("Expected strata: %d\n", nrow(source_strata)))

  state <- build_matching_state(d)
  rm(d)
  gc(verbose = FALSE)

  prepared <- prepare_strata_contract(source_strata, state)
  source_strata <- prepared$source_strata
  manifest_strata <- prepared$manifest_strata
  rm(prepared)

  actual_drop_set <- matching_drop_table(manifest_strata)
  if (!matching_drop_set_is_expected(actual_drop_set)) {
    actual_text <- if (nrow(actual_drop_set) == 0L) {
      "none"
    } else {
      paste(
        paste(actual_drop_set$stratum_id, actual_drop_set$covariate, sep = ":"),
        collapse = ", "
      )
    }
    error_text <- sprintf(
      "Near-constant covariate set requires review; observed: %s",
      actual_text
    )
    manifest_strata$error_text <- error_text
    manifest <- build_matching_manifest(
      input_fingerprint,
      manifest_strata,
      errors = error_text
    )
    atomic_save_rds(manifest, matching_manifest_path)
    stop(error_text)
  }

  expected_files <- manifest_strata$output_file
  existing_files <- list.files(
    matching_output_dir,
    pattern = "\\.Rds$",
    full.names = TRUE
  )
  extra_files <- setdiff(existing_files, expected_files)
  if (length(extra_files) > 0L) {
    error_text <- sprintf(
      "Unexpected match-result files: %s",
      paste(basename(extra_files), collapse = ", ")
    )
    manifest_strata$error_text <- error_text
    manifest <- build_matching_manifest(
      input_fingerprint,
      manifest_strata,
      errors = error_text
    )
    atomic_save_rds(manifest, matching_manifest_path)
    stop(error_text)
  }

  prior_proof <- matching_prior_manifest_proof(
    matching_manifest_path,
    input_fingerprint,
    manifest_strata
  )
  cat(sprintf(
    "Prior complete same-contract manifest: %s | %s\n",
    if (prior_proof$valid) "available" else "unavailable",
    prior_proof$reason
  ))

  rerun <- integer()
  for (i in seq_len(nrow(source_strata))) {
    g <- source_strata[i, ]
    effective <- effective_stratum_rows(g, state)
    validation <- validate_match_file(
      manifest_strata$output_file[i],
      g,
      state,
      effective
    )
    reuse_decision <- matching_reuse_decision(
      manifest_strata[i, , drop = FALSE],
      validation,
      prior_proof
    )
    manifest_strata$reuse_reason[i] <- reuse_decision$reason

    if (reuse_decision$reuse) {
      manifest_strata$pair_count[i] <- validation$pair_count
      manifest_strata$status[i] <- "success"
      manifest_strata$action[i] <- "reused"
      manifest_strata$validation_status[i] <- "valid_reused"
      manifest_strata$output_md5[i] <- validation$md5
      manifest_strata$output_size_bytes[i] <- validation$size_bytes
      manifest_strata$validation_errors[[i]] <- character()
    } else {
      rerun <- c(rerun, i)
      manifest_strata$validation_errors[[i]] <- if (validation$valid) {
        character()
      } else {
        validation$errors
      }
    }
  }

  cat(sprintf(
    paste0(
      "Validated reusable strata: %d | To match/rematch: %d ",
      "(affected contract strata: %d)\n"
    ),
    nrow(source_strata) - length(rerun),
    length(rerun),
    sum(manifest_strata$stratum_id[rerun] %in% matching_affected_strata)
  ))

  worker_errors <- character()
  if (length(rerun) > 0L) {
    core_override <- Sys.getenv("PATHREAT_MATCHING_CORES", unset = "")
    if (nzchar(core_override)) {
      detected_cores <- suppressWarnings(as.integer(core_override))
      if (length(detected_cores) != 1L || is.na(detected_cores) ||
          detected_cores < 1L) {
        stop("PATHREAT_MATCHING_CORES must be a positive integer")
      }
    } else {
      detected_cores <- suppressWarnings(parallel::detectCores())
      if (length(detected_cores) != 1L ||
          is.na(detected_cores) || detected_cores < 1L) {
        # The production machine is configured for the existing four-worker
        # matching design, but sandboxed macOS sessions can return NA here.
        detected_cores <- 4L
      }
    }
    no_cluster <- min(4L, as.integer(detected_cores), length(rerun))
    cluster <- if (no_cluster > 1L) no_cluster else NULL
    if (!is.null(cluster)) {
      cat(sprintf("Running matching with %d cores\n", cluster))
    }

    raw_results <- tryCatch(
      pblapply(
        rerun,
        function(i) {
          match_one_stratum(i, source_strata, manifest_strata, state)
        },
        cl = cluster
      ),
      error = identity
    )
    normalized_results <- normalize_matching_worker_results(
      raw_results,
      rerun
    )
    results <- normalized_results$results
    worker_errors <- normalized_results$errors

    for (result in results) {
      i <- result$index
      if (result$success) {
        manifest_strata$pair_count[i] <- result$pair_count
        manifest_strata$status[i] <- "success"
        manifest_strata$error_text[i] <- NA_character_
        manifest_strata$action[i] <- "rematched"
        manifest_strata$validation_status[i] <- "new_output_valid"
        manifest_strata$output_md5[i] <- result$output_md5
        manifest_strata$output_size_bytes[i] <- result$output_size_bytes
        manifest_strata$validation_errors[[i]] <- character()
      } else {
        manifest_strata$status[i] <- "failure"
        manifest_strata$error_text[i] <- result$error_text
        manifest_strata$action[i] <- "failed"
        manifest_strata$validation_status[i] <- "failed"
        manifest_strata$validation_errors[[i]] <- result$error_text
      }
    }
  }

  input_attestation <- recheck_matching_input(
    paths$data_unmatched,
    input_fingerprint
  )
  output_attestation <- recheck_matching_outputs(
    manifest_strata,
    matching_output_dir
  )
  manifest_strata <- output_attestation$strata

  global_errors <- unique(c(
    worker_errors,
    input_attestation$errors,
    output_attestation$errors,
    matching_failed_strata_errors(manifest_strata)
  ))
  if (!matching_manifest_complete(manifest_strata) &&
      length(global_errors) == 0L) {
    global_errors <- paste0(
      "Final matching manifest contract failed without a ",
      "stratum-specific diagnostic"
    )
  }

  manifest <- build_matching_manifest(
    input_fingerprint,
    manifest_strata,
    errors = global_errors
  )
  atomic_save_rds(manifest, matching_manifest_path)

  if (!manifest$complete) {
    stop(sprintf(
      "Matching incomplete: %s",
      paste(manifest$errors, collapse = "; ")
    ))
  }

  cat(sprintf(
    paste0(
      "Matching complete: %d/%d strata, %s pairs. ",
      "Manifest: %s. Total time: %.1f min\n"
    ),
    nrow(manifest_strata),
    matching_expected_strata,
    format(sum(manifest_strata$pair_count), big.mark = ","),
    matching_manifest_path,
    as.numeric(difftime(Sys.time(), t0_total, units = "mins"))
  ))
}


helper_only <- identical(Sys.getenv("PATHREAT_MATCHING_HELPERS_ONLY"), "1") ||
  "--helpers-only" %in% commandArgs(trailingOnly = TRUE)

if (sys.nframe() == 0L && !helper_only) {
  main()
}
