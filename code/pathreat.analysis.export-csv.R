##########################################
### pathreat.analysis.export-csv.R ######
### CSV twins of the shipped .Rds     ####
### estimate tables + results manifest ###
##########################################
# Writes a same-root .csv next to each canonical .Rds estimate table so the
# public release is readable without R. Data frames become one CSV; named
# lists of data frames (e.g. bytaxa.est) become one CSV per element
# (<stem>.<element>.csv); anything else is skipped with a message.
# Also exports the per-PA estimates from data/store/ into results/ and writes
# results/MANIFEST.csv (file, rows, cols, md5) for the shipped tables.

library(tools)
library(dplyr)

source("code/pathreat.analysis.config.R")

############################
### files to export ########
############################

rds_files <- c(
  "results/pathreat.global.est.Rds",
  "results/pathreat.global.group_means.est.Rds",
  "results/pathreat.global.fe-adjusted-levels.est.Rds",
  "results/pathreat.global-threat.est.Rds",
  "results/pathreat.bycountry.est.Rds",
  "results/pathreat.bybiome.est.Rds",
  "results/pathreat.hotspots.est.Rds",
  "results/pathreat.analysis.hotspots.Rds",
  "results/pathreat.PA-type.est.Rds",
  "results/pathreat.PA-type.group_means.est.Rds",
  "results/pathreat.byspecies.est.Rds",
  "results/pathreat.bytaxa.est.Rds",
  "results/pathreat.sfa.est.Rds",
  "results/pathreat.sfa.country.est.Rds",
  "results/pathreat.sfa.aggregates.Rds",
  "results/pathreat.randomrob.est.Rds",
  "results/pathreat.robustness.est.Rds",
  "results/pathreat.sensitivity.evalues.Rds",
  "results/pathreat.sensitivity.rosenbaum.Rds",
  "results/pathreat.covbalance.by_country.Rds",
  "results/pathreat.threat_composite_weights.Rds",
  "data/store/pathreat.analysis.byPA.est.Rds"
)

# WDPA attribute columns copied from the source database; dropped from the
# CSV export of the PA-level SFA table (wdpaid remains for rejoining)
wdpa_attribute_cols <- c("size", "log_size", "size_class", "iucn_class")

# CSV files written by their own producer script; never overwritten here
producer_owned_csv <- c(
  "results/pathreat.bycountry.est.csv",
  "results/pathreat.global.fe-adjusted-levels.est.csv",
  "results/pathreat.threat_composite_weights.csv",
  "results/pathreat.sfa.aggregates.csv"
)

############################
### helpers ################
############################

csv_target <- function(rds_file, element = NULL) {
  stem <- sub("\\.Rds$", "", basename(rds_file))
  if (!is.null(element)) {
    stem <- paste0(stem, ".", element)
  }
  paste0(paths$results_dir, stem, ".csv")
}

flatten_list_cols <- function(df) {
  # list columns cannot be written to CSV; collapse to "a; b; c"
  is_list <- vapply(df, is.list, logical(1))
  for (col in names(df)[is_list]) {
    df[[col]] <- vapply(
      df[[col]],
      function(x) paste(as.character(unlist(x)), collapse = "; "),
      character(1)
    )
  }
  df
}

write_one <- function(df, out_file, drop_cols = character(0)) {
  if (out_file %in% producer_owned_csv) {
    cat(sprintf("  %-70s producer-owned, kept as is\n", out_file))
    return(out_file)
  }
  df <- as.data.frame(df)
  df <- df[, setdiff(names(df), drop_cols), drop = FALSE]
  df <- flatten_list_cols(df)
  write.csv(df, out_file, row.names = FALSE)
  cat(sprintf("  %-70s %6d rows %3d cols\n", out_file, nrow(df), ncol(df)))
  out_file
}

export_rds <- function(rds_file) {
  if (!file.exists(rds_file)) {
    cat("  MISSING:", rds_file, "\n")
    return(character(0))
  }
  obj <- readRDS(rds_file)
  drop_cols <- character(0)
  if (basename(rds_file) == "pathreat.sfa.est.Rds") {
    drop_cols <- wdpa_attribute_cols
  }
  if (is.data.frame(obj)) {
    return(write_one(obj, csv_target(rds_file), drop_cols))
  }
  if (is.list(obj) && !is.null(names(obj))) {
    df_elements <- names(obj)[vapply(obj, is.data.frame, logical(1))]
    if (length(df_elements) == 0) {
      cat("  SKIP (list without data frames):", rds_file, "\n")
      return(character(0))
    }
    written <- lapply(df_elements, function(el) {
      write_one(obj[[el]], csv_target(rds_file, el))
    })
    skipped <- setdiff(names(obj), df_elements)
    if (length(skipped) > 0) {
      cat("  (non-table elements not exported:", paste(skipped, collapse = ", "), ")\n")
    }
    return(unlist(written))
  }
  cat("  SKIP (not a table):", rds_file, "\n")
  character(0)
}

############################
### export #################
############################

cat("Exporting CSV twins...\n")
written <- lapply(rds_files, export_rds)
written <- unlist(written)

############################
### manifest ###############
############################

manifest_files <- sort(unique(c(
  rds_files[file.exists(rds_files)],
  written
)))

manifest_rows <- lapply(manifest_files, function(f) {
  rows <- NA_integer_
  cols <- NA_integer_
  if (grepl("\\.csv$", f)) {
    df <- read.csv(f)
    cols <- ncol(df)
    rows <- nrow(df)
  }
  data.frame(
    file = f,
    rows = rows,
    cols = cols,
    bytes = file.size(f),
    md5 = unname(md5sum(f)),
    stringsAsFactors = FALSE
  )
})
manifest <- do.call(rbind, manifest_rows)

manifest_file <- paste0(paths$results_dir, "MANIFEST.csv")
write.csv(manifest, manifest_file, row.names = FALSE)
cat(sprintf("\nWrote %s (%d files)\n", manifest_file, nrow(manifest)))
