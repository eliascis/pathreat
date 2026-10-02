##################################################
### pathreat.data.species-count.R ################
### Pixel-level threatened species counts ########
### and phylogenetic score sums ##################
##################################################
#
# Counts threatened Red List species coverage at the pixel level using the
# rasterized species ranges. Produces one total count and one count per taxon,
# plus pixel-level sums of TBL, ED, and EDGE scores across all present species,
# row-aligned with merge.fst for downstream cbind merges into merge.unmatched.
#
# Species with NA phylogenetic scores (unmatched to EDGE dataset) contribute
# to species counts but are skipped in the score sums.
#
# Two-phase approach:
#   Phase 1 (parallel): Extract presence indices per species raster, save
#           sparse index vectors to data/tmp/species_count/sp/.
#   Phase 2 (sequential): Aggregate saved indices into pixel-level counts
#           and score sums.
#
# Inputs:
#   data/store/pathreat.data.merge.fst
#   data/store/pathreat.data.species.Rds
#   data/store/pathreat.data.rasterbase.tif
#   data/store/species_rasters/sp_{taxon}_{id}.tif
#
# Output:
#   data/store/pathreat.data.species-count.Rds
#     Columns: pixel_id, n_redlist_species, n_redlist_species_{taxon},
#              tbl_sum, ed_sum, edge_sum
#

library(fst)
library(terra)
library(parallel)

source("code/pathreat.analysis.config.R")

taxa <- c("amphibian", "bird", "mammal", "reptile")
count_cols <- c(
  total = "n_redlist_species",
  amphibian = "n_redlist_species_amphibian",
  bird = "n_redlist_species_bird",
  mammal = "n_redlist_species_mammal",
  reptile = "n_redlist_species_reptile"
)
score_cols <- c("tbl_sum", "ed_sum", "edge_sum")

n_cores <- min(detectCores() - 1, 8)
cat(sprintf("Using %d cores for parallel extraction\n", n_cores))

############################
### load pixel geometry ####
############################
{
cat("=== Loading merged pixel geometry ===\n")

d_xy <- read_fst(
  "data/store/pathreat.data.merge.fst",
  columns = c("pixel_id", "x", "y")
)

n_pixels <- nrow(d_xy)
cat(sprintf("  Pixels: %s\n", format(n_pixels, big.mark = ",")))

template <- rast("data/store/pathreat.data.rasterbase.tif")
cell_idx <- cellFromXY(template, cbind(d_xy$x, d_xy$y))
valid_idx <- which(!is.na(cell_idx))

cat(sprintf("  Valid raster cells: %s\n", format(length(valid_idx), big.mark = ",")))
if (length(valid_idx) < n_pixels) {
  cat(sprintf("  Pixels with NA cell index: %s\n",
              format(n_pixels - length(valid_idx), big.mark = ",")))
}

d_xy$x <- NULL
d_xy$y <- NULL
rm(template)
gc()
}

############################
### load species metadata ###
############################
{
cat("\n=== Loading species metadata ===\n")

sp_meta <- readRDS("data/store/pathreat.data.species.Rds")
sp_meta <- sp_meta[, c("species_id", "taxon", "raster_file",
                        "tbl_median", "ed_median", "edge_median")]
sp_meta <- sp_meta[sp_meta$taxon %in% taxa, ]
sp_meta <- sp_meta[order(match(sp_meta$taxon, taxa), sp_meta$species_id), ]
rownames(sp_meta) <- NULL

cat(sprintf("  Species: %s\n", format(nrow(sp_meta), big.mark = ",")))
for (tx in taxa) {
  cat(sprintf("    %-12s %s\n",
              tx,
              format(sum(sp_meta$taxon == tx), big.mark = ",")))
}

n_phylo_na <- sum(is.na(sp_meta$tbl_median))
cat(sprintf("  Species with phylogenetic scores: %s / %s (%.1f%%)\n",
            format(nrow(sp_meta) - n_phylo_na, big.mark = ","),
            format(nrow(sp_meta), big.mark = ","),
            100 * (1 - n_phylo_na / nrow(sp_meta))))

sp_dir <- "data/store/species_rasters"
sp_meta$raster_path <- file.path(sp_dir, sp_meta$raster_file)

missing_rasters <- sp_meta$raster_path[!file.exists(sp_meta$raster_path)]
if (length(missing_rasters) > 0) {
  stop(sprintf(
    "Missing %d species rasters; first missing file: %s",
    length(missing_rasters),
    missing_rasters[1]
  ))
}
}

###########################################
### phase 1: parallel raster extraction ###
###########################################
{
cat("\n=== Phase 1: Parallel raster extraction ===\n")

tmp_sp_dir <- "data/tmp/species_count/sp"
dir.create(tmp_sp_dir, showWarnings = FALSE, recursive = TRUE)

already_done <- list.files(tmp_sp_dir, pattern = "^sp_\\d+\\.rds$")
done_ids <- as.integer(sub("^sp_(\\d+)\\.rds$", "\\1", already_done))
todo <- sp_meta[!sp_meta$species_id %in% done_ids, ]

cat(sprintf("  Already extracted: %s / %s\n",
            format(nrow(sp_meta) - nrow(todo), big.mark = ","),
            format(nrow(sp_meta), big.mark = ",")))

if (nrow(todo) > 0) {
  cat(sprintf("  Extracting %s species on %d cores...\n",
              format(nrow(todo), big.mark = ","),
              n_cores))

  cell_idx_valid <- cell_idx[valid_idx]

  extract_one <- function(i) {
    sp_id <- todo$species_id[i]
    out_file <- file.path(tmp_sp_dir, paste0("sp_", sp_id, ".rds"))

    tryCatch({
      sp_rast <- terra::rast(todo$raster_path[i])
      sp_vals <- sp_rast[cell_idx_valid][, 1]
      present_idx <- valid_idx[sp_vals == 1]
      saveRDS(present_idx, out_file)
      TRUE
    }, error = function(e) {
      warning(sprintf("  ERROR species %d: %s", sp_id, e$message))
      FALSE
    })
  }

  batch_size <- 500
  n_batches <- ceiling(nrow(todo) / batch_size)

  for (b in seq_len(n_batches)) {
    idx_start <- (b - 1) * batch_size + 1
    idx_end <- min(b * batch_size, nrow(todo))
    batch_idx <- idx_start:idx_end

    results <- mclapply(
      batch_idx,
      extract_one,
      mc.cores = n_cores
    )

    n_ok <- sum(unlist(results), na.rm = TRUE)
    cat(sprintf("    Batch %d/%d: %d/%d species extracted\n",
                b,
                n_batches,
                n_ok,
                length(batch_idx)))
  }

  rm(cell_idx_valid)
  gc()
}

extracted_files <- list.files(tmp_sp_dir, pattern = "^sp_\\d+\\.rds$", full.names = TRUE)
extracted_ids <- as.integer(sub("^sp_(\\d+)\\.rds$", "\\1", basename(extracted_files)))

missing_sp <- setdiff(sp_meta$species_id, extracted_ids)
if (length(missing_sp) > 0) {
  stop(sprintf("Missing extraction results for %d species. Check errors above.",
               length(missing_sp)))
}

cat(sprintf("  All %s species extracted successfully\n",
            format(length(extracted_files), big.mark = ",")))
}

##########################################
### phase 2: aggregate into pixel sums ###
##########################################
{
cat("\n=== Phase 2: Aggregating pixel-level counts and scores ===\n")

count_by_taxon <- lapply(taxa, function(tx) integer(n_pixels))
names(count_by_taxon) <- taxa

tbl_total <- numeric(n_pixels)
ed_total <- numeric(n_pixels)
edge_total <- numeric(n_pixels)

sp_lookup <- sp_meta[, c("species_id", "taxon", "tbl_median", "ed_median", "edge_median")]
rownames(sp_lookup) <- sp_lookup$species_id

for (i in seq_len(nrow(sp_meta))) {
  sp_id <- sp_meta$species_id[i]
  tx <- sp_meta$taxon[i]

  present_idx <- readRDS(file.path(tmp_sp_dir, paste0("sp_", sp_id, ".rds")))

  if (length(present_idx) > 0) {
    count_by_taxon[[tx]][present_idx] <- count_by_taxon[[tx]][present_idx] + 1L

    if (!is.na(sp_meta$tbl_median[i])) {
      tbl_total[present_idx] <- tbl_total[present_idx] + sp_meta$tbl_median[i]
    }
    if (!is.na(sp_meta$ed_median[i])) {
      ed_total[present_idx] <- ed_total[present_idx] + sp_meta$ed_median[i]
    }
    if (!is.na(sp_meta$edge_median[i])) {
      edge_total[present_idx] <- edge_total[present_idx] + sp_meta$edge_median[i]
    }
  }

  if (i %% 500 == 0 || i == nrow(sp_meta)) {
    cat(sprintf("    Aggregated %s / %s species\n",
                format(i, big.mark = ","),
                format(nrow(sp_meta), big.mark = ",")))
  }
}

rm(sp_lookup)
gc()
}

############################
### assemble + verify ######
############################
{
cat("\n=== Assembling output ===\n")

count_total <- count_by_taxon[[taxa[1]]]
if (length(taxa) > 1) {
  for (tx in taxa[-1]) {
    count_total <- count_total + count_by_taxon[[tx]]
  }
}

d_count <- data.frame(
  pixel_id = d_xy$pixel_id,
  n_redlist_species = count_total,
  n_redlist_species_amphibian = count_by_taxon$amphibian,
  n_redlist_species_bird = count_by_taxon$bird,
  n_redlist_species_mammal = count_by_taxon$mammal,
  n_redlist_species_reptile = count_by_taxon$reptile,
  tbl_sum = tbl_total,
  ed_sum = ed_total,
  edge_sum = edge_total
)

sum_check <- d_count$n_redlist_species_amphibian +
  d_count$n_redlist_species_bird +
  d_count$n_redlist_species_mammal +
  d_count$n_redlist_species_reptile

if (!all(d_count$n_redlist_species == sum_check)) {
  stop("Total red-listed species count does not equal the sum of taxon counts.")
}

cat("  Count summaries:\n")
for (v in count_cols) {
  cat(sprintf("    %-28s min=%d median=%d mean=%.2f max=%d\n",
              v,
              min(d_count[[v]], na.rm = TRUE),
              stats::median(d_count[[v]], na.rm = TRUE),
              mean(d_count[[v]], na.rm = TRUE),
              max(d_count[[v]], na.rm = TRUE)))
}

cat("  Score summaries:\n")
for (v in score_cols) {
  cat(sprintf("    %-28s min=%.2f median=%.2f mean=%.2f max=%.2f\n",
              v,
              min(d_count[[v]], na.rm = TRUE),
              stats::median(d_count[[v]], na.rm = TRUE),
              mean(d_count[[v]], na.rm = TRUE),
              max(d_count[[v]], na.rm = TRUE)))
}
}

############################
### save ###################
############################
{
outfile <- "data/store/pathreat.data.species-count.Rds"
saveRDS(d_count, outfile)

cat(sprintf("\nSaved: %s (%.1f MB)\n",
            outfile,
            file.size(outfile) / 1e6))
}
