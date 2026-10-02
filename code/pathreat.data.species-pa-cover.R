##################################################
### pathreat.data.species-pa-cover.R #############
### Species x PA overlap area (intersection) #####
### All four vertebrate taxa #####################
##################################################

library(fst)
library(sf)
library(terra)
source("code/pathreat.analysis.config.R")


############################
### load data ##############
############################
{
cat("=== Loading pixel coordinates and WDPA ===\n")

d_xy <- read_fst("data/store/pathreat.data.prepdata.fst",
                  columns = c("x", "y"))
d_wdpa <- read_fst("data/store/pathreat.data.wdpa.fst",
                    columns = c("wdpaid"))

stopifnot(nrow(d_xy) == nrow(d_wdpa))
cat("  Pixels:", format(nrow(d_xy), big.mark = ","), "\n")

# NA wdpaid → 0 (unprotected)
wdpaid_raw <- d_wdpa$wdpaid
wdpaid_raw[is.na(wdpaid_raw)] <- 0L

rm(d_wdpa)
gc()
}


############################
### cell indices ###########
############################
{
cat("\n=== Computing cell indices ===\n")

template <- rast("data/store/pathreat.data.rasterbase.tif")
cell_idx <- cellFromXY(template, cbind(d_xy$x, d_xy$y))

# drop invalid cell indices
valid <- !is.na(cell_idx)
if (any(!valid)) {
  cat(sprintf("  Dropped %d pixels with NA cell index\n", sum(!valid)))
  cell_idx   <- cell_idx[valid]
  wdpaid_vec <- wdpaid_raw[valid]
} else {
  wdpaid_vec <- wdpaid_raw
}

rm(d_xy, wdpaid_raw, valid)
gc()

cat("  Valid pixels:", format(length(cell_idx), big.mark = ","), "\n")
cat("  Unique wdpaids:", format(length(unique(wdpaid_vec)), big.mark = ","),
    "(incl. 0 = unprotected)\n")
}


############################
### species metadata #######
############################
{
cat("\n=== Loading species metadata ===\n")
sp_meta <- readRDS("data/store/pathreat.data.species.Rds")
cat(sprintf("  Loaded %d species (%s)\n",
            nrow(sp_meta),
            paste(names(table(sp_meta$taxon)), collapse = ", ")))
}


############################
### enumerate rasters ######
############################
{
cat("\n=== Enumerating species rasters ===\n")

sp_dir   <- "data/store/species_rasters"
sp_files <- list.files(sp_dir, pattern = "^sp_.*\\.tif$", full.names = TRUE)

# Parse taxon and id from filename: sp_{taxon}_{id}.tif
bn <- basename(sp_files)
parsed <- regmatches(bn, regexec("^sp_(mammal|amphibian|reptile|bird)_(\\d+)\\.tif$", bn))

# Keep only validly named files
valid_parse <- sapply(parsed, length) == 3
if (any(!valid_parse)) {
  cat(sprintf("  WARNING: %d files with unrecognized naming pattern (skipped)\n",
              sum(!valid_parse)))
}
sp_files <- sp_files[valid_parse]
parsed   <- parsed[valid_parse]

all_sp_taxon <- sapply(parsed, `[`, 2)
all_sp_ids   <- as.integer(sapply(parsed, `[`, 3))

# tmp directory for resumable per-species results
tmp_dir <- "data/tmp/species_pa_coverage"
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

# check which species already done (use taxon_id as key)
done_files <- list.files(tmp_dir, pattern = "^sp_.*\\.rds$")
done_keys  <- sub("\\.rds$", "", done_files)
all_keys   <- paste0("sp_", all_sp_taxon, "_", all_sp_ids)
todo_idx   <- which(!all_keys %in% done_keys)

cat(sprintf("  Found %d species rasters (%d already done, %d to process)\n",
            length(sp_files), length(done_files), length(todo_idx)))

# Breakdown by taxon
for (tx in c("mammal", "amphibian", "reptile", "bird")) {
  n_total <- sum(all_sp_taxon == tx)
  n_todo  <- sum(all_sp_taxon[todo_idx] == tx)
  cat(sprintf("    %-12s %d total, %d to process\n", tx, n_total, n_todo))
}
}


############################
### parallel processing ####
############################
if (length(todo_idx) > 0) {

n_cores <- no_cluster
cat(sprintf("\n=== Processing %d species on %d cores ===\n",
            length(todo_idx), n_cores))

parallel::mclapply(seq_along(todo_idx), function(j) {
  # j<-seq_along(todo_idx)[1]
  i     <- todo_idx[j]
  sp_id <- all_sp_ids[i]
  sp_tx <- all_sp_taxon[i]
  tmp_file <- file.path(tmp_dir, paste0("sp_", sp_tx, "_", sp_id, ".rds"))

  sp_rast    <- rast(sp_files[i])
  sp_present <- sp_rast[cell_idx][, 1]

  in_range <- which(sp_present == 1)
  n_range  <- length(in_range)

  if (n_range == 0) {
    saveRDS(data.frame(
      species_id  = integer(0),
      taxon       = character(0),
      wdpaid      = integer(0),
      overlap_area = numeric(0),
      range_total = numeric(0),
      stringsAsFactors = FALSE
    ), tmp_file)
    return(NULL)
  }

  # tabulate wdpaid within species range
  wdpa_in_range <- wdpaid_vec[in_range]
  tab <- as.data.frame(table(wdpa_in_range), stringsAsFactors = FALSE)
  names(tab) <- c("wdpaid", "overlap_area")
  tab$wdpaid       <- as.integer(tab$wdpaid)
  tab$overlap_area <- as.numeric(tab$overlap_area)

  tab$species_id  <- sp_id
  tab$taxon       <- sp_tx
  tab$range_total <- as.numeric(n_range)

  tab <- tab[, c("species_id", "taxon", "wdpaid", "overlap_area", "range_total")]

  saveRDS(tab, tmp_file)

  if (j %% 100 == 0)
    cat(sprintf("  [%d/%d] sp_%s_%d: %s pixels, %d PAs\n",
                j, length(todo_idx), sp_tx, sp_id,
                format(n_range, big.mark = ","), nrow(tab) - 1L))
  return(NULL)
}, mc.cores = n_cores, mc.preschedule = FALSE)

rm(cell_idx, wdpaid_vec, template)
gc()

} # end if (todo_idx > 0)


############################
### collect + assemble #####
############################
{
cat("\n=== Collecting results ===\n")

tmp_files <- list.files(tmp_dir, pattern = "^sp_.*\\.rds$", full.names = TRUE)
species_pa <- do.call(rbind, lapply(tmp_files, readRDS))
rownames(species_pa) <- NULL

# handle legacy tmp files with old column name
if ("area" %in% names(species_pa) && !"overlap_area" %in% names(species_pa)) {
  names(species_pa)[names(species_pa) == "area"] <- "overlap_area"
}

# drop empty species (zero-row tmp files)
cat(sprintf("  Raw rows: %s\n", format(nrow(species_pa), big.mark = ",")))
species_pa <- species_pa[!is.na(species_pa$species_id), ]

# Handle old tmp files that lack taxon column (from mammal-only runs)
if (!"taxon" %in% names(species_pa)) {
  species_pa$taxon <- "mammal"
}
species_pa$taxon[is.na(species_pa$taxon)] <- "mammal"

# join species names and category from metadata
species_pa <- merge(species_pa,
                    sp_meta[, c("species_id", "taxon", "species_name", "category")],
                    by = c("species_id", "taxon"),
                    all.x = TRUE)

# reorder columns and sort
species_pa <- species_pa[, c("species_id", "species_name", "taxon", "category",
                              "wdpaid", "overlap_area", "range_total")]
species_pa <- species_pa[order(species_pa$taxon,
                               species_pa$species_id,
                               -species_pa$overlap_area), ]
rownames(species_pa) <- NULL
}


############################
### verify + save ##########
############################
{
cat("\n=== Verification ===\n")

# check area sums = range_total for every species
area_check <- species_pa %>%
  group_by(species_id, taxon) %>%
  dplyr::summarize(
    sum_area    = sum(overlap_area),
    range_total = first(range_total),
    .groups = "drop"
  ) %>%
  mutate(diff = abs(sum_area - range_total))

n_mismatch <- sum(area_check$diff > 0)
if (n_mismatch > 0) {
  cat(sprintf("  WARNING: %d species with area sum != range_total\n", n_mismatch))
} else {
  cat("  OK: area sums match range_total for all species\n")
}

# species count by taxon
sp_counts <- species_pa %>%
  distinct(species_id, taxon) %>%
  group_by(taxon) %>%
  dplyr::summarize(n = n(), .groups = "drop")
cat("  Species by taxon:\n")
for (i in seq_len(nrow(sp_counts))) {
  cat(sprintf("    %-12s %d\n", sp_counts$taxon[i], sp_counts$n[i]))
}
cat(sprintf("    %-12s %d\n", "TOTAL", sum(sp_counts$n)))

# missing names
n_missing_name <- species_pa %>%
  distinct(species_id, taxon, species_name) %>%
  filter(is.na(species_name)) %>%
  nrow()
if (n_missing_name > 0) {
  cat(sprintf("  WARNING: %d species with missing names\n", n_missing_name))
} else {
  cat("  OK: all species have names\n")
}

# missing category
n_missing_cat <- species_pa %>%
  distinct(species_id, taxon, category) %>%
  filter(is.na(category)) %>%
  nrow()
if (n_missing_cat > 0) {
  cat(sprintf("  WARNING: %d species with missing IUCN category\n", n_missing_cat))
} else {
  cat("  OK: all species have IUCN category\n")
}

# summary
cat(sprintf("  Total rows: %s\n", format(nrow(species_pa), big.mark = ",")))
cat(sprintf("  Unique PAs (excl. 0): %d\n",
            length(unique(species_pa$wdpaid[species_pa$wdpaid != 0]))))


d<-species_pa
d$species_name<-NULL
d$species_taxon<-NULL
d$taxon<-NULL
d$category<-NULL


# save
outfile <- "data/store/pathreat.data.species-pa-cover.Rds"
saveRDS(d, outfile)
cat(sprintf("\nSaved: %s (%.1f MB)\n",
            outfile, file.size(outfile) / 1e6))
}

head(d)


