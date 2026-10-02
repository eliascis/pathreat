##########################################
### pathreat.data.species.raster.R ######
### Read IUCN species range shapefiles ##
### All four vertebrate taxa ############
##########################################

library(sf)
library(terra)
library(dplyr)
library(readxl)

# --- Paths ---
ext_dir <- "data/external"
out_dir <- "data/store/species_rasters"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Grid template (1-km Mollweide) ---
grid_template <- rast("data/store/pathreat.data.rasterbase.tif")
mollweide_crs <- crs(grid_template)


############################
### migrate mammal names ###
############################
{
# One-time rename: sp_{id}.tif → sp_mammal_{id}.tif
old_mammal_files <- list.files(out_dir, pattern = "^sp_\\d+\\.tif$", full.names = TRUE)
if (length(old_mammal_files) > 0) {
  cat("=== Migrating", length(old_mammal_files), "mammal rasters to new naming convention ===\n")
  for (f in old_mammal_files) {
    new_name <- sub("^sp_(\\d+)\\.tif$", "sp_mammal_\\1.tif", basename(f))
    file.rename(f, file.path(out_dir, new_name))
  }
  cat("  Done.\n")
}
}


############################
### taxa configurations ####
############################
{
# Define how to read each taxon's data
# For shapefiles: read PART1 + PART2, filter by category column
# For birds: read GPKG, join Red List category from Excel checklist

taxa <- list(
  mammal = list(
    type = "shp",
    dir = "MAMMALS",
    parts = c("MAMMALS_PART1.shp", "MAMMALS_PART2.shp"),
    id_col = "id_no",
    name_col = "sci_name",
    cat_col = "category"
  ),
  amphibian = list(
    type = "shp",
    dir = "AMPHIBIANS",
    parts = c("AMPHIBIANS_PART1.shp", "AMPHIBIANS_PART2.shp"),
    id_col = "id_no",
    name_col = "sci_name",
    cat_col = "category"
  ),
  reptile = list(
    type = "shp",
    dir = "REPTILES",
    parts = c("REPTILES_PART1.shp", "REPTILES_PART2.shp"),
    id_col = "id_no",
    name_col = "sci_name",
    cat_col = "category"
  ),
  bird = list(
    type = "gpkg",
    dir = "BIRDS",
    file = "BOTW_2025.gpkg",
    layer = "all_species",
    id_col = "sisid",
    name_col = "sci_name",
    checklist = "Handbook of the Birds of the World and BirdLife International Digital Checklist of the Birds of the World_Version_10.xlsx"
  )
)
}


############################
### process each taxon #####
############################

for (taxon_name in names(taxa)) {

  cfg <- taxa[[taxon_name]]
  prefix <- paste0("sp_", taxon_name, "_")

  cat(sprintf("\n========================================\n"))
  cat(sprintf("=== Processing: %s\n", toupper(taxon_name)))
  cat(sprintf("========================================\n"))

  ############################
  ### read + filter ##########
  ############################

  if (cfg$type == "shp") {
    # Read shapefile parts and bind
    parts <- lapply(cfg$parts, function(p) {
      st_read(file.path(ext_dir, cfg$dir, p), quiet = TRUE)
    })
    sp_data <- bind_rows(parts)
    rm(parts)

    cat(sprintf("  Combined features: %d\n", nrow(sp_data)))
    cat(sprintf("  Categories: %s\n",
                paste(sort(unique(sp_data[[cfg$cat_col]])), collapse = ", ")))

    # Filter to threatened (CR, EN, VU)
    sp_data <- sp_data %>%
      filter(.data[[cfg$cat_col]] %in% c("CR", "EN", "VU"))

    cat(sprintf("  Threatened features: %d\n", nrow(sp_data)))

    # Standardize columns
    sp_data$species_id <- sp_data[[cfg$id_col]]
    sp_data$sci_name <- sp_data[[cfg$name_col]]

  } else if (cfg$type == "gpkg") {
    # Read GeoPackage
    sp_data <- st_read(file.path(ext_dir, cfg$dir, cfg$file),
                       layer = cfg$layer, quiet = TRUE)
    cat(sprintf("  GPKG features: %d\n", nrow(sp_data)))

    # Read Red List categories from Excel checklist
    checklist <- read_excel(
      file.path(ext_dir, cfg$dir, cfg$checklist),
      sheet = 1,
      skip = 3
    )
    checklist <- checklist[!is.na(checklist$SISRecID), ]
    checklist$category <- checklist[["2025 IUCN Red List category"]]
    checklist$SISRecID <- as.integer(checklist$SISRecID)

    # Map CR(PE) and CR(PEW) → CR
    checklist$category <- sub("^CR \\(PE\\)$", "CR", checklist$category)
    checklist$category <- sub("^CR \\(PEW\\)$", "CR", checklist$category)

    cat(sprintf("  Checklist species: %d\n", nrow(checklist)))

    # Join category to spatial data
    sp_data <- merge(sp_data,
                     checklist[, c("SISRecID", "category")],
                     by.x = cfg$id_col,
                     by.y = "SISRecID",
                     all.x = TRUE)
    rm(checklist)

    cat(sprintf("  Categories: %s\n",
                paste(sort(unique(sp_data$category)), collapse = ", ")))

    # Filter to threatened
    sp_data <- sp_data %>%
      filter(category %in% c("CR", "EN", "VU"))

    cat(sprintf("  Threatened features: %d\n", nrow(sp_data)))

    # Standardize columns
    sp_data$species_id <- sp_data[[cfg$id_col]]
    sp_data$sci_name <- sp_data[[cfg$name_col]]
  }

  # Summary
  n_species <- length(unique(sp_data$species_id))
  cat(sprintf("  Unique threatened species: %d\n", n_species))

  ############################
  ### reproject ###############
  ############################

  sp_vect <- vect(sp_data)
  rm(sp_data)
  gc()
  sp_vect <- project(sp_vect, mollweide_crs)
  cat("  Reprojected to Mollweide\n")

  ############################
  ### check existing #########
  ############################

  species_ids <- unique(sp_vect$species_id)
  pattern <- paste0("^", prefix, "\\d+\\.tif$")
  existing <- list.files(out_dir, pattern = pattern)
  existing_ids <- as.integer(sub(
    paste0("^", prefix, "(\\d+)\\.tif$"), "\\1", existing
  ))
  todo_ids <- setdiff(species_ids, existing_ids)

  cat(sprintf("  Total: %d | Already done: %d | Remaining: %d\n",
              length(species_ids), length(existing_ids), length(todo_ids)))

  ############################
  ### rasterize ##############
  ############################

  if (length(todo_ids) == 0) {
    cat("  All species already rasterized. Skipping.\n")
  } else {
    n_cores <- 8L
    cat(sprintf("  Rasterizing %d species on %d cores ...\n",
                length(todo_ids), n_cores))

    results <- parallel::mclapply(seq_along(todo_ids), function(i) {
      sid <- todo_ids[i]
      out_file <- file.path(out_dir, paste0(prefix, sid, ".tif"))
      sp_polys <- sp_vect[sp_vect$species_id == sid, ]
      sp_rast <- rasterize(sp_polys, grid_template, field = 1, background = 0)
      writeRaster(sp_rast, out_file, datatype = "INT1U", overwrite = TRUE,
                  gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
      if (i %% 50 == 0 || i == length(todo_ids))
        cat(sprintf("  [%d/%d] %s id=%s (%s)\n",
                    i, length(todo_ids), taxon_name,
                    sid, sp_polys$sci_name[1]))
      return(TRUE)
    }, mc.cores = n_cores, mc.preschedule = FALSE)

    n_errors <- sum(sapply(results, inherits, "try-error"))
    cat(sprintf("  Done. %d species rasterized.%s\n",
                length(todo_ids) - n_errors,
                if (n_errors > 0) sprintf(" %d errors.", n_errors) else ""))
  }

  rm(sp_vect)
  gc()

} # end taxon loop


############################
### summary ################
############################
{
cat("\n=== Final summary ===\n")
all_rasters <- list.files(out_dir, pattern = "^sp_.*\\.tif$")
for (tx in names(taxa)) {
  pat <- paste0("^sp_", tx, "_\\d+\\.tif$")
  n <- sum(grepl(pat, all_rasters))
  cat(sprintf("  %-12s %d rasters\n", tx, n))
}
cat(sprintf("  %-12s %d rasters\n", "TOTAL", length(all_rasters)))

# Check for zero-byte files
sizes <- file.size(file.path(out_dir, all_rasters))
n_zero <- sum(sizes == 0, na.rm = TRUE)
if (n_zero > 0) {
  cat(sprintf("  WARNING: %d zero-byte raster files\n", n_zero))
} else {
  cat("  OK: no zero-byte files\n")
}
}
