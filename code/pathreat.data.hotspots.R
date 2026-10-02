##########################################
### pathreat.data.hotspots.R #############
### Rasterize biodiversity hotspots   ####
### to pixel-level dataset            ####
### (row-aligned with prepdata.fst)   ####
##########################################

library(sf)
library(terra)
library(fst)
library(dplyr)
library(conflicted)
  conflict_prefer("select", "dplyr")
  conflict_prefer("filter", "dplyr")

# --- Paths ---
gpkg_path   <- "data/external/hotspots_wgs84.gpkg"
grid_path   <- "data/store/pathreat.data.rasterbase.tif"
prep_path   <- "data/store/pathreat.data.prepdata.fst"
out_path    <- "data/store/pathreat.data.hotspots.fst"
raster_path <- "data/store/pathreat.data.hotspots.raster.tif"

# ============================================================
# Step 1: Read hotspot polygons
# ============================================================
cat("Reading hotspot polygons...\n")
hotspots <- st_read(gpkg_path, layer = "hotspots")
cat("  Features:", nrow(hotspots), "\n")

# ============================================================
# Step 2: Create hotspot_id and lookup table
# ============================================================
cat("Building hotspot lookup...\n")
hotspots$hotspot_id <- seq_len(nrow(hotspots))

hotspot_lookup <- st_drop_geometry(hotspots) %>%
  select(
    hotspot_id,
    hotspot_name = NAME
  )

cat("  Hotspots:", nrow(hotspot_lookup), "\n")
cat("  Names:\n")
print(hotspot_lookup)

# ============================================================
# Step 3: Reproject to Mollweide & rasterize
# ============================================================
if (file.exists(raster_path)) {
  cat("Raster exists — skipping rasterization...\n")
  hotspot_rast <- rast(raster_path)

} else {
  cat("Loading grid template...\n")
  grid_template <- rast(grid_path)
  mollweide_crs <- crs(grid_template)

  cat("Converting to SpatVector and reprojecting to Mollweide...\n")
  hotspots_vect <- vect(hotspots)
  hotspots_vect <- project(hotspots_vect, mollweide_crs)

  rm(hotspots)
  gc()

  cat("Rasterizing hotspot_id to 1-km grid...\n")
  t0 <- proc.time()
  hotspot_rast <- rasterize(hotspots_vect, grid_template, field = "hotspot_id")
  cat("  Rasterization took", round((proc.time() - t0)[3] / 60, 1), "min\n")

  writeRaster(hotspot_rast, raster_path)
  rm(hotspots_vect)
  gc()
}

# ============================================================
# Step 4: Extract at prepdata coordinates (row-aligned)
# ============================================================
cat("Reading prepdata coordinates...\n")
coords <- read_fst(prep_path, columns = c("pixel_id", "x", "y"))
n_pixels <- nrow(coords)
cat("  Pixels:", n_pixels, "\n")

cat("Extracting hotspot_id at pixel locations...\n")
t0 <- proc.time()
hotspot_values <- terra::extract(hotspot_rast, cbind(coords$x, coords$y))[, 1]
cat("  Extraction took", round((proc.time() - t0)[3] / 60, 1), "min\n")

pixel_ids <- coords$pixel_id
rm(coords, hotspot_rast)
gc()

# ============================================================
# Step 5: Build output and join hotspot names
# ============================================================
cat("Building output data.frame...\n")
d <- data.frame(
  pixel_id = pixel_ids,
  hotspot_id = as.integer(hotspot_values)
)
rm(pixel_ids, hotspot_values)
gc()

d <- left_join(d, hotspot_lookup, by = "hotspot_id")

# ============================================================
# Step 6: Verify & save
# ============================================================
cat("\n--- Verification ---\n")
cat("Output rows:", nrow(d), "\n")
cat("Prepdata rows:", n_pixels, "\n")
stopifnot(nrow(d) == n_pixels)

cat("\nHotspot coverage:\n")
cat("  In hotspot:", sum(!is.na(d$hotspot_id)),
    sprintf("(%.1f%%)\n", 100 * mean(!is.na(d$hotspot_id))))
cat("  Not in hotspot:", sum(is.na(d$hotspot_id)),
    sprintf("(%.1f%%)\n", 100 * mean(is.na(d$hotspot_id))))
cat("  Unique hotspots:", length(unique(na.omit(d$hotspot_id))), "\n")

cat("\nHotspot counts:\n")
print(table(d$hotspot_name, useNA = "always"))

cat("\nWriting to", out_path, "...\n")
write_fst(d, out_path)
cat("Done.\n")
