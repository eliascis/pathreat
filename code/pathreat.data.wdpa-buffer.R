##########################################
### pathreat.data.wdpa-buffer.R #########
### 5 km buffer zones around PAs     ####
### Uses wdpa.fst + prepdata.fst     ####
### (no WDPA geodatabase needed)     ####
##########################################

library(terra)
library(dplyr)
library(fst)

# --- Paths ---
grid_path  <- "data/store/pathreat.data.rasterbase.tif"
prep_path  <- "data/store/pathreat.data.prepdata.fst"
wdpa_path  <- "data/store/pathreat.data.wdpa.fst"
out_path   <- "data/store/pathreat.data.wdpa-buffer.fst"

buffer_dist <- 5000  # meters
tmp_dir     <- "data/tmp"

# ============================================================
# Step 1: Load grid template for CRS and dimensions
# ============================================================
cat("Loading grid template...\n")
grid_template <- rast(grid_path)
mollweide_crs <- crs(grid_template)

# ============================================================
# Step 2: Read prepdata (x, y) + wdpa (wdpaid, pa_year_designated)
# ============================================================
cat("Reading prepdata coordinates...\n")
coords <- read_fst(prep_path, columns = c("pixel_id", "x", "y"))
n_pixels <- nrow(coords)
cat("  Pixels:", format(n_pixels, big.mark = ","), "\n")

cat("Reading WDPA data...\n")
d_wdpa <- read_fst(wdpa_path, columns = c("wdpaid", "pa_year_designated"))
stopifnot(nrow(d_wdpa) == n_pixels)

# ============================================================
# Step 3: Filter to PA pixels with pa_year_designated <= 2020
#          (includes 0 = unknown year — conservative, prevents
#           near-PA contamination of controls)
# ============================================================
cat("Filtering to study-period PA pixels...\n")
is_pa <- !is.na(d_wdpa$wdpaid) & d_wdpa$pa_year_designated <= 2020
cat("  PA pixels (year <= 2020 or unknown):", format(sum(is_pa), big.mark = ","),
    sprintf("(%.1f%%)\n", 100 * mean(is_pa)))

pa_coords <- cbind(coords$x[is_pa], coords$y[is_pa])
rm(d_wdpa); gc()

# ============================================================
# Step 4: Create SpatVector from PA pixel coords → rasterize
# ============================================================
cat("Creating PA SpatVector...\n")
pa_points <- vect(pa_coords, type = "points", crs = mollweide_crs)
rm(pa_coords); gc()

cat("Rasterizing PA pixels to binary mask...\n")
t0 <- proc.time()
pa_mask <- rasterize(pa_points, grid_template, field = 1)
cat("  Rasterization took", round((proc.time() - t0)[3] / 60, 1), "min\n")
rm(pa_points); gc()

# ============================================================
# Step 5: Compute distance from each pixel to nearest PA pixel
# ============================================================
# Write to disk with BIGTIFF — distance raster is Float64 on a 20k×38k grid
# (~5.7 GB), exceeding the 4 GB standard TIFF limit
dist_file <- file.path(tmp_dir, "pa_distance.tif")
cat("Computing distance to nearest PA pixel (this may take a while)...\n")
cat("  Writing distance raster to", dist_file, "\n")
t0 <- proc.time()
dist_rast <- distance(pa_mask, filename = dist_file, overwrite = TRUE,
                      gdal = c("BIGTIFF=YES"))
cat("  Distance computation took", round((proc.time() - t0)[3] / 60, 1), "min\n")
rm(pa_mask); gc()

# ============================================================
# Step 6: Buffer raster: distance > 0 & distance <= 5000 m
# ============================================================
cat("Creating buffer indicator (distance <= ", buffer_dist, " m)...\n")
buffer_rast <- (dist_rast > 0 & dist_rast <= buffer_dist)
rm(dist_rast); gc()
# clean up temp distance raster
unlink(dist_file)

# ============================================================
# Step 7: Extract buffer_5km at all prepdata (x, y) coordinates
# ============================================================
cat("Extracting buffer values at pixel locations...\n")
t0 <- proc.time()
buffer_values <- extract(buffer_rast, cbind(coords$x, coords$y))[, 1]
cat("  Extraction took", round((proc.time() - t0)[3] / 60, 1), "min\n")
pixel_ids <- coords$pixel_id
rm(buffer_rast, coords); gc()

# Convert to integer (0/1), NAs → 0
buffer_values <- as.integer(ifelse(is.na(buffer_values), 0L, buffer_values))

# ============================================================
# Step 8: Verify & save
# ============================================================
cat("\n--- Verification ---\n")
cat("Output rows:", format(length(buffer_values), big.mark = ","), "\n")
cat("Prepdata rows:", format(n_pixels, big.mark = ","), "\n")
stopifnot(length(buffer_values) == n_pixels)

cat("\nbuffer_5km distribution:\n")
print(table(buffer_5km = buffer_values, useNA = "always"))
cat(sprintf("  Buffer pixels: %s (%.1f%% of all pixels)\n",
            format(sum(buffer_values == 1), big.mark = ","),
            100 * mean(buffer_values == 1)))

# Save as single-column fst
d_out <- data.frame(pixel_id = pixel_ids, buffer_5km = buffer_values)
cat("\nWriting to", out_path, "...\n")
write_fst(d_out, out_path)
cat("Done.\n")
