##########################################
### pathreat.data.wdpa.points.R ##########
### Buffer WDPA point PAs and rasterize ##
### to pixel-level flags (row-aligned  ###
### with prepdata.fst)                 ###
##########################################

library(sf)
library(terra)
library(fst)
library(dplyr)
library(conflicted)
  conflict_prefer("select", "dplyr")
  conflict_prefer("filter", "dplyr")
  conflict_prefer("lag", "dplyr")

############################
### paths ##################
############################
union_pts_path <- "data/store/pathreat.data.wdpa.union.points.gpkg"
grid_path      <- "data/store/pathreat.data.rasterbase.tif"
prep_path      <- "data/store/pathreat.data.prepdata.fst"
out_path       <- "data/store/pathreat.data.wdpa.points.fst"
raster_path    <- "data/store/pathreat.data.wdpa.points.raster.tif"

############################
### rasterize (or skip) ####
############################
if (file.exists(raster_path)) {
  # ============================================================
  # Fast path: raster exists — skip rasterization
  # ============================================================
  cat("Points raster exists — skipping rasterization, rebuilding pt_attrs only...\n")
  pts_raw <- st_read(union_pts_path, quiet = TRUE)
  cat("  Points from union gpkg:", nrow(pts_raw), "\n")

  pt_attrs <- st_drop_geometry(pts_raw) %>%
    select(
      SITE_ID   = SITE_ID,
      STATUS_YR = STATUS_YR,
      REP_AREA  = REP_AREA
    ) %>%
    mutate(
      SITE_ID   = as.integer(SITE_ID),
      STATUS_YR = as.integer(STATUS_YR)
    ) %>%
    distinct(SITE_ID, .keep_all = TRUE)
  cat("  Unique point PAs:", nrow(pt_attrs), "\n")
  rm(pts_raw)
  gc()
  pts_rast <- rast(raster_path)

} else {
  # ============================================================
  # Step 1: Read & filter WDPA point layer
  # ============================================================
  cat("Reading union WDPA points...\n")
  pts <- st_read(union_pts_path, quiet = TRUE)
  cat("  Points from union gpkg:", nrow(pts), "\n")

  # ============================================================
  # Step 2: Build attribute lookup table
  # ============================================================
  cat("Building attribute lookup...\n")
  pt_attrs <- st_drop_geometry(pts) %>%
    select(
      SITE_ID   = SITE_ID,
      STATUS_YR = STATUS_YR,
      REP_AREA  = REP_AREA
    ) %>%
    mutate(
      SITE_ID   = as.integer(SITE_ID),
      STATUS_YR = as.integer(STATUS_YR)
    ) %>%
    distinct(SITE_ID, .keep_all = TRUE)

  cat("  Unique point PAs:", nrow(pt_attrs), "\n")
  cat("  STATUS_YR == 0 (unknown):", sum(pt_attrs$STATUS_YR == 0, na.rm = TRUE), "\n")
  cat("  REP_AREA range: [",
      min(pts$REP_AREA, na.rm = TRUE), ",",
      max(pts$REP_AREA, na.rm = TRUE), "] km²\n")

  # ============================================================
  # Step 3: Compute per-point buffer radius from REP_AREA
  # ============================================================
  # REP_AREA is in km²; convert to metres for st_buffer
  # pmax guards zero/negative values (1e-3 km² → ~18 m radius)
  cat("Computing buffer radii from REP_AREA...\n")
  pts <- pts %>%
    mutate(radius_m = sqrt(pmax(REP_AREA, 1e-3) * 1e6 / pi))
  cat("  Radius range (m): [",
      round(min(pts$radius_m), 1), ",",
      round(max(pts$radius_m), 1), "]\n")

  # ============================================================
  # Step 4: Sort for overlap priority (oldest PA wins)
  # ============================================================
  # terra::rasterize() with fun="last" keeps the last feature written.
  # Sort newest first (STATUS_YR==0 treated as newest → -Inf) so older PAs overwrite.
  cat("Sorting points (newest first, oldest last)...\n")
  pts <- pts %>%
    mutate(sort_yr = ifelse(STATUS_YR == 0, -Inf, STATUS_YR)) %>%
    arrange(desc(sort_yr)) %>%
    select(-sort_yr)

  # ============================================================
  # Step 5: Buffer in geographic CRS, then reproject & rasterize
  # ============================================================
  cat("Loading grid template...\n")
  grid_template <- rast(grid_path)
  mollweide_crs <- crs(grid_template)

  # Buffer in WGS84 (st_buffer with per-feature dist vector)
  cat("Buffering point PAs (per-feature radius)...\n")
  pts_buf <- st_buffer(pts, dist = pts$radius_m)
  rm(pts)
  gc()

  cat("Converting to terra SpatVector...\n")
  pts_vect <- vect(pts_buf)
  rm(pts_buf)
  gc()

  cat("Reprojecting to Mollweide...\n")
  pts_vect <- project(pts_vect, mollweide_crs)

  cat("Rasterizing SITE_ID to 1-km grid...\n")
  t0 <- proc.time()
  pts_rast <- rasterize(pts_vect, grid_template, field = "SITE_ID")
  cat("  Rasterization took", round((proc.time() - t0)[3] / 60, 1), "min\n")

  writeRaster(pts_rast, raster_path, overwrite = TRUE)
  cat("  Points raster written to", raster_path, "\n")
  rm(pts_vect)
  gc()
}

############################
### extract at pixels ######
############################
cat("\nReading prepdata coordinates...\n")
coords <- read_fst(prep_path, columns = c("pixel_id", "x", "y"))
n_pixels <- nrow(coords)
cat("  Pixels:", format(n_pixels, big.mark = ","), "\n")

cat("Extracting point-PA SITE_ID at pixel locations...\n")
t0 <- proc.time()
site_id_vals <- terra::extract(pts_rast, cbind(coords$x, coords$y))[, 1]
cat("  Extraction took", round((proc.time() - t0)[3] / 60, 1), "min\n")
cat("  Covered pixels (non-NA):", sum(!is.na(site_id_vals)), "\n")

pixel_ids <- coords$pixel_id
rm(coords, pts_rast)
gc()

############################
### join attrs & derive ####
############################
cat("Joining STATUS_YR from point PA attributes...\n")
d_pts <- data.frame(
  pixel_id = pixel_ids,
  SITE_ID  = as.integer(site_id_vals)
)
rm(pixel_ids, site_id_vals)
gc()

d_pts <- left_join(
  d_pts,
  pt_attrs %>% select(SITE_ID, STATUS_YR),
  by = "SITE_ID"
)

d_pts <- d_pts %>%
  mutate(
    pa_from_point                  = !is.na(SITE_ID),
    pa_from_point_designation_year = if_else(pa_from_point, as.integer(STATUS_YR), NA_integer_)
  ) %>%
  select(pixel_id, pa_from_point, pa_from_point_designation_year)

############################
### verify & save ##########
############################
cat("\n--- Verification ---\n")
cat("Output rows:", format(nrow(d_pts), big.mark = ","), "\n")
cat("Prepdata rows:", format(n_pixels, big.mark = ","), "\n")
stopifnot(nrow(d_pts) == n_pixels)

cat("\npa_from_point:\n")
print(table(d_pts$pa_from_point, useNA = "always"))
cat("\npa_from_point_designation_year (non-NA):\n")
print(summary(d_pts$pa_from_point_designation_year[d_pts$pa_from_point]))

cat("\nWriting to", out_path, "...\n")
write_fst(d_pts, out_path)
cat("Done.\n")

############################
### summary ################
############################
cat("\n--- Output summary ---\n")
cat("Columns:", paste(names(d_pts), collapse = ", "), "\n")
cat("  pa_from_point TRUE:", sum(d_pts$pa_from_point, na.rm = TRUE), "\n")
cat("  pa_from_point FALSE:", sum(!d_pts$pa_from_point, na.rm = TRUE), "\n")
