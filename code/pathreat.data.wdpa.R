##########################################
### pathreat.data.wdpa.R ################
### Rasterize WDPA polygons to pixel  ####
### level PA dataset (row-aligned     ####
### with prepdata.fst)                ####
##########################################

library(sf)
library(terra)
library(fst)
library(dplyr)
library(conflicted)
  conflict_prefer("select", "dplyr")
  conflict_prefer("filter", "dplyr")
  conflict_prefer("lag", "dplyr")
  
# --- Paths ---
union_poly_path <- "data/store/pathreat.data.wdpa.union.poly.gpkg"
grid_path       <- "data/store/pathreat.data.rasterbase.tif"
prep_path       <- "data/store/pathreat.data.prepdata.fst"
out_path        <- "data/store/pathreat.data.wdpa.fst"
raster_path     <- "data/store/pathreat.data.wdpa.raster.tif"

if (file.exists(raster_path)) {
  # ============================================================
  # Steps 1-4 (fast path): raster exists — skip rasterization
  # ============================================================
  cat("Raster exists — skipping rasterization, rebuilding pa_attrs only...\n")
  wdpa_raw <- st_read(union_poly_path, quiet = TRUE)
  cat("  Polygons from union gpkg:", nrow(wdpa_raw), "\n")

  pa_attrs <- st_drop_geometry(wdpa_raw) %>%
    select(
      wdpaid             = SITE_ID,
      pa_iucn_cat        = IUCN_CAT,
      pa_year_designated = STATUS_YR,
      gis_area           = GIS_AREA,
      iso3               = ISO3,
      desig_type         = DESIG_TYPE,
      desig_eng          = DESIG_ENG,
      gov_type           = GOV_TYPE,
      own_type           = OWN_TYPE
    ) %>%
    mutate(
      wdpaid             = as.integer(wdpaid),
      pa_year_designated = as.integer(pa_year_designated)
    ) %>%
    distinct(wdpaid, .keep_all = TRUE)
  cat("  Unique PAs:", nrow(pa_attrs), "\n")
  rm(wdpa_raw)
  gc()
  wdpa_rast <- rast(raster_path)

} else {
  # ============================================================
  # Step 1: Read & filter WDPA polygons
  # ============================================================
  cat("Reading union WDPA polygons...\n")
  wdpa <- st_read(union_poly_path, quiet = TRUE)
  cat("  Polygons from union gpkg:", nrow(wdpa), "\n")

  # ============================================================
  # Step 2: Build attribute lookup table
  # ============================================================
  cat("Building attribute lookup...\n")
  pa_attrs <- st_drop_geometry(wdpa) %>%
    select(
      wdpaid             = SITE_ID,
      pa_iucn_cat        = IUCN_CAT,
      pa_year_designated = STATUS_YR,
      gis_area           = GIS_AREA,
      iso3               = ISO3,
      desig_type         = DESIG_TYPE,
      desig_eng          = DESIG_ENG,
      gov_type           = GOV_TYPE,
      own_type           = OWN_TYPE
    ) %>%
    mutate(
      wdpaid             = as.integer(wdpaid),
      pa_year_designated = as.integer(pa_year_designated)
    ) %>%
    distinct(wdpaid, .keep_all = TRUE)

  cat("  Unique PAs:", nrow(pa_attrs), "\n")
  cat("  pa_year_designated == 0 (unknown):", sum(pa_attrs$pa_year_designated == 0, na.rm = TRUE), "\n")

  # ============================================================
  # Step 3: Sort for overlap priority (oldest PA wins)
  # ============================================================
  # terra::rasterize() with fun="last" keeps the last polygon written.
  # Sort newest first (+ STATUS_YR==0 as newest) so older PAs overwrite newer ones.
  cat("Sorting polygons (newest first, oldest last)...\n")
  wdpa <- wdpa %>%
    mutate(sort_yr = ifelse(STATUS_YR == 0, -Inf, STATUS_YR)) %>%
    arrange(desc(sort_yr)) %>%
    select(-sort_yr)

  # ============================================================
  # Step 4: Reproject to Mollweide & rasterize
  # ============================================================
  cat("Loading grid template...\n")
  grid_template <- rast(grid_path)
  mollweide_crs <- crs(grid_template)

  cat("Converting to terra SpatVector...\n")
  wdpa_vect <- vect(wdpa)

  cat("Reprojecting to Mollweide...\n")
  wdpa_vect <- project(wdpa_vect, mollweide_crs)

  # Free the sf object
  rm(wdpa)
  gc()

  cat("Rasterizing SITE_ID to 1-km grid (this may take 30-90 min)...\n")
  t0 <- proc.time()
  wdpa_rast <- rasterize(wdpa_vect, grid_template, field = "SITE_ID")
  cat("  Rasterization took", round((proc.time() - t0)[3] / 60, 1), "min\n")

  # Free the SpatVector
  writeRaster(wdpa_rast, raster_path)
  rm(wdpa_vect)
  gc()
}

# ============================================================
# Step 5: Extract at prepdata coordinates (row-aligned)
# ============================================================
cat("Reading prepdata coordinates...\n")
coords <- read_fst(prep_path, columns = c("pixel_id", "x", "y"))
n_pixels <- nrow(coords)
cat("  Pixels:", n_pixels, "\n")

cat("Extracting wdpaid at pixel locations...\n")
t0 <- proc.time()
wdpaid_values <- terra::extract(wdpa_rast, cbind(coords$x, coords$y))[, 1]
length(unique(wdpaid_values))

cat("  Extraction took", round((proc.time() - t0)[3] / 60, 1), "min\n")

pixel_ids <- coords$pixel_id
rm(coords, wdpa_rast)
gc()

# ============================================================
# Step 6: Join with attribute lookup & derive variables
# ============================================================
cat("Joining attributes...\n")
d <- data.frame(pixel_id = pixel_ids, wdpaid = as.integer(wdpaid_values))
rm(pixel_ids)
rm(wdpaid_values); gc()

d <- left_join(d, pa_attrs, by = "wdpaid")

# Derive iucn_class (matching existing convention in prepdata)
d$iucn_class <- case_when(
  d$pa_iucn_cat %in% c("Ia", "Ib", "II", "III", "IV") ~ "Strict protection",
  d$pa_iucn_cat %in% c("V", "VI")                      ~ "Less strict",
  !is.na(d$pa_iucn_cat)                                 ~ "Not Reported",
  TRUE                                               ~ NA_character_
)

# Derive size_class (matching existing bins: <10, 10-500, 500-2000, >=2000)
d$size_class <- case_when(
  d$gis_area <   10                      ~ "class 1",
  d$gis_area >=  10 & d$gis_area <  500  ~ "class 2",
  d$gis_area >= 500 & d$gis_area < 2000  ~ "class 3",
  d$gis_area >= 2000                     ~ "class 4",
  TRUE                                   ~ NA_character_
)

# ============================================================
# Step 7: Verify & save
# ============================================================
cat("\n--- Verification ---\n")
cat("Output rows:", nrow(d), "\n")
cat("Prepdata rows:", n_pixels, "\n")
stopifnot(nrow(d) == n_pixels)

cat("\nwdpaid coverage:\n")
cat("  Non-NA:", sum(!is.na(d$wdpaid)), sprintf("(%.1f%%)\n", 100 * mean(!is.na(d$wdpaid))))
cat("  NA:    ", sum(is.na(d$wdpaid)),  sprintf("(%.1f%%)\n", 100 * mean(is.na(d$wdpaid))))

cat("\niucn_class:\n")
print(table(d$iucn_class, useNA = "always"))

cat("\nsize_class:\n")
print(table(d$size_class, useNA = "always"))

cat("\nWriting to", out_path, "...\n")
write_fst(d, out_path)
cat("Done.\n")

# ============================================================
# Step 8: Summary of output
# ============================================================
cat("\n--- Output summary ---\n")
cat("Columns:", paste(names(d), collapse = ", "), "\n")
cat("  PA pixels (wdpaid non-NA):", sum(!is.na(d$wdpaid)), "\n")
cat("  Control pixels (wdpaid NA):", sum(is.na(d$wdpaid)), "\n")
