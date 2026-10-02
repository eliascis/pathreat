##########################################
### pathreat.data.admin1.R ##############
### Rasterize admin1 boundaries      ####
### from Natural Earth (global)      ####
### to pixel-level dataset           ####
### (row-aligned with prepdata.fst)  ####
##########################################

library(sf)
library(terra)
library(rnaturalearth)
library(fst)
library(dplyr)

# ne_states() at its default scale needs the rnaturalearthhires data package,
# which is distributed via rOpenSci's r-universe rather than CRAN.
if (!requireNamespace("rnaturalearthhires", quietly = TRUE)) {
  stop(
    "ne_states() requires 'rnaturalearthhires'. Install it with:\n",
    "  install.packages('rnaturalearthhires', ",
    "repos = c('https://ropensci.r-universe.dev', 'https://cloud.r-project.org'))"
  )
}

# --- Paths ---
grid_path   <- "data/store/pathreat.data.rasterbase.tif"
prep_path   <- "data/store/pathreat.data.prepdata.fst"
out_path    <- "data/store/pathreat.data.admin1.fst"
raster_path <- "data/store/pathreat.data.admin1.raster.tif"

############################
### load admin1 boundaries #
############################
{
cat("Loading Natural Earth admin1 boundaries...\n")
admin1 <- ne_states(returnclass = "sf")
cat("  Features:", nrow(admin1), "\n")
cat("  Countries:", length(unique(admin1$iso_a2)), "\n")

admin1$admin1_id <- seq_len(nrow(admin1))
}

############################
### rasterize ##############
############################
{
if (file.exists(raster_path)) {
  cat("Raster exists — skipping rasterization...\n")
  admin1_rast <- rast(raster_path)

} else {
  cat("Loading grid template...\n")
  grid_template <- rast(grid_path)
  mollweide_crs <- crs(grid_template)

  cat("Reprojecting to Mollweide...\n")
  admin1_vect <- vect(admin1)
  admin1_vect <- project(admin1_vect, mollweide_crs)

  rm(admin1)
  gc()

  cat("Rasterizing admin1_id to 1-km grid...\n")
  t0 <- proc.time()
  admin1_rast <- rasterize(admin1_vect, grid_template, field = "admin1_id")
  cat("  Rasterization took", round((proc.time() - t0)[3] / 60, 1), "min\n")

  writeRaster(
    admin1_rast,
    raster_path,
    datatype  = "INT2U",
    gdal      = c("COMPRESS=DEFLATE", "BIGTIFF=YES"),
    overwrite = TRUE
  )

  rm(admin1_vect, grid_template)
  gc()
}
}

############################
### extract at prepdata ####
############################
{
cat("Reading prepdata coordinates...\n")
coords <- read_fst(prep_path, columns = c("pixel_id", "x", "y"))
n_pixels <- nrow(coords)
cat("  Pixels:", format(n_pixels, big.mark = ","), "\n")

cat("Extracting admin1_id at pixel locations...\n")
t0 <- proc.time()
admin1_values <- terra::extract(admin1_rast, cbind(coords$x, coords$y))[, 1]
cat("  Extraction took", round((proc.time() - t0)[3] / 60, 1), "min\n")

pixel_ids <- coords$pixel_id
rm(coords, admin1_rast)
gc()
}

############################
### build & save ###########
############################
{
cat("Building output data.frame...\n")
d <- data.frame(
  pixel_id  = pixel_ids,
  admin1_id = as.integer(admin1_values)
)
rm(pixel_ids, admin1_values)
gc()

cat("\n--- Verification ---\n")
cat("Output rows:", format(nrow(d), big.mark = ","), "\n")
cat("Prepdata rows:", format(n_pixels, big.mark = ","), "\n")
stopifnot(nrow(d) == n_pixels)

cat("\nAdmin1 coverage:\n")
cat("  Assigned:", format(sum(!is.na(d$admin1_id)), big.mark = ","),
    sprintf("(%.1f%%)\n", 100 * mean(!is.na(d$admin1_id))))
cat("  NA:", format(sum(is.na(d$admin1_id)), big.mark = ","),
    sprintf("(%.1f%%)\n", 100 * mean(is.na(d$admin1_id))))
cat("  Unique admin1 units:", length(unique(na.omit(d$admin1_id))), "\n")

cat("\nWriting to", out_path, "...\n")
write_fst(d, out_path)
cat("Done.\n")
}


# d<-read_fst("data/store/pathreat.data.admin1.fst")
# head(d)
