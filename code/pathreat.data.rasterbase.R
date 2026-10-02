##########################################
### pathreat.data.rasterbase.R ##########
### Create lightweight shell raster   ####
### defining the project grid         ####
### (1-km Mollweide, 20k x 38k)      ####
##########################################

library(terra)
library(fst)

# --- Paths ---
out_path  <- "data/store/pathreat.data.rasterbase.tif"
prep_path <- "data/store/pathreat.data.prepdata.fst"

# --- Grid parameters (from PA_MAP_2022_b.tif header) ---
grid_nrows <- 20000L
grid_ncols <- 38000L
grid_xmin  <- -19000000
grid_xmax  <-  19000000
grid_ymin  <- -10000000
grid_ymax  <-  10000000
grid_crs   <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

############################
### create shell raster ####
############################
{
  cat("Creating shell raster...\n")
  cat("  Dimensions:", grid_nrows, "rows x", grid_ncols, "cols\n")
  cat("  Resolution: 1000 1000 (m)\n")
  cat("  Extent:    ", grid_xmin, grid_xmax, grid_ymin, grid_ymax, "\n")
  cat("  CRS:        World_Mollweide\n")

  shell <- rast(
    nrows = grid_nrows,
    ncols = grid_ncols,
    xmin  = grid_xmin,
    xmax  = grid_xmax,
    ymin  = grid_ymin,
    ymax  = grid_ymax,
    crs   = grid_crs
  )
  # row-based gradient (0-255) — compresses well, visible in QGIS
  values(shell) <- rep(seq(0L, 255L, length.out = grid_nrows), each = grid_ncols)

  cat("\nWriting shell raster to:", out_path, "\n")
  writeRaster(
    shell,
    out_path,
    datatype = "INT1U",
    gdal     = c("COMPRESS=DEFLATE"),
    overwrite = TRUE
  )
  cat("  File size:", round(file.size(out_path) / 1024, 1), "KB\n")
  rm(shell)
}

############################
### verify grid props ######
############################
{
  cat("\nVerification — reading back shell raster:\n")
  check <- rast(out_path)
  cat("  Dimensions:", nrow(check), "rows x", ncol(check), "cols\n")
  cat("  Resolution:", res(check), "\n")
  cat("  Extent:    ", as.vector(ext(check)), "\n")
  cat("  CRS:       ", crs(check, describe = TRUE)$name, "\n")
}

############################
### verify alignment #######
############################
{
  # prepdata (x, y) are upper-left cell corners, not centroids.
  # Centroids are offset by (+res/2, -res/2) from corners.
  # Verify: prepdata corner -> cell -> centroid should equal corner + (500, -500).
  cat("\nAlignment check against prepdata (upper-left corners)...\n")
  prep <- read_fst(prep_path, columns = c("x", "y"), as.data.table = FALSE)

  set.seed(42)
  idx <- sample(nrow(prep), min(100000, nrow(prep)))
  xy_corner <- as.matrix(prep[idx, c("x", "y")])
  rm(prep)

  # shift corners to centroids for cell lookup
  xy_centroid <- xy_corner
  xy_centroid[, 1] <- xy_corner[, 1] + 500
  xy_centroid[, 2] <- xy_corner[, 2] - 500

  # round-trip: centroid -> cell -> centroid
  cells <- cellFromXY(check, xy_centroid)
  xy_rt <- xyFromCell(check, cells)

  max_dev_x <- max(abs(xy_centroid[, 1] - xy_rt[, 1]))
  max_dev_y <- max(abs(xy_centroid[, 2] - xy_rt[, 2]))

  cat("  Sampled points:", length(idx), "\n")
  cat("  Max |x deviation|:", max_dev_x, "\n")
  cat("  Max |y deviation|:", max_dev_y, "\n")

  if (max_dev_x == 0 & max_dev_y == 0) {
    cat("  PASS: grid perfectly aligned with prepdata coordinates.\n")
  } else {
    stop("FAIL: alignment mismatch! Max deviation: x=",
         max_dev_x, " y=", max_dev_y)
  }
  rm(check, xy_corner, xy_centroid, xy_rt, cells)
}

cat("\nDone.\n")
