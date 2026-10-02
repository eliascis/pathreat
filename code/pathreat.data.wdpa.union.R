##########################################
### pathreat.data.wdpa.union.R ##########
### Merge public + restricted WDPA    ###
### layers into unified GeoPackages   ###
### (polygons + points)               ###
##########################################

library(sf)
library(dplyr)
library(conflicted)
  conflict_prefer("select", "dplyr")
  conflict_prefer("filter", "dplyr")

# --- Paths ---
gdb_path            <- "data/external/WDPA_Feb2026_Public/WDPA_Feb2026_Public.gdb"
restricted_poly_path <- "data/external/Restricted_PAs/WDPA_poly_restricted.shp"
restricted_pts_path  <- "data/external/Restricted_PAs/WDPA_points_restricted.shp"
out_poly_path        <- "data/store/pathreat.data.wdpa.union.poly.gpkg"
out_pts_path         <- "data/store/pathreat.data.wdpa.union.points.gpkg"

############################
### polygon layer ##########
############################
{
  cat("=== Polygon layer ===\n")
  cat("Reading WDPA polygons from GDB...\n")
  wdpa_poly <- st_read(gdb_path, layer = "WDPA_poly_Feb2026") %>%
    filter(REALM == "Terrestrial", STATUS == "Designated")
  cat("  Public polygons (Terrestrial & Designated):", nrow(wdpa_poly), "\n")

  # Remove CHN/IND from public layer
  n_removed <- sum(wdpa_poly$ISO3 %in% c("CHN", "IND"))
  wdpa_poly <- wdpa_poly %>% filter(!ISO3 %in% c("CHN", "IND"))
  cat("  CHN/IND polygons removed:", n_removed, "\n")

  # Read restricted replacements
  restricted_poly <- st_read(restricted_poly_path, quiet = TRUE) %>%
    filter(
      ISO3 %in% c("CHN", "IND"),
      MARINE %in% c("0", "2"),
      STATUS == "Designated"
    ) %>%
    rename(SITE_ID = WDPAID)
  cat("  CHN/IND polygons from restricted:", nrow(restricted_poly), "\n")

  # Union and select needed columns
  wdpa_poly <- bind_rows(wdpa_poly, restricted_poly) %>%
    select(
      SITE_ID,
      IUCN_CAT,
      STATUS_YR,
      GIS_AREA,
      ISO3,
      DESIG_TYPE,
      DESIG_ENG,
      GOV_TYPE,
      OWN_TYPE
    )
  rm(restricted_poly)
  cat("  Total polygons after union:", nrow(wdpa_poly), "\n")

  cat("Writing", out_poly_path, "...\n")
  st_write(wdpa_poly, out_poly_path, delete_dsn = TRUE, quiet = TRUE)
  cat("  Done.\n")
  rm(wdpa_poly)
  gc()
}

############################
### point layer ############
############################
{
  cat("\n=== Point layer ===\n")
  cat("Reading WDPA points from GDB...\n")
  wdpa_pts <- st_read(gdb_path, layer = "WDPA_point_Feb2026") %>%
    filter(REALM == "Terrestrial", STATUS == "Designated")
  cat("  Public points (Terrestrial & Designated):", nrow(wdpa_pts), "\n")

  # Remove CHN/IND from public layer
  n_removed <- sum(wdpa_pts$ISO3 %in% c("CHN", "IND"))
  wdpa_pts <- wdpa_pts %>% filter(!ISO3 %in% c("CHN", "IND"))
  cat("  CHN/IND points removed:", n_removed, "\n")

  # Read restricted replacements
  restricted_pts <- st_read(restricted_pts_path, quiet = TRUE) %>%
    filter(
      ISO3 %in% c("CHN", "IND"),
      MARINE %in% c("0", "2"),
      STATUS == "Designated"
    ) %>%
    rename(SITE_ID = WDPAID)
  cat("  CHN/IND points from restricted:", nrow(restricted_pts), "\n")

  # Union and select needed columns
  wdpa_pts <- bind_rows(wdpa_pts, restricted_pts) %>%
    select(
      SITE_ID,
      STATUS_YR,
      REP_AREA,
      ISO3
    )
  rm(restricted_pts)
  cat("  Total points after union:", nrow(wdpa_pts), "\n")

  cat("Writing", out_pts_path, "...\n")
  st_write(wdpa_pts, out_pts_path, delete_dsn = TRUE, quiet = TRUE)
  cat("  Done.\n")
  rm(wdpa_pts)
  gc()
}

cat("\n=== Union complete ===\n")
cat("  Polygons:", out_poly_path, "\n")
cat("  Points:  ", out_pts_path, "\n")
