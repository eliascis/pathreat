############################
### countries ##############
############################
#
# 1. Downloads GADM level-0 boundaries for all countries → countries.gpkg (Mollweide)
# 2. Rasterizes onto 1-km grid, extracts at prepdata coordinates → countries.fst
#
# Unique identifier: GID_0 (GADM country code, unique for every territory)
#
# Output: data/store/pathreat.data.countries.gpkg
#         data/store/pathreat.data.countries.fst
#
# Dependencies: sf, terra, fst, geodata, countrycode, dplyr

library(sf)
library(terra)
library(fst)
library(geodata)
library(countrycode)
library(dplyr)
nowrun=0

grid_template <- rast("data/store/pathreat.data.rasterbase.tif")
mollweide_crs <- crs(grid_template)

gadm_dir <- "data/external/gadm"
out_path <- "data/store/pathreat.data.countries.gpkg"
dir.create(gadm_dir, recursive = TRUE, showWarnings = FALSE)

####################################
### download GADM level-0 #########
####################################

if (file.exists(out_path)) {
  cat("Loading existing", out_path, "...\n")
  countries <- st_read(out_path, quiet = TRUE)
} else {
  iso3_codes <- countrycode::codelist$iso3c
  iso3_codes <- iso3_codes[!is.na(iso3_codes)]
  cat("Downloading GADM level-0 for", length(iso3_codes), "countries...\n")

  poly_list <- list()
  for (i in seq_along(iso3_codes)) {
    iso <- iso3_codes[i]
    adm <- tryCatch(
      gadm(country = iso, level = 0, path = gadm_dir),
      error = function(e) NULL
    )
    if (!is.null(adm)) {
      poly_list[[iso]] <- st_as_sf(adm)
    }
    if (i %% 50 == 0) {
      cat(sprintf("  %d/%d downloaded\n", i, length(iso3_codes)))
    }
  }
  cat(sprintf("  Downloaded %d countries\n", length(poly_list)))

  countries <- bind_rows(poly_list)
  rm(poly_list)

  # Add region attributes from countrycode
  countries <- countries %>%
    mutate(
      iso_a2 = countrycode(GID_0, "iso3c", "iso2c", warn = FALSE),
      name = COUNTRY,
      continent = countrycode(GID_0, "iso3c", "continent", warn = FALSE),
      region_un = countrycode(GID_0, "iso3c", "un.region.name", warn = FALSE),
      subregion = countrycode(GID_0, "iso3c", "un.regionsub.name", warn = FALSE)
    ) %>%
    rename(gid_0 = GID_0, country = COUNTRY) %>%
    select(gid_0, country, iso_a2, name, continent, region_un, subregion, geometry)

  countries <- st_transform(countries, mollweide_crs)
  countries$country_rast <- as.integer(factor(countries$country))
  countries <- countries[, c("gid_0", "country", "country_rast", "iso_a2", "name", "continent", "region_un", "subregion", "geom")]

  cat("Rows:", nrow(countries), "\n")
  cat("Unique GID_0:", length(unique(countries$gid_0)), "\n")
  cat("CRS:", st_crs(countries)$input, "\n")

  st_write(countries, out_path, delete_dsn = TRUE)
  cat("Saved to", out_path, "\n")
}

####################################
### rasterize countries ############
####################################
{
rast_path <- "data/store/pathreat.data.countries.raster.tif"
cat("Rasterizing countries onto grid template...\n")
t0 <- proc.time()
cty_vect <- vect(countries)
cty_rast <- rasterize(cty_vect, grid_template, field = "country_rast")
cat("  Rasterization took", round((proc.time() - t0)[3] / 60, 1), "min\n")

writeRaster(cty_rast, rast_path,datatype="INT2U",overwrite =TRUE,gdal = "BIGTIFF=YES")
cat("Saved to", rast_path, "\n")

rm(cty_vect, cty_rast)
gc()
}

####################################
### extract at prepdata coords #####
####################################
{
prep_path <- "data/store/pathreat.data.prepdata.fst"
fst_path  <- "data/store/pathreat.data.countries.fst"

cat("Reading prepdata coordinates...\n")
coords <- read_fst(prep_path, columns = c("pixel_id", "x", "y"))
cat("  Pixels:", nrow(coords), "\n")

cat("Extracting country_rast at pixel locations...\n")
cty_rast <- rast(rast_path)
t0 <- proc.time()
cty_values <- terra::extract(cty_rast, cbind(coords$x, coords$y))[, 1]
cat("  Extraction took", round((proc.time() - t0)[3] / 60, 1), "min\n")

rm(cty_rast)
gc()

# Build lookup from countries vector (one row per country_rast)
lookup <- st_drop_geometry(countries) %>%
  select(country_rast, gid_0, country) %>%
  distinct(country_rast, .keep_all = TRUE)

d <- data.frame(
  pixel_id = coords$pixel_id,
  country_rast = as.integer(cty_values)
)

rm(coords, cty_values)
gc()

n_before <- nrow(d)
d <- left_join(d, lookup, by = "country_rast")
stopifnot(nrow(d) == n_before)

cat("\n--- Verification ---\n")
cat("Rows:", nrow(d), "\n")
cat("Unique country_rast:", length(unique(na.omit(d$country_rast))), "\n")
cat("Unique countries:", length(unique(na.omit(d$country))), "\n")
cat("NA pixels:", sum(is.na(d$country_rast)),
    sprintf("(%.1f%%)\n", 100 * mean(is.na(d$country_rast))))

head(d)
nrow(d)
table(is.na(d$country))/nrow(d)
table(is.na(d$country_rast))/nrow(d)


cat("\nWriting to", fst_path, "...\n")
write_fst(d, fst_path)
cat("Done.\n")
}

# ####################################
# ### pixel-level data frame #########
# ####################################
# {
# # Build lookup: row index → country attributes
# countries$cty_idx <- seq_len(nrow(countries))
# lookup <- st_drop_geometry(countries)
# 
# # Rasterize country index onto 1-km grid
# cat("Rasterizing countries onto grid template...\n")
# t0 <- proc.time()
# cty_vect <- vect(countries)
# cty_rast <- rasterize(cty_vect, grid_template, field = "cty_idx")
# cat("  Rasterization took", round((proc.time() - t0)[3] / 60, 1), "min\n")
# 
# rm(cty_vect, countries)
# gc()
# 
# # Convert raster to data frame (non-NA cells only, with coordinates)
# cat("Converting raster to data frame...\n")
# t0 <- proc.time()
# d <- as.data.frame(cty_rast, xy = TRUE, na.rm = TRUE)
# names(d) <- c("x", "y", "cty_idx")
# d$cty_idx <- as.integer(d$cty_idx)
# d$pixel_id <- seq_len(nrow(d))
# cat("  Conversion took", round((proc.time() - t0)[3] / 60, 1), "min\n")
# cat("  Land pixels:", nrow(d), "\n")
# 
# rm(cty_rast)
# gc()
# 
# # Join country attributes
# d <- left_join(d, lookup, by = "cty_idx")
# d$cty_idx <- NULL
# }
# 
# ######################
# ### Verify & save ####
# ######################
# {
# cat("\n--- Verification ---\n")
# cat("Rows:", nrow(d), "\n")
# cat("Unique countries:", length(unique(d$country)), "\n")
# 
# 
# length(unique(d$gid_0))
# length(unique(d$country))
# length(unique(d$name))
# length(unique(d$country_rast))
# head(d)
# d$gid_0<-NULL
# d$iso_a2<-NULL
# d$name<-NULL
# d$continent<-NULL
# d$region_un<-NULL
# d$subregion<-NULL
# head(d)
# 
# fst_path <- "data/store/pathreat.data.countries.fst"
# cat("Writing to", fst_path, "...\n")
# write_fst(d, fst_path)
# cat("Done.\n")
# 
# }
# ############
# ### test ###
# ############
# if (nowrun==1) {
#   d <- st_read("data/store/pathreat.data.countries.gpkg") %>% data.frame()
#   head(d)
#   nrow(d)
#   length(unique(d$gid_0))
#   length(unique(d$country_rast))
# 
#   d <- read_fst("data/store/pathreat.data.countries.fst")
#   head(d)
#   length(unique(d$gid_0))
#   length(unique(d$country_rast))
# 
# }
