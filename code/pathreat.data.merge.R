library(fst)

source("code/pathreat.analysis.config.R")

# verify pixel_id alignment with random sample (fast check for row-aligned cbinds)
check_pixel_alignment <- function(ids_a, ids_b, n = 1000) {
  idx <- sample(length(ids_a), min(n, length(ids_a)))
  stopifnot(all(ids_a[idx] == ids_b[idx]))
}

#################
### read data ###
#################
{
##pixel pre-match data
d.um <- read_fst("data/store/pathreat.data.prepdata.fst")

}




################################
### merge deforestation data ###
################################
{
##reload main data
cat("Reloading main data for merge...\n")
d.def<-read_fst("data/store/pathreat.data.deforestation.fst")

##verify pixel_id alignment
cat("Main data rows:", nrow(d.um), "\n")
cat("Deforestation data rows:", nrow(d.def), "\n")
cat("Matching pixel_ids:", sum(d.um$pixel_id %in% d.def$pixel_id), "\n")

##merge by pixel_id
if (nrow(d.um)==nrow(d.def)){
  # d.um<-d.um[order(d.um$pixel_id),]
  # d.def<-d.def[order(d.def$pixel_id),]
  check_pixel_alignment(d.um$pixel_id, d.def$pixel_id)
  d.def$pixel_id <- NULL
  d<-cbind(d.um,d.def)
}
# d <- merge(d.um, d.def, by = "pixel_id", all.x = TRUE)
cat("Merged data rows:", nrow(d), "\n")

##clean up
rm(d.um, d.def)
gc()

head(d)
}




######################################
### merge land-cover vulnerability ###
######################################
{
cat("\n=== Merging land-cover vulnerability ===\n")
d_vuln <- read_fst("data/store/pathreat.data.land-vulnerability.fst", columns = c("pixel_id", "lu.vulnerable"))
cat("  Vulnerability rows:", nrow(d_vuln), "  Data rows:", nrow(d), "\n")
if (nrow(d) == nrow(d_vuln)) {
  check_pixel_alignment(d$pixel_id, d_vuln$pixel_id)
  d <- cbind(d, lu.vulnerable = d_vuln$lu.vulnerable)
  cat("  lu.vulnerable NAs:", sum(is.na(d$lu.vulnerable)), "\n")
} else {
  stop("Row count mismatch. Re-run land-vulnerability.R.")
}
rm(d_vuln)
gc()
}



########################
### merge WDPA data ####
########################
{
cat("\n=== Merging WDPA data ===\n")
d_wdpa <- read_fst("data/store/pathreat.data.wdpa.fst")
cat("  WDPA rows:", nrow(d_wdpa), "  Data rows:", nrow(d), "\n")

# row-aligned with prepdata (same row count & order before any filtering)
if (nrow(d) == nrow(d_wdpa)) {
  check_pixel_alignment(d$pixel_id, d_wdpa$pixel_id)
  d_wdpa$pixel_id <- NULL
  wdpa_vars <- names(d_wdpa)
  d <- cbind(d, d_wdpa)
  names(d)[names(d) == "gis_area"] <- "size"
  cat("  Columns added:", paste(wdpa_vars, collapse = ", "), "\n")
} else {
  stop("Row count mismatch. Re-run pathreat.data.wdpa.R.")
}
cat("  wdpaid NAs:", sum(is.na(d$wdpaid)), "\n")

rm(d_wdpa)
gc()
}


################################
### merge WDPA point flag ######
################################
{
cat("\n=== Merging WDPA point PA flags ===\n")
pts_flag_path <- "data/store/pathreat.data.wdpa.points.fst"
if (file.exists(pts_flag_path)) {
  d_wdpa_pts <- read_fst(pts_flag_path)
  cat("  WDPA points rows:", nrow(d_wdpa_pts), "  Data rows:", nrow(d), "\n")
  if (nrow(d) == nrow(d_wdpa_pts)) {
    check_pixel_alignment(d$pixel_id, d_wdpa_pts$pixel_id)
    d <- cbind(d,
      pa_from_point                  = d_wdpa_pts$pa_from_point,
      pa_from_point_designation_year = d_wdpa_pts$pa_from_point_designation_year
    )
    cat("  pa_from_point TRUE:", sum(d$pa_from_point, na.rm = TRUE), "\n")
  } else {
    stop("Row count mismatch. Re-run pathreat.data.wdpa.points.R.")
  }
  rm(d_wdpa_pts)
  gc()
} else {
  cat("  WDPA points file not found — skipping (run pathreat.data.wdpa.points.R first)\n")
}
}


################################
### merge WDPA buffer data #####
################################
{
cat("\n=== Merging WDPA buffer data ===\n")
d_buf <- read_fst("data/store/pathreat.data.wdpa-buffer.fst")
cat("  Buffer rows:", nrow(d_buf), "  Data rows:", nrow(d), "\n")

# row-aligned with prepdata (same row count & order before any filtering)
if (nrow(d) == nrow(d_buf)) {
  check_pixel_alignment(d$pixel_id, d_buf$pixel_id)
  d <- cbind(d, buffer_5km = d_buf$buffer_5km)
  cat("  buffer_5km distribution:\n")
  print(table(d$buffer_5km, useNA = "always"))
} else {
  stop("Row count mismatch. Re-run pathreat.data.wdpa-buffer.R.")
}

rm(d_buf)
gc()
}


################################
### merge country data #########
################################
{
cat("\n=== Merging country data ===\n")
cty_path <- "data/store/pathreat.data.countries.fst"
if (file.exists(cty_path)) {
  d_cty <- read_fst(cty_path)
  cat("  Country rows:", nrow(d_cty), "  Data rows:", nrow(d), "\n")
  if (nrow(d) == nrow(d_cty)) {
    check_pixel_alignment(d$pixel_id, d_cty$pixel_id)
    d_cty$pixel_id <- NULL
    # drop overlapping columns from d (prepdata may still carry them)
    overlap <- intersect(names(d), names(d_cty))
    if (length(overlap) > 0) {
      cat("  Replacing overlapping columns:", paste(overlap, collapse = ", "), "\n")
      d[overlap] <- NULL
    }
    d <- cbind(d, d_cty)
    cat("  Columns added:", paste(names(d_cty), collapse = ", "), "\n")
    cat("  country_rast NAs:", sum(is.na(d$country_rast)), "\n")
    cat("  gid_0 NAs:", sum(is.na(d$gid_0)), "\n")
  } else {
    stop("Row count mismatch. Re-run pathreat.data.countries.R.")
  }
  rm(d_cty)
  gc()
} else {
  cat("  Country file not found — skipping (run pathreat.data.countries.R first)\n")
}
}


################################
### merge admin1 data ##########
################################
{
cat("\n=== Merging admin1 data ===\n")
admin1_path <- "data/store/pathreat.data.admin1.fst"
if (file.exists(admin1_path)) {
  d_admin1 <- read_fst(admin1_path, columns = c("pixel_id", "admin1_id"))
  cat("  Admin1 rows:", nrow(d_admin1), "  Data rows:", nrow(d), "\n")
  if (nrow(d) == nrow(d_admin1)) {
    check_pixel_alignment(d$pixel_id, d_admin1$pixel_id)
    d <- cbind(d, admin1_id = d_admin1$admin1_id)
    cat("  admin1_id assigned:", sum(!is.na(d$admin1_id)),
        sprintf("(%.1f%%)\n", 100 * mean(!is.na(d$admin1_id))))
  } else {
    stop("Row count mismatch. Re-run pathreat.data.admin1.R.")
  }
  rm(d_admin1)
  gc()
} else {
  cat("  Admin1 file not found — skipping (run pathreat.data.admin1.R first)\n")
}
}




################################
### merge hotspots #############
################################
{
cat("\n=== Merging biodiversity hotspots ===\n")
hotspot_path <- "data/store/pathreat.data.hotspots.fst"
d_hs <- read_fst(hotspot_path)
cat("  Hotspot rows:", nrow(d_hs), "  Data rows:", nrow(d), "\n")
if (nrow(d) == nrow(d_hs)) {
  check_pixel_alignment(d$pixel_id, d_hs$pixel_id)
  d <- cbind(d,
    hotspot_id   = d_hs$hotspot_id,
    hotspot_name = d_hs$hotspot_name
  )
  cat("  In hotspot:", sum(!is.na(d$hotspot_id)),
      sprintf("(%.1f%%)\n", 100 * mean(!is.na(d$hotspot_id))))
} else {
  stop("Row count mismatch — skipping. Re-run pathreat.data.hotspots.R.")
}
rm(d_hs)
gc()
}


#############################
### remove country NAs #####
#############################
{
cat("\n=== Removing country NAs ===\n")
cat("  country NAs:", sum(is.na(d$country)), "\n")
cat("  country_rast NAs:", sum(is.na(d$country_rast)), "\n")
i <- which(is.na(d$country_rast))
if (length(i) > 0) {
  d <- d[-i, ]
}
cat("  Rows after country_rast removal:", nrow(d), "\n")

cat("  admin1_id NAs:", sum(is.na(d$admin1_id)), "\n")
i <- which(is.na(d$admin1_id))
if (length(i) > 0) {
  d <- d[-i, ]
}
cat("  Rows after admin1_id removal:", nrow(d), "\n")
gc()
}


###########################
### rescale percentages ###
###########################
{
cat("\n=== Rescaling percentage variables (×100) ===\n")
pct_vars <- c(
  "cropland",
  "planted",
  "oil",
  "mining",
  "renewables",
  "roads",
  "powerlines",
  "fires",
  "dams"
)
for (v in pct_vars) {
  d[[v]] <- d[[v]] * 100
}
cat("Rescaled:", paste(pct_vars, collapse = ", "), "\n")
}

###########################
#### new variables ########
###########################
{
  
  
# treat: 1 = study-period PA (2001-2020), 0 = everything else
d$treat <- as.integer(!is.na(d$wdpaid) & d$pa_year_designated >= 2001 & d$pa_year_designated <= 2020)
cat("  treat distribution:\n")
print(table(d$treat, useNA = "always"))

# ever_pa: 1 if pixel is inside any PA regardless of designation year
d$ever_pa <- as.integer(!is.na(d$wdpaid))
cat("  ever_pa distribution:\n")
print(table(d$ever_pa, useNA = "always"))
}




############
### save ###
############
{
cat("\n=== Saving data ===\n")
cat("Total pixels:", format(nrow(d), big.mark = ","), "\n")

write_fst(d, "data/store/pathreat.data.merge.fst", compress = 50)

}

###
# d<-read_fst("data/store/pathreat.data.merge.fst")



