###################################################
### pathreat.data.merge.unmatched.R ###############
### Apply country-level filters to merged data ####
###################################################
#
# Reads the full merged pixel dataset (merge.fst) and applies
# country-level exclusions to produce the analysis-ready
# unmatched sample (merge.unmatched.fst).
#
# Filters applied:
#   1. SIDS (small island developing states)
#   2. PAs ≤1 km² (ever_pa pixels)
#   3. UNESCO-MAB Biosphere Reserves
#   4. Countries with >=30% biodiversity data missing
#   5. Countries with mean biodiversity <15
#   6. Point-PA pixels (pa_from_point == TRUE)
#
# Note: NA country removal is applied upstream in merge.R.
#
# Requires: data/store/pathreat.data.merge.fst
# Produces: data/store/pathreat.data.merge.unmatched.fst

library(dplyr)
library(fst)

source("code/pathreat.analysis.config.R")

############################
### load data ##############
############################
{
cat("Loading merged data...\n")
d <- read_fst("data/store/pathreat.data.merge.fst")
cat(sprintf("Loaded: %s obs, %d countries\n",
            format(nrow(d), big.mark = ","),
            length(unique(d$country_rast))))
}

# verify pixel_id alignment with random sample (fast check for row-aligned cbinds)
check_pixel_alignment <- function(ids_a, ids_b, n = 1000) {
  idx <- sample(length(ids_a), min(n, length(ids_a)))
  stopifnot(all(ids_a[idx] == ids_b[idx]))
}

#############################
### merge species counts ####
#############################
{
cat("\n=== Merging threatened species counts ===\n")
d_sc <- readRDS("data/store/pathreat.data.species-count.Rds")
cat("  Species-count rows:", nrow(d_sc), "  Data rows:", nrow(d), "\n")
if (nrow(d) == nrow(d_sc)) {
  check_pixel_alignment(d$pixel_id, d_sc$pixel_id)
  d_sc$pixel_id <- NULL
  d <- cbind(d, d_sc)
  cat("  Columns added:", paste(names(d_sc), collapse = ", "), "\n")
} else {
  stop("Row count mismatch. Re-run pathreat.data.species-count.R.")
}
rm(d_sc)
gc()
}

#############################
### merge threat indices ####
#############################
{
cat("\n=== Merging threat indices ===\n")
d_ti <- read_fst("data/store/pathreat.data.threat-indices.fst")
cat("  Threat indices rows:", nrow(d_ti), "  Data rows:", nrow(d), "\n")
if (nrow(d) == nrow(d_ti)) {
  check_pixel_alignment(d$pixel_id, d_ti$pixel_id)
  d_ti$pixel_id <- NULL
  d <- cbind(d, d_ti)
  cat("  Columns added:", paste(names(d_ti), collapse = ", "), "\n")
} else {
  warning("Row count mismatch — skipping. Re-run pathreat.data.threat-indices.R.")
}
rm(d_ti)
gc()
}

######################################
### merge land-use pressure ##########
######################################
{
cat("\n=== Merging land-use pressure ===\n")
d_lp <- read_fst("data/store/pathreat.data.lu.pressure.fst",
                  columns = c("pixel_id", "lu.pressure.est", "lu.pressure.comp", "lu.pressure.access"))
cat("  Pressure rows:", nrow(d_lp), "  Data rows:", nrow(d), "\n")
if (nrow(d) == nrow(d_lp)) {
  check_pixel_alignment(d$pixel_id, d_lp$pixel_id)
  d <- cbind(d,
    lu.pressure.est    = d_lp$lu.pressure.est,
    lu.pressure.comp   = d_lp$lu.pressure.comp,
    lu.pressure.access = d_lp$lu.pressure.access
  )
  cat("  lu.pressure.est NAs:", sum(is.na(d$lu.pressure.est)), "\n")
  cat("  lu.pressure.comp NAs:", sum(is.na(d$lu.pressure.comp)), "\n")
  cat("  lu.pressure.access NAs:", sum(is.na(d$lu.pressure.access)), "\n")
} else {
  warning("Row count mismatch — skipping. Re-run pathreat.data.lu.pressure.calc.R.")
}
rm(d_lp)
gc()
}


#############################
### remove SIDS #############
#############################
{
sids_names <- c(
  "Anguilla", "Antigua and Barbuda", "Aruba", "Bahamas", "Barbados",
  "Belize", "Bermuda", "British Virgin Islands", "Cayman Islands",
  "Cuba", "Curacao", "Dominica", "Dominican Republic", "Grenada",
  "Guadeloupe", "Haiti", "Jamaica", "Martinique", "Montserrat",
  "Puerto Rico", "Saint Kitts and Nevis", "Saint Lucia",
  "Saint Vincent and the Grenadines", "St. Maarten", "Saint-Martin",
  "Saint-Barthelemy", "Trinidad and Tobago", "Turks and Caicos Islands",
  "United States Virgin Islands",
  "American Samoa", "Cook Islands", "Fiji", "French Polynesia",
  "Guam", "Kiribati", "Marshall Islands", "Micronesia",
  "Nauru", "New Caledonia", "Niue", "Northern Mariana Islands",
  "Palau", "Papua New Guinea", "Pitcairn Islands", "Samoa",
  "Solomon Islands", "Timor-Leste", "Tokelau", "Tonga", "Tuvalu",
  "Vanuatu", "Wallis and Futuna Islands",
  "Cabo Verde", "Comoros", "Guinea-Bissau", "Maldives", "Mauritius",
  "Sao Tome and Principe", "Seychelles", "Singapore",
  "Aland Islands", "Faroe Islands", "Gibraltar", "Guernsey",
  "Isle of Man", "Jersey", "Malta", "Monaco", "San Marino",
  "Norfolk Island", "Saint Helena", "Saint Pierre and Miquelon",
  "British Indian Ocean Territory", "Heard I. and McDonald Islands",
  "United States Minor Outlying Islands", "Hong Kong", "Macao"
)
n_before <- nrow(d)
d <- d[!d$country %in% sids_names, ]
cat(sprintf("SIDS removal: %s → %s obs (%d removed)\n",
            format(n_before, big.mark = ","),
            format(nrow(d), big.mark = ","),
            n_before - nrow(d)))
rm(sids_names)
gc()
}



##############################################################
### remove countries with >=30% biodiversity data missing ####
##############################################################
{
bio_missing_threshold <- 30 # percentage
country_lookup <- d %>% distinct(country_rast, country)

bio_stats <- d %>%
  group_by(country_rast) %>%
  dplyr::summarise(
    pct_bio_missing = sum(is.na(biodiversity2024)) / n() * 100,
    .groups = "drop"
  ) %>%
  filter(pct_bio_missing >= bio_missing_threshold) %>%
  left_join(country_lookup, by = "country_rast")

no_bio_ids <- bio_stats$country_rast
n_before <- nrow(d)
d <- d[!d$country_rast %in% no_bio_ids, ]
cat(sprintf("Bio-missing removal (>=%d%% missing): %s → %s obs (%d removed, %d countries)\n",
            bio_missing_threshold,
            format(n_before, big.mark = ","),
            format(nrow(d), big.mark = ","),
            n_before - nrow(d),
            length(no_bio_ids)))
cat("  Excluded:", paste(sprintf("%s (%.0f%%)", bio_stats$country, bio_stats$pct_bio_missing), collapse = ", "), "\n")
rm(bio_stats, no_bio_ids, bio_missing_threshold, country_lookup)
gc()
}


########################################################################
### remove countries with mean biodiversity < 15 #### i.e. Greenland ###
########################################################################
{
bio_low_threshold <- 15 # percentage
country_lookup <- d %>% distinct(country_rast, country)

bio_stats <- d %>%
  group_by(country_rast) %>%
  dplyr::summarise(
    mean_bio = mean(biodiversity2024, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(mean_bio < bio_low_threshold) %>%
  left_join(country_lookup, by = "country_rast")

no_bio_ids <- bio_stats$country_rast
n_before <- nrow(d)
d <- d[!d$country_rast %in% no_bio_ids, ]
cat(sprintf("Bio-low removal (mean <%.0f): %s → %s obs (%d removed, %d countries)\n",
            bio_low_threshold,
            format(n_before, big.mark = ","),
            format(nrow(d), big.mark = ","),
            n_before - nrow(d),
            length(no_bio_ids)))
cat("  Excluded:", paste(sprintf("%s (mean=%.1f)", bio_stats$country, bio_stats$mean_bio), collapse = ", "), "\n")
rm(bio_stats, no_bio_ids, bio_low_threshold, country_lookup)
gc()
}



############################
### remove small PAs #######
############################
{
n_before <- nrow(d)
i <- which(d$ever_pa == 1 & d$size <= 1)
d <- d[-i, ]
cat(sprintf("Small PA removal (size ≤1 km²): %s → %s obs (%d removed)\n",
            format(n_before, big.mark = ","),
            format(nrow(d), big.mark = ","),
            length(i)))
rm(i)
gc()
}


############################################
### remove UNESCO-MAB Biosphere Reserves ###
############################################
{
table(d$desig_eng[grep("UNESCO",d$desig_eng)])

n_before <- nrow(d)
i <- which(d$desig_eng == "UNESCO-MAB Biosphere Reserve" |  d$desig_eng=="Biosphere Reserve-National & UNESCO-MAB")
d <- d[-i, ]
cat(sprintf("UNESCO-MAB Biosphere Reserve removal: %s → %s obs (%d removed)\n",
            format(n_before, big.mark = ","),
            format(nrow(d), big.mark = ","),
            length(i)))
rm(i)
gc()
}


##################################
### remove point-PA pixels #######
##################################
{
n_before <- nrow(d)
i <- which(d$pa_from_point)
d <- d[-i, ]
cat(sprintf("Point-PA removal (pa_from_point == TRUE): %s → %s obs (%d removed)\n",
            format(n_before, big.mark = ","),
            format(nrow(d), big.mark = ","),
            length(i)))
rm(i)
gc()
}


###########################
### some PA analysis ######
###########################
table(d$desig_eng)
x<-d %>% 
  filter(!is.na(wdpaid)) %>% 
  group_by(wdpaid) %>% 
  summarise(
    country = first(country),
    coutnry_rast = first(country_rast),
    gov_type = first(gov_type),
    own_type = first(own_type),
    desig_type = first(desig_type),
    desig_eng = first(desig_eng),
    iucn_class = first(iucn_class)
  ) %>% data.frame()
nrow(x)

table(x$desig_eng,useNA="always")
table(x$gov_type,useNA="always")
table(x$own_type,useNA="always")
table(x$desig_type,useNA="always")
table(x$iucn_class,useNA="always")
table(x$own_type,x$desig_type,useNA="always")
table(x$own_type,x$iucn_class,useNA="always")


#######################
### exclude columns ###
#######################
{
d$pa_from_point                  <- NULL
d$pa_from_point_designation_year <- NULL
}

############
### save ###
############
{
cat("\n=== Saving unmatched analysis sample ===\n")
cat("Total pixels:", format(nrow(d), big.mark = ","), "\n")
cat("Countries:", length(unique(d$country_rast)), "\n")

write_fst(d, "data/store/pathreat.data.merge.unmatched.fst", compress = 50)
cat("Saved to: data/store/pathreat.data.merge.unmatched.fst\n")
}
