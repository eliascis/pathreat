library(fst)

source("code/pathreat.analysis.config.R")

###############################
### add deforestation data ####
###############################
{
##load deforestation data
cat("Loading deforestation data...\n")
d <- readRDS("data/store/pre-process/global_prematch_deforestation.Rds")
cat("Class:", class(d), "\n")
cat("Dimensions:", dim(d), "\n")

##extract only needed columns (pixel_id + deforestation years)
def_cols <- c(
  "pixel_id",
  names(d)[grep("^deforestation_[0-9]{4}$", names(d))]
)
cat("Selecting", length(def_cols), "columns\n")
d <- d[, def_cols, drop = FALSE]
d <- as.data.frame(d)

##rename deforestation columns
# summary(d$deforestation_2005) 
names(d) <- sub("deforestation_", "def_pixel", names(d))
cat("Final columns:", length(names(d)), "- class:", class(d), "\n")

  cat("Deforestation data processed.\n")
}



##########################################
### create forest area & deforestation ####
##########################################
{
# conversion factor: 30m pixels to km² area
# each 30m pixel = 30*30 m² = 900 m²
# 1 km² = 1,000,000 m²
# so: pixel_count * 30 * 30 / 1000^2 = area in km²
area_conversion <- 30 * 30 / 1000^2  # 0.0009 km² per pixel


# convert yearly deforestation from pixel counts to area (km²)
cat("\nConverting yearly deforestation from pixel counts to area (km²)...\n")
for (yr in 2001:2020) {
  d[[paste0("def_pixel",yr)]] <- d[[paste0("def_pixel",yr)]] * area_conversion
  gc()
}
names(d) <- sub("def_pixel", "def_parea", names(d))
names(d)
summary(d$def_parea2005)


# sum yearly deforestation to get total deforested area 2001-2020
cat("Creating def0120_parea (deforestation as share of pixel area)...\n")
d$def0120_parea <- 0
def_years<-paste0("def_parea",2001:2020)
for (yr in def_years) {
  d$def0120_parea <- d$def0120_parea + ifelse(is.na(d[[yr]]), 0, d[[yr]])
}
print(summary(d$def0120_parea))
gc()

#deforesation area
d$def0120_area<-d$def0120_parea*1000^2
summary(d$def0120_area)
}



#############
### saving ###
##############
head(d)
dim(d)
vlist<-c(
  "pixel_id",
  "def0120_parea",
  "def0120_area"
)
d<-d[,vlist]
write_fst(d, "data/store/pathreat.data.deforestation.fst", compress = 50)
