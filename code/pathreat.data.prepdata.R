library(dplyr)
library(fst)

nowrun=0

################################
### read and manipulate base data ###
################################
{
##pixel pre-match data
d <- readRDS("data/store/pre-process/global_prematch_data.Rds")
class(d)
dim(d)
head(d)
names(d)

##id
nrow(d)
nrow(d) - length(unique(d$pixel_id))
length(unique(d$row_id))
# table(d$pixel_id == d$row_id)
d$row_id <- NULL
head(d)

##rename
d <- rename(d, biodiversity2024 = biodiversity)
head(d)
sum(is.na(d$biodiversity2024)) / nrow(d)


## exclude re-generated data
d$size<-NULL
d$pa_iucn_cat<-NULL
d$pa_year_designated<-NULL
d$wdpaid<-NULL
d$size_class<-NULL
d$iucn_class<-NULL
d$treat<-NULL
d$ever_pa<-NULL
head(d)


## efficient data classes
d$pixel_id <- as.integer(d$pixel_id)
d$biome_raster <- as.integer(d$biome_raster)
# d$treat <- as.integer(d$treat)
d$country_rast <- as.integer(d$country_rast)
d$elevation <- as.integer(d$elevation)
d$annual_total_precipitation <- as.integer(d$annual_total_precipitation)
d$access <- as.integer(d$access)
d$cropsuit2 <- as.integer(d$cropsuit2)
d$country<-NULL
d$country_rast<-NULL
# d$pa_year_designated <- as.integer(d$pa_year_designated)
# d$pa_iucn_cat <- as.integer(d$pa_iucn_cat)
# d$wdpaid <- as.integer(d$wdpaid)

# ##treatment definitions
# table(d$ever_pa,d$treat) #absolute congruent
# # d$treat<-NULL
# class(d$pa_year_designated)
# table(d$pa_year_designated,d$ever_pa,useNA="always")
# i<-which(d$pa_year_designated==9999)
# d[i,"pa_year_designated"]<-NA
# table(is.na(d$wdpaid))
# table(is.na(d$wdpaid),d$ever_pa)


##check missing values
if(nowrun==1){
vlist <- names(d)
for (v in vlist) {
  n <- round(sum(is.na(d[, v])) / nrow(d), 3) * 100
  print(paste0(v, ": ", n))
}
}


## area conversions
area_conversion <- 30 * 30 / 1000^2  # 0.0009 km² per pixel
# create forest2000_parea (forest area in km² per 1km² cell, i.e., area share 0-1)
cat("Creating forest2000_parea from treecover_2000...\n")
d$forest2000_parea <- d$treecover_2000 * area_conversion
print(summary(d$forest2000_parea))
d$treecover_2000<-NULL

# keep x and y coordinates (needed for threat propensity spatial analysis)
# d$x <- NULL
# d$y <- NULL

##save intermediate file to free memory for deforestation loading
# saveRDS(d, "data/tmp/pathreat.data.1_unmatched_temp.Rds")
d.um<-d
rm(d)
gc()



}

################
####saving ####
##############

write_fst(d.um,"data/store/pathreat.data.prepdata.fst")
