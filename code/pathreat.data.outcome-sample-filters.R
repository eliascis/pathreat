###################################################
### pathreat.data.outcome-sample-filters.R ########
### Treated-side outcome target filters ##########
###################################################
#
# Builds treated-side target-filter expressions for each outcome.
# The analysis helper evaluates these only on the protected member of each pair,
# then requires the outcome to be observed for both protected and control members.
# Output: named list of 16 expression strings (one per outcome), NULL where no filter needed.
#
# Country exclusions (NA country, SIDS, no biodiversity) applied upstream in merge.unmatched.R
#
# Requires: data/store/pathreat.data.merge.unmatched.fst
# Produces: data/store/pathreat.data.outcome-sample-filters.Rds

library(dplyr)
library(fst)

source("code/pathreat.analysis.config.R")

####################################
### outcome filters ################
####################################
{
cat("=== Building per-outcome sample filters ===\n")

# load unmatched data (only columns needed for filters)
d <- read_fst(paths$data_unmatched,
              columns = c("country_rast", "country",
                          "planted", "dams", "fires",
                          "mining", "oil", "renewables"))

outcome_flags <- d %>%
  group_by(country_id = country_rast) %>%
  summarise(
    include_planted = as.integer(!any(is.na(planted))),
    include_dams = as.integer(any(dams > 0, na.rm = TRUE)),
    include_fires = as.integer(any(fires > 0, na.rm = TRUE)),
    include_mining = as.integer(any(mining > 0, na.rm = TRUE)),
    include_oil = as.integer(any(oil > 0, na.rm = TRUE)),
    include_renewables = as.integer(any(renewables > 0, na.rm = TRUE)),
    .groups = "drop"
  )

# correct singular cases: US planted data exists but is all NA in source
us_id <- unique(d$country_rast[d$country == "United States"])
outcome_flags$include_planted[outcome_flags$country_id == us_id] <- 1L

cat("Countries:", nrow(outcome_flags), "\n")

# summary
for (v in c("planted", "dams", "fires", "mining", "oil", "renewables")) {
  col <- paste0("include_", v)
  cat(sprintf("  %s: %d included, %d excluded\n", v,
              sum(outcome_flags[[col]], na.rm = TRUE),
              sum(1 - outcome_flags[[col]], na.rm = TRUE)))
}

# build "country_rast %in% c(...)" expression strings from inclusion flags
country_filter <- function(flag_col) {
  ids <- outcome_flags$country_id[outcome_flags[[flag_col]] == 1]
  if (length(ids) == 0) return(NULL)
  paste0("country_rast %in% c(", paste(ids, collapse = ", "), ")")
}

# build treated-side target-filter expressions for each outcome
dep_sample_filters_base <- list(
  built            = NULL,
  cropland         = NULL,
  planted          = country_filter("include_planted"),
  pasture          = NULL,
  oil              = country_filter("include_oil"),
  mining           = country_filter("include_mining"),
  renewables       = country_filter("include_renewables"),
  roads            = NULL,
  powerlines       = NULL,
  fires            = country_filter("include_fires"),
  swu              = NULL,
  dams             = country_filter("include_dams"),
  light            = NULL,
  def0120_parea    = "forest2000_parea > 0.01",
  any_threat       = NULL,
  threat_composite = NULL
)

saveRDS(dep_sample_filters_base, "data/store/pathreat.data.outcome-sample-filters.Rds")
cat("Saved to: data/store/pathreat.data.outcome-sample-filters.Rds\n")

# report non-NULL filters
non_null <- names(Filter(Negate(is.null), dep_sample_filters_base))
cat(sprintf("Filters defined for %d outcomes: %s\n",
            length(non_null), paste(non_null, collapse = ", ")))

rm(d, outcome_flags, us_id, dep_sample_filters_base)
gc()
}
