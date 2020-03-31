library(tidycensus)
library(tidyverse)

options(tigris_use_cache = TRUE)
api_key <- Sys.getenv("CENSUS_API_KEY")
census_api_key(api_key)

msa_url <- "https://www2.census.gov/programs-surveys/metro-micro/geographies/reference-files/2018/delineation-files/list1_Sep_2018.xls"
msa_file <- "data-raw/msa-list.xls"
download.file(msa_url, msa_file, mode = "wb")


# download acs ------------------------------------------------------------
acs_vars <- c(
  acs_total_pop = "B25026_001",
  acs_median_income = "B07011_001",
  acs_median_age = "B01002_001"
)

acsdf <- get_acs(
  geography = "county", 
  variables = acs_vars,
  state = NULL, # get everystate
  geometry = TRUE, 
  keep_geo_vars = TRUE,
  output = "wide"
) %>%
  janitor::clean_names() %>%
  select(
    geoid = geoid,
    geoid_full = affgeoid,
    state_fips = statefp,
    county_fips = countyfp,
    county_name = name_x,
    county_name_full = name_y,
    starts_with("acs_")
  )


# join counties to msa data -----------------------------------------------
msadf <- readxl::read_xls(msa_file, skip = 2) %>%
  janitor::clean_names() %>%
  mutate(
    geoid = paste0(fips_state_code, fips_county_code)
  )

acsdf <- left_join(acsdf, msadf, by = "geoid") %>%
  select(
    -(county_county_equivalent:fips_county_code)
  )



# save cleaned up data ----------------------------------------------------
saveRDS(acsdf, "data/us-acs.RData")
