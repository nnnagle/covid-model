# Extract summaries
# INPUTS: 
#   a. the county-level geodf
#   b. the NYT data
#   c. a tibble with one row per state and a <list><tibble> column called "good_files" with a one-column tibble of good files to use (as from diagnostic)
#
# Use run_diagnostics to filter out the files for bad chains
#
# OUTPUT:
#   A sf object, with one row per county x date, with the following summaries
#    - Number of new cases
#    - Number of total cases
#    - modeled new cases that day
#    - q05, q15, q25, q50, q75, q85, q95 of modeled new cases

# NYT_FILE <- '../data/2020-04-03-covid19-nyt.csv'
# ACS_FILE <- '../data/us-acs.RData'
# DIAG_DF_LOC <- '../tmp/2020-04-04/diagnostic.Rdata'

# DATE = '2020-04-04' # Date of run
# DATE_0 <- '2020-03-01' # First date to use
# WARMUP = 500
# SAMPLES_DIR <- file.path('../tmp', DATE)

source('analysis/00-PARAMS.R')
source('analysis/00-functions.R')
# Uses zoom_stan()

library(tidyverse)

########################################################
# Set up a county/date sf object
geodf <- readr::read_rds(ACS_FILE)

coviddf <- read_csv(NYT_FILE) %>%
  filter(new_cases >= 0) %>%
  filter(date >= as.Date(DATE_0)) %>%
  mutate(t = as.numeric(date) - as.numeric(as.Date(DATE_0))+1)

# get a tibble of unique dates, and cartesian product with geodf
date_df <- coviddf %>% dplyr::select(date, t) %>% unique() %>% arrange(t)

out_df <- crossing(geodf, date_df) %>% sf::st_as_sf()
# How do I make this a sf again?


########################################################
FIPS='01 02 04 05 06 08 09 10 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56'

load(DIAG_DF_LOC)

summary_df <- diagnostic_df %>% select(State, good_files)
# A little bit of overhead: counties may be missing from states, so 
# we need to create a crosswalk from the id in the samles to the FIPS
# The necessary info is in standata.Rdata inside each tmp/{DATE}/{STATE} folder

geo_crosswalk_df <- summary_df %>%
  mutate(crosswalk = map(.x=State, function(x){
    path <- file.path(DATA_DIR, as.character(x), 'standata.RData')
    load(path)
    return(stan_dat$Xdf %>% select(geoid, i))
  })) %>%
  select(State, crosswalk) %>%
  unnest(crosswalk)



######################################################################
#options(error=browser)
#options(error=NULL)

summary_df <- summary_df %>%
  mutate(mu_summary = map(good_files, function(x){
    if(length(x)==0){
      return(tibble(mu_mean = NA,
                    mu_q05 = NA,
                    mu_q15 = NA,
                    mu_q25 = NA,
                    mu_q50 = NA,
                    mu_q75 = NA,
                    mu_q85 = NA,
                    mu_q95 = NA))
    }
    samples <- read_stan_draws(x, par_select=c(starts_with("log_mu")), warmup = WARMUP)
    summary <- samples %>% group_by(.variable) %>%
      summarize(mu_mean = mean(exp(.value)),
                mu_q05 = exp(quantile(.value, .05)),
                mu_q15 = exp(quantile(.value, .15)),
                mu_q25 = exp(quantile(.value, .25)),
                mu_q50 = exp(quantile(.value, .5)),
                mu_q75 = exp(quantile(.value, .75)),
                mu_q85 = exp(quantile(.value, .85)),
                mu_q95 = exp(quantile(.value, .95))) 
    summary <- summary %>%
      separate(col = .variable, into=c('.variable','i','t'), sep = '\\.') %>%
      mutate(i = as.integer(i), t=as.integer(t)) %>%
      mutate(.variable = 'mu')
    return(summary)
  }))

mu_out <- summary_df %>%
  select(State, mu_summary) %>%
  unnest(mu_summary)   %>%
  left_join(geo_crosswalk_df, by=c('State','i')) %>%
  mutate(date = as.Date(DATE_0) + t - 1) %>%
  select(geoid, date, everything()) %>%
  arrange(geoid, date)


out_df <- out_df %>% 
  left_join(mu_out, by=c('geoid', 'date'))

save(out_df, file=file.path(RESULTS_DIR, paste0('results_', DATE,'.RData')))
