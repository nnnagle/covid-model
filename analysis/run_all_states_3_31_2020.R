

DATE = '2020-04-03' # Date of run
NYT_FILE <- '../data/2020-04-03-covid19-nyt.csv'
ACS_FILE <- '../data/us-acs.RData'
DATE_0 <- '2020-03-01' # First date to use
SAMPLES_ROOT <- '../tmp' # All samples will be stored under {SAMPLES_ROOT}/{DATE}/{STATE}
NITER = 3000
NTHIN =3
NCHAINS = 8

################################################################################

library(sf)
library(tidyverse)
library(rstan)
library(tidybayes)

# We need a list of states.
# Check the coviddf to make sure that there are enough cases statewide to model.
# Check for only one county, in which case we'll skip.
# Except for DC (11), which we'll process with Maryland (24)
state_list = read_csv(NYT_FILE) %>%
  select(geoid) %>%
  mutate(state = substr(geoid, start=1, stop=2)) %>% 
  dplyr::filter(!(state %in%  c('11'))) %>% # DC
  pull(state) %>% 
  unique() %>%
  sort()

#state_list <- '53'

for(i in 1:length(state_list)){

#STATE_FIPS = '53' # 53 = Washington
#STATE_FIPS = '47' 
STATE_FIPS = state_list[i]
SAMPLES_DIR <- file.path(SAMPLES_ROOT, DATE, STATE_FIPS)
if(!dir.exists(SAMPLES_DIR)) dir.create(SAMPLES_DIR, recursive=TRUE)

geodf <- readr::read_rds(ACS_FILE) %>%
  filter(state_fips == STATE_FIPS)
coviddf <- read_csv(NYT_FILE) %>% 
  filter(str_starts(geoid, STATE_FIPS)) %>%
  filter(new_cases >= 0) %>%
  filter(date >= as.Date(DATE_0))

# Add DC(11) to Maryland (24)
if(STATE_FIPS=='24'){
  coviddf.11 <-  read_csv(NYT_FILE) %>% 
    filter(str_starts(geoid, '11')) %>%
    filter(new_cases >= 0) %>%
    filter(date >= as.Date(DATE_0))
  coviddf <- rbind(coviddf, coviddf.11)
  rm(coviddf.11)
  geodf.11 <- readr::read_rds(ACS_FILE) %>%
    filter(state_fips == '11')
  geodf <- rbind(geodf, geodf.11)
  rm(geodf.11)
}



# Negative Binomial with Gaussian Random Walk on slope

stan_mod <- stan_model(file='nb_grw.stan', model_name='nb_grw')

# extract unique counties from coviddf:
covid_cty <- coviddf %>% select(geoid) %>% arrange() %>% unique() %>%
  mutate(i = as.integer(as.factor(geoid))) %>%
  left_join(geodf %>% st_drop_geometry, by='geoid') %>%
  arrange(i)
Xdf <- covid_cty %>%
  select(geoid, i, county_name, csa_code, csa_title, pop=acs_total_pop_e, inc=acs_median_income_e) %>%
  mutate(pop_s = scale(pop),
         inc_s = scale(inc),
         lpop_s = scale(log(pop)))
# Create metro code
metro_df <- Xdf %>% select(csa_code) %>%
  arrange(csa_code) %>%
  unique() %>%
  mutate(msa_j = as.integer(as.factor(csa_code)))
# Insert that id back into the main dataframe, and create a code for NA:
Xdf <- Xdf %>% left_join(metro_df, by='csa_code')  %>%
  mutate(msa_j = ifelse(is.na(msa_j), max(msa_j,na.rm=TRUE)+1, msa_j))

# Give the i and t indexes to coviddata
coviddf <- coviddf %>%
  left_join(Xdf %>% select(geoid, i), by='geoid') %>%
  mutate(t = as.numeric(date) - as.numeric(as.Date(DATE_0))+1)

stan_data <- list(Pop    = Xdf %>% pull(pop),
                  Y      = coviddf %>% pull(new_cases),
                  Y_i    = coviddf %>% pull(i),
                  Y_t    = coviddf %>% pull(t),
                  I      = max(coviddf$i),
                  T      = max(coviddf$t),
                  N      = nrow(coviddf),
                  K      = 2,
                  X      = Xdf %>% select(lpop_s, inc_s),
                  J1     = if(all(is.finite(Xdf$msa_j)))
                                   max(Xdf$msa_j) else 0,
                  group1 = if(all(is.finite(Xdf$msa_j)))  
                                   Xdf %>% pull(msa_j) else 
                                   rep(0,max(coviddf$i)),
                  taub_scale = .1,
                  zero_scale=.001 # scaling for the sum-to-zero consraint
                  )

any.na <- sapply(stan_data, function(x) any(is.na(x)))
if(any(any.na)) {
  save.image(file=paste0(SAMPLES_DIR, 'error.RData'))
  warning('na exists in data. Saving workspace to SAMPLES_DIR and exiting')
  next()
}

stan_fit <- try(sampling(object =          stan_mod,
                     data =            stan_data,
                     iter =            NITER, 
                     thin =            NTHIN,
                     chains =          NCHAINS,
                     cores =           NCHAINS,
                     sample_file =     file.path(SAMPLES_DIR, 'samples_grw.csv'),
                     diagnostic_file = file.path(SAMPLES_DIR, 'diagnostic_grw.csv'),
                     append_samples =  FALSE,
                     #init = 0,
                     control =         list(adapt_delta = 0.99, max_treedepth=15)) )
# Save county and covid data to directory
save(Xdf, coviddf, file=file.path(SAMPLES_DIR, 'R_workspace.RData'))
}
