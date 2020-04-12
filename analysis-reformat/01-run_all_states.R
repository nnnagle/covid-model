library(sf)
library(tidyverse)
library(rstan)
library(tidybayes)
library(foreach)
library(doParallel)
library(covidmodeldata)
################################################################################

source('analysis-reformat/00-PARAMS.R')

if(CLEAN_DIR) unlink(DATA_DIR, recursive=TRUE)

if (is.null(NYT_FILE)) {
  
  nyt_data <- get_nyt() %>%
    format_nyt(skip_assignment = c("44")) # don't assign Rhode Island cases
  
} else nyt_data <- read_csv(NYT_FILE)


if (is.null(ACS_FILE)) {
  acs_data <- acs_data # from the covidmodeldata package 
  } else acs_data <- readr::read_rds(ACS_FILE)






######################################################
# COMPILE MODEL
# Negative Binomial with Gaussian Random Walk on slope

stan_mod <- stan_model(file='analysis-reformat/nb_grw.stan', model_name='nb_grw')


######################################################
# GET DATA


# We need a list of states.
# Check the coviddf to make sure that there are enough cases statewide to model.
# Check for only one county, in which case we'll skip.
# Except for DC (11), which we'll process with Maryland (24)

# Step 1. get a list of states from the NYT file
state_list <- nyt_data %>%
  distinct(state_fips) %>%
  filter(!(state_fips %in%  c('11'))) %>% # DC
  pull(state_fips) %>%
  sort()

#state_list <- '53'



#####################################################
# DATE PREP
# CREATE A LIST WITH ONE ENTRY PER CHAIN (50 * NITER)
# 
stan_fit_list <- vector(mode='list', length=length(state_list)*NCHAINS)

slot = 1 # A counter
for(i in 1:length(state_list)){

  #STATE_FIPS = '53' # 53 = Washington
  #STATE_FIPS = '47' 
  # STATE_FIPS='23' Minnesota i=23
  STATE_FIPS = state_list[i]
  # Create sample directory
  STATE_SAMPLES_DIR <- file.path(SAMPLES_ROOT, DATE, STATE_FIPS)
  if(!dir.exists(STATE_SAMPLES_DIR)) dir.create(STATE_SAMPLES_DIR, recursive=TRUE)
  
  geodf <- acs_data %>%
    filter(state_fips == STATE_FIPS)
  
  coviddf <- nyt_data %>% 
    filter(
      state_fips == STATE_FIPS,
      new_cases_mdl >= 0, # using the modeled cases
      date >= as.Date(DATE_0)
      )
  
  # Add DC(11) to Maryland (24)
  if(STATE_FIPS == '24') {
    
    coviddf.11 <-  nyt_data %>% 
      filter(
        state_fips == '11',
        new_cases_mdl >= 0, # using the modeled cases
        date >= as.Date(DATE_0)
      )
    
    coviddf <- rbind(coviddf, coviddf.11)
    rm(coviddf.11)
    
    geodf.11 <- acs_data %>%
      filter(state_fips == '11')
    
    geodf <- rbind(geodf, geodf.11)
    rm(geodf.11)
  }
  
  # extract unique counties from coviddf, create id, and add to geo df:
  covid_cty <- coviddf %>% 
    select(geoid) %>% 
    arrange() %>% 
    unique() %>%
    mutate(
      i = as.integer(as.factor(geoid))
      ) %>%
    left_join(
      geodf %>% st_drop_geometry, by = 'geoid') %>%
    arrange(i)
  
  Xdf <- covid_cty %>%
    select(
      geoid, 
      i, 
      county_name, 
      csa_code, 
      csa_title, 
      pop = acs_total_pop_e, 
      inc = acs_median_income_e
      ) %>%
    mutate(
      pop_s = scale(pop),
      inc_s = scale(inc),
      lpop_s = scale(log(pop))
      )
  
  # Create metro code
  metro_df <- Xdf %>% 
    arrange(csa_code) %>%
    group_by(csa_code) %>%
    summarize(
      ncounties_j = n()
      ) %>%
    mutate(
      msa_j = 1:n()
      )
  
  # If there is a metro with two counties, we can create a level for it.
  # Otherwise, we have to eliminate the metro-level.
  level_1_bool <- metro_df %>% 
    filter(
      !is.na(csa_code)
      ) %>% 
    summarize(
      bool = max(ncounties_j) > 1
      ) %>% 
    pull(bool)
  
  # Insert that id back into the main dataframe, and create a code for NA:
  if(level_1_bool){
    Xdf <- Xdf %>% left_join(metro_df %>% select(csa_code, msa_j), by='csa_code') 
  } else {
    Xdf$msa_j = NA
  }

  # Give the i and t indexes to coviddata
  coviddf <- coviddf %>%
    left_join(
      Xdf %>% select(geoid, i), by='geoid'
      ) %>%
    mutate(
      t = as.numeric(date) - as.numeric(as.Date(DATE_0))+1
      )
  
  stan_data <- list(Pop    = Xdf %>% pull(pop),
                    Y      = coviddf %>% pull(new_cases_mdl), # using modeled cases
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
  
  for(j in 1:NCHAINS){
    stan_fit_list[[slot]] <- list(stan_data=stan_data,
                                  sample_file = file.path(STATE_SAMPLES_DIR, paste0('samples_grw_',j,'.csv')),
                                  diagnostic_file = file.path(STATE_SAMPLES_DIR, paste0('diagnostic_grw_',j,'.csv')),
                                  Xdf = Xdf,
                                  coviddf = coviddf)
    slot = slot + 1
  }
  stan_dat <- stan_fit_list[[slot-1]]
  save(stan_dat, 
       file = file.path(dirname(stan_fit_list[[slot-1]]$sample_file), 'standata.RData'))
}



#################################################################
# Fit model here
######################################################
# Set up cluster here
cl <- makeCluster(NNODES)
registerDoParallel(cl)

r <- foreach(i = 1:length(stan_fit_list),
             .verbose=TRUE, .packages='rstan') %dopar% {
  stan_fit <- try(sampling(object =          stan_mod,
                           data =            stan_fit_list[[i]]$stan_data,
                           iter =            NITER, 
                           thin =            NTHIN,
                           chains =          1,
                           cores =           1,
                           sample_file =     stan_fit_list[[i]]$sample_file,
                           diagnostic_file = stan_fit_list[[i]]$diagnostic_file,
                           append_samples =  FALSE,
                           #init = 0,
                           control =         list(adapt_delta = 0.99, max_treedepth=15)) ) 
}
stopCluster(cl)


