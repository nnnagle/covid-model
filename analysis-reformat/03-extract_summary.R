# OUTPUT:
#   A data_frme with one row per county x date, with the following summaries
#    - Number of new cases
#    - Number of total cases
#    - lambda: modeled new cases that day
#    - lambda_q05, _q15, _q25, _q50, _q75, _q85, _q95L quantiles of modeled new cases
#    - Ysim - simulated cases and quantiles


#source('analysis-reformat/00-PARAMS.R')
source('/data/covid/tmp/2020-07-10/00-PARAMS.R')


source('analysis-reformat/00-functions.R')
# Uses zoom_stan()

library(tidyverse)
library(sf)
library(future)
library(furrr)


cl <- makeClusterPSOCK(NNODES)
plan(cluster, workers = cl)

load(file=file.path(DATA_DIR, 'data_frames.Rdata'))

FIPS='01 02 04 05 06 08 09 10 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56'

load(DIAG_DF_LOC)

summary_df <- diagnostic_df %>% select(State, good_files)


#  Collect summary stats for lambda: cases per person
summary_lambda <- 
  summary_df %>%
  mutate(
    lambda_summary = future_map(good_files, function(x){
      if(length(x)==0){
        return(tibble(lambda_mean = NA,
                      lambda_q05 = NA,
                      lambda_q15 = NA,
                      lambda_q25 = NA,
                      lambda_q50 = NA,
                      lambda_q75 = NA,
                      lambda_q85 = NA,
                      lambda_q95 = NA))
      }
      samples <- read_stan_draws(x, par_select=c(starts_with("log_lambda")), warmup = WARMUP)
      summary <- samples %>% 
        group_by(.variable) %>%
        summarize(lambda_mean = mean(exp(.value)),
                  lambda_q05 = exp(quantile(.value, .05)),
                  lambda_q15 = exp(quantile(.value, .15)),
                  lambda_q25 = exp(quantile(.value, .25)),
                  lambda_q50 = exp(quantile(.value, .5)),
                  lambda_q75 = exp(quantile(.value, .75)),
                  lambda_q85 = exp(quantile(.value, .85)),
                  lambda_q95 = exp(quantile(.value, .95))) 
      summary <- summary %>%
        separate(col = .variable, into=c('.variable','i','t'), sep = '\\.') %>%
        mutate(
          i = as.integer(i), 
          t=as.integer(t)) %>%
        mutate(
          .variable = 'lambda')
      return(summary)
    }
    )
  ) %>%
  select(
    State, 
    lambda_summary) %>%
  unnest(
    lambda_summary)   %>%
  select(-.variable) 


#  Collect summary stats for lambda: cases per person
summary_Ysim <- 
  summary_df %>%
  mutate(
    Ysim_summary = future_map(good_files, function(x){
      if(length(x)==0){
        return(tibble(Ysim_mean = NA,
                      Ysim_q05 = NA,
                      Ysim_q15 = NA,
                      Ysim_q25 = NA,
                      Ysim_q50 = NA,
                      Ysim_q75 = NA,
                      Ysim_q85 = NA,
                      Ysim_q95 = NA))
      }
      samples <- read_stan_draws(x, par_select=c(starts_with("Y_sim")), warmup = WARMUP)
      summary <- samples %>% 
        group_by(.variable) %>%
        summarize(Ysim_mean = mean(exp(.value)),
                  Ysim_q05 = as.integer(quantile(.value, .05)),
                  Ysim_q15 = as.integer(quantile(.value, .15)),
                  Ysim_q25 = as.integer(quantile(.value, .25)),
                  Ysim_q50 = as.integer(quantile(.value, .5)),
                  Ysim_q75 = as.integer(quantile(.value, .75)),
                  Ysim_q85 = as.integer(quantile(.value, .85)),
                  Ysim_q95 = as.integer(quantile(.value, .95))) 
      summary <- summary %>%
        separate(col = .variable, into=c('.variable','i','t'), sep = '\\.') %>%
        mutate(
          i = as.integer(i), 
          t=as.integer(t)) %>%
        mutate(
          .variable = 'Ysim')
      return(summary)
    }
    )
  ) %>%
  select(
    State, 
    Ysim_summary) %>%
  unnest(
    Ysim_summary) %>%
  select(-.variable)



##############################################
# Create a tibble with one row per county/date
########################################################
data_out <- 
  crossing(
    county_df %>% select(geoid, mygeoid, i, mystate, state_name, county_name, pop), 
    date_df)   %>%
  left_join(
    summary_lambda,
    by = c('mystate'='State', 'i', 't')
  ) %>%
  left_join(
    summary_Ysim,
    by = c('mystate'='State', 'i', 't')
  ) %>%
  left_join(
    covid_df %>% select(-pop),
    by = c('geoid', 'date','t')
  ) %>%
  select(-mygeoid, -mystate, -i, -t)


# Save the county data.
save(data_out, file=file.path(RESULTS_DIR, paste0('results_', DATE,'.RData')))


#################################################################################
# create state summary
summary_mu <- 
  summary_df %>%
  mutate(
    mu_summary = future_map2(.x=State, .y=good_files, function(x,y){
      if(length(y)==0){
        return(tibble(mu_mean = NA,
                      mu_q05 = NA,
                      mu_q15 = NA,
                      mu_q25 = NA,
                      mu_q50 = NA,
                      mu_q75 = NA,
                      mu_q85 = NA,
                      mu_q95 = NA))
      }
      samples <- read_stan_draws(y, par_select=c(starts_with("log_lambda")), warmup = WARMUP)
      summary <- samples %>% 
        separate(.variable, into=c('.variable','i','t'), sep = '\\.') %>%
        mutate(i=as.integer(i)) %>%
        mutate(State=x) %>%
        left_join(
          county_df %>% select(mystate, state_name, i, pop),   
          by = c('State'='mystate', 'i')) %>%
        group_by(.draw, t, state_name) %>%
        summarize(.value = sum(exp(.value)*pop, na.rm=TRUE)) %>%
        drop_na() %>%
        ungroup() %>% 
        group_by(t, state_name) %>%
        summarize(
                  mu_q05 = (quantile(.value, .05)),
                  mu_q15 = (quantile(.value, .15)),
                  mu_q25 = (quantile(.value, .25)),
                  mu_q50 = (quantile(.value, .5)),
                  mu_q75 = (quantile(.value, .75)),
                  mu_q85 = (quantile(.value, .85)),
                  mu_q95 = (quantile(.value, .95))) %>%
        ungroup()
      summary <- summary %>%
        mutate(
          t=as.integer(t)) 
      return(summary)
    }
    )
  ) %>%
  select(
    State, 
    mu_summary) %>%
  unnest(
    mu_summary)

# Fix the FIPS code for DC
summary_mu %>%
  mutate(State=if_else(state_name=='District of Columbia','11', State ))

state_out <- 
  summary_mu %>%
  left_join(county_df %>%
              select(state_name, pop) %>% 
              group_by(state_name) %>%
              summarize(pop = sum(pop)),
            by=c('state_name')) %>%
  left_join(date_df, by='t') %>%
  arrange(State, t) %>%
  select(-t)

# Save the state data.
save(state_out, file=file.path(RESULTS_DIR, paste0('results_state_', DATE,'.RData')))
