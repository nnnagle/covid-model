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

source('analysis-reformat/00-PARAMS.R')
source('analysis-reformat/00-functions.R')
# Uses zoom_stan()

library(tidyverse)
library(future)
library(furrr)
library(covidmodeldata)

cl <- makeClusterPSOCK(NNODES)
plan(cluster, workers = cl)


if (is.null(NYT_FILE)) {
  
  nyt_data <- get_nyt() %>%
    format_nyt(skip_assignment = c("44")) # don't assign Rhode Island cases
  
} else nyt_data <- read_csv(NYT_FILE)

if (is.null(ACS_FILE)) {
  acs_data <- acs_data # from the covidmodeldata package 
} else acs_data <- readr::read_rds(ACS_FILE)

geodf <- acs_data

coviddf <- nyt_data %>%
  filter(new_cases >= 0) %>%
  filter(date >= as.Date(DATE_0)) %>%
  mutate(t = as.numeric(date) - as.numeric(as.Date(DATE_0))+1)


# get a tibble of unique dates, and cartesian product with geodf
T <- max(coviddf$t)
date_df <- tibble(t=seq(1,T+TPRED)) %>%
  mutate(
    date = as.Date(DATE_0)+t-1)

#date_df <- coviddf %>% dplyr::select(date, t) %>% unique() %>% arrange(t)

########################################################
# Set up a county/date sf object
out_df <- crossing(geodf, date_df) %>% sf::st_as_sf() %>%  left_join(coviddf)



########################################################
# DC is intentionally missing from the following FIPS b/c it is with Maryland
FIPS='01 02 04 05 06 08 09 10 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56'

load(DIAG_DF_LOC)

summary_df <- diagnostic_df %>% select(State, good_files)
# A little bit of overhead: counties may be missing from states, so 
# we need to create a crosswalk from the id in the samles to the FIPS
# The necessary info is in standata.Rdata inside each tmp/{DATE}/{STATE} folder

geo_crosswalk_df <- summary_df %>%
  mutate(
    crosswalk = map(.x=State, function(x){
      path <- file.path(DATA_DIR, as.character(x), 'standata.RData')
      load(path)
      return(stan_dat$Xdf %>% select(geoid, i))
    })) %>%
  select(
    State, 
    crosswalk) %>%
  unnest(
    crosswalk)



######################################################################
#options(error=browser)
#options(error=NULL)
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
  )

##############################################
# Create an tibble with one row per county/date
lambda_out <- summary_lambda %>%
  select(
    State, 
    lambda_summary) %>%
  unnest(
    lambda_summary)   %>%
  left_join(
    geo_crosswalk_df, by=c('State','i')) %>%
  mutate(
    date = as.Date(DATE_0) + t - 1) %>%
  select(
    geoid, 
    date, 
    everything()) %>%
  arrange(
    geoid, 
    date)

# Merge into sf object
out_df <- crossing(geodf, date_df) %>% sf::st_as_sf() %>%  left_join(coviddf)
out_df <- out_df %>% 
  left_join(lambda_out, by=c('geoid', 'date'))

# normalize lambda
out_df <-
  out_df %>%
  group_by(geoid) %>%
  mutate(fudge = sum(new_cases_mdl,na.rm=TRUE)/(sum(lambda_q50,na.rm=TRUE)*acs_total_pop_e)) %>%
  ungroup()


########################################################
# Create a panel plot of all counties for a state.
ggplot(
  data= out_df %>% 
    filter(state_fips=='54'), 
  aes( 
    x=date, 
    y=fudge*lambda_q50*10000,
    ymax = fudge*lambda_q85*10000,
    ymin = fudge*lambda_q15*10000)) + 
  geom_ribbon( 
    fill='grey50') +
  geom_line() + 
  geom_point(aes(y=(new_cases_mdl)/acs_total_pop_e*10000),alpha=.1,color='red') +
  facet_wrap(
    ~county_name,
    scales='free_y') +
  scale_y_continuous('New Cases per 10,000 people per day') +
  scale_x_date(
    limits = c(min(date_df$date+7), max(coviddf$date))) + 
  theme(
    axis.text.x = element_text(angle=45, hjust = 1, vjust=1)
  )



##############################################
# Create a state_boundary layer for reference
state <- geodf %>%
  group_by(state_fips, state_name) %>%
  summarize()






#################################################
# Calculate doubling time.
# Calcualate mu today, mu one week ago, and then rate, and then convert to doubling time
# rate = log(mu[t])-log(mu[t-1])
# double time = log(2)  / rate

lambda_out %>% 
  filter( 
    (date == as.Date(DATE_0)+T-1) | (date == as.Date(DATE_0)+T-1-7) ) %>%
  group_by(
    geoid, State, i) %>%
  arrange(
    date) %>%
  summarize(
    gr = log(lambda_q50[2]) - log(lambda_q50[1]),
    rate = lambda_q50[2]) %>%
  mutate(
    double = log(2)/gr) %>%
  mutate(
    double_c = cut(gr, 
                   breaks=c(-Inf, -log(2)/1, -log(2)/4,  log(2)/4, log(2)/1,Inf),
                   labels=c('Halve in < 1 week', 'Halve in 1-4 weeks', 'No Change', 'Double in 1-4 weeks', 'Double in <1 week')),
    rate_c = cut(rate*10000,
                 breaks=c(0, .01, .1, 1,Inf),
                 labels = c('Less than 1 in 1,000,000','1 in 1,000,000', '1 in 100,000','Greater than 1 in 10,000'))
  ) %>%
  left_join(geodf) %>%
  st_as_sf() -> map_temp

##########################################################################
# Cross-tab the doubling time by the rate.
tab <- table(map_temp$double_c, map_temp$rate_c)
tab

#######################################################
# Create the doubling time map.
ggplot(
  data = map_temp,
  aes = aes( 
    fill = double_c)) + 
  geom_sf(aes(fill=double_c), color=NA) +
  scale_fill_brewer('doubling time',
                    direction = -1,
                    palette = 'RdYlBu') + 
  geom_sf(data=state, color='black', fill=NA) +
  theme( 
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank(),
    axis.title = ggplot2::element_blank(),
    panel.grid = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank()) + 
  labs(title='Number of weeks for new cases to double/halve (red/blue)',
       subtitle = 'Rate calculated Apr 8 - Apr 15')
  
#######################################################
# Save the data.
save(out_df, file=file.path(RESULTS_DIR, paste0('results_', DATE,'.RData')))
