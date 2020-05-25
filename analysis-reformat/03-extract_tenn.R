library(rstan)
library(tidyverse)
library(shinystan)

source('analysis-reformat/00-functions.R')

DATA_DIR <- '/data/covid/tmp/TN/2020-05-24'


# Load in the standata
load(file.path(DATA_DIR, '47', 'standata.RData'))
# Load in ancillary data (shapefiles and what not)
load(file=file.path(DATA_DIR, 'data_frames.Rdata'))
acs_data <- covidmodeldata::acs_data %>% sf::st_drop_geometry()


# Load draws
chain_files <- list.files(file.path(DATA_DIR, '47'), pattern='*.csv',full.names = TRUE)
draws <- read_stan_draws(chain_files,  par_select=c(everything()), warmup=250)


############################################################
# Create county summaries
draws_lambda <- 
  draws %>%
  filter(
    str_detect(.variable, 'log_lambda')) %>%
  separate(
    .variable, 
    sep='\\.', 
    into=c('variable', 'i','t')) %>%
  mutate(
    state_name='Tennessee',
    i = as.integer(i), 
    t=as.integer(t)
  )   %>%
  left_join(
    county_df %>%
      select(state_name, i, county_name, geoid, pop, group1,group2, group1_name, group2_name)
  ) %>%
  left_join(
    date_df %>%
      select(
        t,
        date
      )
  )


county_lambda <- 
  draws_lambda %>%
  group_by(
    variable, 
    county_name,
    geoid,
    date
    ) %>%
  summarize(
    pop = pop[1],
    q50 = exp(quantile(.value, .5)),
    q05 = exp(quantile(.value, .05)),
    q15 = exp(quantile(.value, .15)),
    q25 = exp(quantile(.value, .25)),
    q35 = exp(quantile(.value, .35)),
    q45 = exp(quantile(.value, .45)),
    q55 = exp(quantile(.value, .55)),
    q65 = exp(quantile(.value, .65)),
    q75 = exp(quantile(.value, .75)),
    q85 = exp(quantile(.value, .85)),
    q95 = exp(quantile(.value, .95))) %>%
  ungroup() %>%
  mutate(q50 = q50*pop,
         q05 = q05*pop,
         q15 = q15*pop,
         q25 = q25*pop,
         q35 = q35*pop,
         q45 = q45*pop,
         q55 = q55*pop,
         q65 = q65*pop,
         q75 = q75*pop,
         q85 = q85*pop,
         q95 = q95*pop) %>%
  left_join(
    covid_df %>%
      select(
        geoid,
        date,
        Y
      )
  )

# Merge with case data and set to NA before first period
county_first_case <- county_lambda %>%
  filter(!is.na(Y)) %>%
  group_by(geoid) %>%
  summarize(
    first_case = min(date)
  ) %>%
  ungroup()

county_out <- 
  county_lambda %>%
  left_join(county_first_case) %>%
  filter(date >= first_case) %>%
  select(-first_case)

alpha = .1
county_out %>%
  filter(county_name=='Knox') %>%
  mutate(new_cases = pmax(0, Y)) %>%
  ggplot(aes(x=date)) +
  geom_ribbon(aes(ymin = q05, ymax=q95), alpha=alpha) +
  geom_ribbon(aes(ymin = q15, ymax=q85), alpha=alpha) +
  geom_ribbon(aes(ymin = q25, ymax=q75), alpha=alpha) +
  geom_ribbon(aes(ymin = q35, ymax=q65), alpha=alpha) +
  geom_ribbon(aes(ymin = q45, ymax=q55), alpha=3*alpha) +
  geom_point(aes(y=Y), color='red', alpha=.3) +
  geom_line(aes(y=q50)) +
  labs(
    y='Cases', 
    title = 'Actual and Predicted New Cases in Knox County') 
  


# Merge with geoid, population and count and set lambda to NA if before the first time period.
county_out <- 
  county_lambda %>%
  mutate(
    state_name == 'Tennessee') %>%
  left_join(
    county_df %>%
      select(geoid, state_name, county_name, pop)
  )
  


# Calculate probability of decline over 14 days
trend_14 <- draws_lambda %>%
  filter(county_name=='Knox') %>%
  group_by(county_name,.draw) %>%
  arrange(date) %>%
  mutate(lag14 = pop*(exp(.value) - exp(lag(.value,14)))) %>%
  ungroup() %>%
  filter(date==as.Date('2020-05-16'))

ggplot(data=trend_14,aes(x=lag14)) + 
  geom_histogram(aes(y=..density..)) + 
  labs(x = 'Change in Number of Cases over 14 days',
    title='Knox County: Change in expected cases over last 14 days',
        subtitle=sprintf('There is a %s chance that the rate has been declining over this period',
                         as.character(round(mean(trend_14$lag14<0),3))))




############################################################
# Create county summaries
county_lambda <- 
  draws %>%
  group_by(
    .variable) %>%
  summarize(
    q50 = exp(quantile(.value, .5)),
    q05 = exp(quantile(.value, .05)),
    q95 = exp(quantile(.value, .95))) %>%
  separate(
    .variable, 
    sep='\\.', 
    into=c('variable', 'i','t')) %>%
  filter(
    variable == 'log_lambda')  %>%
  mutate(
    i = as.integer(i), 
    t=as.integer(t),
    variable = 'lambda')






