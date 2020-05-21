
################################################################################
# Create a covid_df data_frame with columns
#   geoid
#   date
#   Y (new_cases)
#   new_cases

# The data should extend from ZERO_PAD (7?) days before the first record to the current time.


proc_covid <- function(raw_data, DATE_0, ZERO_PAD){
require(tidyverse)
#library(lubridate)
#library(covidmodeldata)
################################################################################

if(! ('date' %in% colnames(raw_data))) stop('Expected column called date in raw_date')
if(! ('geoid' %in% colnames(raw_data))) stop('Expected column called geoid in raw_date')
if(! ('total_cases' %in% colnames(raw_data))) stop('Expected column called total_cases in raw_date')

raw_data <- raw_data %>%
  filter(!is.na(geoid))

LAST_DAY = max(raw_data$date)
GEOIDS <- unique(raw_data$geoid)

# Create a data_frame with every day, county, and cases
covid_df <- 
  crossing(
    geoid = unique(raw_data$geoid),
    date  = seq.Date(from=as.Date(DATE_0), to=max(raw_data$date), by=1) ) %>%
  left_join(
    nyt_data %>%
      select(geoid,
             date,
             total_cases)
  )


# Filter to only include days after ZERO_PAD
covid_df <- 
  covid_df %>%
  group_by(
    geoid) %>%
  mutate(
    first_day = min(date[!is.na(total_cases)])
  ) %>%
  filter(
    date >= first_day - ZERO_PAD
  ) %>%
  ungroup() %>%
  arrange(
    geoid, 
    date
  ) 


#########################################################
# Remove the negative cases;
# the way I will do this is by making sure  that the confirmed cases is monotonic, then recomputing new cases as difference

covid_df <- 
  covid_df %>%
  group_by(
    geoid) %>%
  arrange(
    geoid,
    desc(date)) %>%
  mutate(
    future_min = cummin(total_cases)) %>%
  arrange(
    geoid,
    date) %>%
  mutate(
    future_min = ifelse(date < first_day, 0, future_min),
    Y = future_min - lag(future_min, n = 1L, default = 0)
  ) %>%
  ungroup() %>%
  arrange(
    geoid,
    date
  ) %>%
  select(
    -first_day,
    -future_min
  ) %>%
  mutate(t = as.integer(date-as.Date(DATE_0) + 1))
return(covid_df)
}

remove_prisons <- function(covid_df, county_df){
  expected_cols <- c('geoid', 'date', 'Y')
  for( i in expected_cols){
    if ( !(i %in% colnames(covid_df))) stop(sprintf('expected column named %s in covid_df', i))
  }

  expected_cols <- c('geoid', 'pop')
  for( i in expected_cols){
    if ( !(i %in% colnames(county_df))) stop(sprintf('expected column named %s in county_df', i))
  }
  
  # Filter for huge jumps that are likely prisons
  covid_df <- covid_df %>%
    left_join(
      county_df %>%
        select(geoid, pop),
      by='geoid') %>%
    arrange(geoid, date) %>%
    mutate(Yold = Y)
  
  MAX <- 0
  for(i in 2:nrow(covid_df)){
    if(covid_df$geoid[i]!=covid_df$geoid[i-1]) MAX = 0
    if(is.na(covid_df$Y[i])) {next}
    if(covid_df$pop[i] > 30000) {next}
    if (covid_df$Y[i] > 20 & 
        sum(covid_df$Yold[(i-1):(i+1)], na.rm=TRUE) > 50 & 
        sum(covid_df$Yold[(i-1):(i+1)], na.rm=TRUE) > 8*MAX
    ) covid_df$Y[i] <- NA
    if (!is.na(covid_df$Y[i]) & covid_df$Y[i]>MAX) MAX = covid_df$Y[i]
  }
  return(
    covid_df %>% select(-Yold)
  )
}


#' Load Balance Counties across N shards
#' 
assign_shards <- function(covid_df, county_df, n_shards = NSHARDS) {
  
  shard_df <- 
    covid_df %>%
    left_join(
      county_df %>%
        select(geoid,
               mygeoid, mystate,
               state_name, 
               i, county_name, 
               group1, group1_name, 
               group2,group2_name),
      by = "geoid"
      ) %>%
    group_by(
      geoid, mygeoid, mystate, state_name, i, county_name, 
      group1, group1_name, group2, group2_name
      ) %>%
    count() %>%
    group_by(state_name) %>% 
    mutate(shard_name = if_else(group1_name == "not_metro", county_name, group1_name)) %>%
    ungroup()
  
  cut_df <-
    shard_df %>% 
    group_by(state_name, shard_name) %>% 
    summarise(n = sum(n, na.rm = TRUE)) %>%
    group_by(state_name) %>% 
    arrange(n, .by_group = TRUE) %>% 
    mutate(
      shard = dense_rank(cut(cumsum(n), breaks = n_shards, labels = 1:n_shards))
      ) %>%
    select(-n) %>%
    ungroup()
  
  out_df <- left_join(
    shard_df, cut_df, 
    by = c("state_name", "shard_name")
    ) 
  
  # I wasn't sure if you wanted the shard id with the geo information or joined back
  # to the covid_df. If you want it joined back to the covid_df just uncomment the
  # line below
  #
  # out_df <- left_join(covid_df, select(out_df, geoid, shard), by = "geoid")
  #
  out_df
}

