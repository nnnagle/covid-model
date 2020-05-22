
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
  
  
  
  group2_df <-
    shard_df %>% 
    group_by(state_name, group2_name) %>% 
    summarise(
      n = sum(n), 
      u = length(na.omit(unique(shard_name)))
      ) %>% 
    group_by(state_name) %>%
    mutate(
      per_n = n / sum(n),
      shard_allocation = if_else(u == 1, 1, as.numeric(NA)),
      shard_allocation = if_else(sum(shard_allocation, na.rm = TRUE) == 1 & is.na(shard_allocation), n_shards - 1, shard_allocation),
      shard_allocation = if_else(is.na(shard_allocation), n_shards * per_n, shard_allocation),
      shard_allocation = if_else(shard_allocation > u, as.numeric(u), shard_allocation),
      shard_allocation = if_else(shard_allocation < 1, as.numeric(1), shard_allocation),
      shard_allocation = if_else(any(shard_allocation == 1, na.rm = TRUE) & shard_allocation != 1, n_shards - 1, shard_allocation),
      shard_allocation = round(shard_allocation)
    ) %>%
    ungroup() %>%
    select(state_name, group2_name, shard_allocation)
  
    shard_df <- left_join(shard_df, group2_df, by = c("state_name", "group2_name"))
  


  cut_df <-
    shard_df %>% 
    group_by(state_name, shard_name, group2_name) %>% 
    summarise(
      n = sum(n, na.rm = TRUE),
      shard_allocation = median(shard_allocation)
    ) 
  
  cut_df_single <- 
    cut_df %>% filter(shard_allocation == 1) %>%
    mutate(
      group2_shard = 1L,
      group2_shard = paste0(group2_name, "-", group2_shard)
    )
  
  cut_df_multi <- 
    cut_df %>% 
    filter(shard_allocation > 1) 
  
  if (nrow(cut_df_multi) > 0) {
  cut_df_multi <-
  cut_df_multi %>%
    group_by(state_name, group2_name) %>% 
    arrange(n, .by_group = TRUE) %>% 
  mutate(
     group2_shard = dense_rank(cut(cumsum(n), breaks = median(shard_allocation))),
     group2_shard = paste0(group2_name, "-", group2_shard)
    )
  }
  
  cut_df <- 
    bind_rows(cut_df_single, cut_df_multi) %>%
    group_by(state_name) %>%
    arrange(group2_shard, .by_group = TRUE) %>%
    mutate(
      shard = cumsum(!duplicated(group2_shard))
    ) %>%
  select(-n, -shard_allocation) %>%
    ungroup()
  
  out_df <- left_join(
   shard_df, cut_df, 
   by = c("state_name","group2_name", "shard_name")
    ) %>% select(-shard_allocation)
  
  # I wasn't sure if you wanted the shard id with the geo information or joined back
  # to the covid_df. If you want it joined back to the covid_df just uncomment the
  # line below
  #
  # out_df <- left_join(covid_df, select(out_df, geoid, shard), by = "geoid")
  #
  out_df
}

