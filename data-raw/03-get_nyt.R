#' Get Latest COVID-19 Data from New York Times Database
#'
#' @param level character. county level data or state level
#'
#' @return df
#' @export
get_nyt <- function(level = "county") {
  
  if (level == "county") {
    url <- "https://github.com/nytimes/covid-19-data/raw/master/us-counties.csv"
  }
  
  if (level == "state") {
    url <- "https://github.com/nytimes/covid-19-data/raw/master/us-states.csv"
  }
  
  df <- readr::read_csv(url)
  df <- janitor::clean_names(df)
  
  df
}


library(tidyverse)

df <- get_nyt() %>%
  transmute(
    geoid = fips,
    date = date,
    county_name = county,
    state_name = state,
    total_cases = cases,
    total_deaths = deaths
  ) %>%
  group_by(
    county_name, state_name, geoid # some geoids are NA so so county and state
  ) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(
    new_cases = total_cases - lag(total_cases, 1, default = 0),
    new_deaths = total_deaths - lag(total_deaths, 1, default = 0),
  )

write_csv(df, paste0("data/", Sys.Date(), "-covid19-nyt.csv"))
  




