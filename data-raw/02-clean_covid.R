library(tidyverse)

casedf <- readxl::read_excel("data-raw/2020-03-25-washington-county-data.xlsx") %>%
  select(
    state_name,
    counties,
    cases,
    deaths,
    update_time
  ) %>%
  rename(
    county_name = counties,
    total_cases = cases,
    total_deaths = deaths
  ) %>%
  mutate(
    date = lubridate::as_date(update_time)
  ) %>%
  filter(
    !(county_name %in% c("Unassigned", "Total"))
  ) %>%
  select(
    -update_time
  )

coviddf <- readRDS("data/washington-acs.RData") %>%
  select(geoid, county_name) %>%
  sf::st_drop_geometry() %>%
  group_by(geoid) %>%
  expand(
    date = full_seq(c(min(casedf$date),max(casedf$date)), 1),
    county_name
  ) %>%
  ungroup() %>%
  left_join(casedf, by = c("county_name", "date")) %>%
  group_by(geoid) %>%
  arrange(
    date, .by_group = TRUE
  ) %>%
  fill(
    total_cases,
    total_deaths
  ) %>%
  mutate(
    total_cases = if_else(is.na(total_cases), 0, total_cases),
    total_deaths = if_else(is.na(total_deaths), 0, total_deaths),
    new_cases = total_cases - lag(total_cases, 1, default = 0),
    new_deaths = total_deaths - lag(total_deaths, 1, default = 0),
    state_name = "Washington"
  ) 


write_csv(coviddf, "data/washington-covid19.csv")
