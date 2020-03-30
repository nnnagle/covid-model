
<!-- README.md is generated from README.Rmd. Please edit that file -->

# covid-model

<!-- badges: start -->

<!-- badges: end -->

An initial repo to test spatio-temporal modeling of covid-19.

## Example Data

There are now 4 files under `data/`

### Newer expanded data

  - `2020-03-30-covid19-nyt.csv`: a csv of the latest pull of the [NYT
    county level data](https://github.com/nytimes/covid-19-data) with
    some additional formating. Most states should be fine, but a couple
    of notes.
    
      - New York reports `New York City` with a `geoid = NA` instead of
        the 5 boroughs. Here are the county fips if you want to do
        anything with them `nyc_county_fips <- c("36061", "36047",
        "36081", "36005", "36085")`
      - `Unknown` county data is reported for almost all states.
      - `Kansas City, Missouri` data is reported seperately from the
        counties it is apart of. I wrote some code to distribute its
        counts into the counties using multinomial based on county
        population. I didn’t include it here. Will clean it up and add
        later if wanted

  - `us-acs.RData`: data.frame including `sf` geometries for the whole
    US. Contains example variables from ACS (now including median age)
    as well as county membership of MSAs

### older washingston data

  - `washington-acs.RData`: data.frame including `sf` geometries.
    Contains example variables from ACS as well as county membership of
    MSAs
  - `washington-covid19.csv`: csv of daily county level covid19 cases
    and deaths.

Both of these files can be recreated with the scripts in `data-raw/`

## Daily Mobility data

I didn’t pull it into this but you might be interested in using
[Descartes Labs Mobility
Data](https://github.com/descarteslabs/DL-COVID-19)

``` r
mobility_df <- readr::read_csv("https://github.com/descarteslabs/DL-COVID-19/raw/master/DL-us-mobility-daterow.csv")
#> Parsed with column specification:
#> cols(
#>   date = col_date(format = ""),
#>   country_code = col_character(),
#>   admin_level = col_double(),
#>   admin1 = col_character(),
#>   admin2 = col_character(),
#>   fips = col_character(),
#>   samples = col_double(),
#>   m50 = col_double(),
#>   m50_index = col_double()
#> )

knitr::kable(head(mobility_df))
```

| date       | country\_code | admin\_level | admin1  | admin2 | fips | samples |    m50 | m50\_index |
| :--------- | :------------ | -----------: | :------ | :----- | :--- | ------: | -----: | ---------: |
| 2020-03-01 | US            |            1 | Alabama | NA     | 01   |  133826 |  8.331 |         79 |
| 2020-03-02 | US            |            1 | Alabama | NA     | 01   |  143632 | 10.398 |         98 |
| 2020-03-03 | US            |            1 | Alabama | NA     | 01   |  146009 | 10.538 |        100 |
| 2020-03-04 | US            |            1 | Alabama | NA     | 01   |  149352 | 10.144 |         96 |
| 2020-03-05 | US            |            1 | Alabama | NA     | 01   |  144109 | 10.982 |        104 |
| 2020-03-06 | US            |            1 | Alabama | NA     | 01   |  141491 | 13.024 |        123 |
