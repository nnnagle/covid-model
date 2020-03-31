library(tidyverse)

# saved as RData b/c it has sf geometries 
geodf <- readRDS("data/us-acs.RData")

coviddf <- read_csv("data/2020-03-30-covid19-nyt.csv")
