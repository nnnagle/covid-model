library(tidyverse)

# saved as RData b/c it has sf geometries 
geodf <- readRDS("data/washingon-acs.RData")

coviddf <- read_csv("data/washington-covid19.csv")
