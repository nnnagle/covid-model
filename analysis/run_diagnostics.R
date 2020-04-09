ROOT <- '/home/nnagle/Dropbox/students/Piburn/covid-model'
DATA_DIR <- file.path(ROOT, 'tmp','2020-04-04')
WARMUP = 500
source('functions.R')

# This uses the future library and furrr package for multithreading.
# I use availableCores()-4 on my own machine
# On my machinge with 32 cores (28 used) it takes about 15 minutes (mostly reading)
#
# The output is a small diagnostic.Rdata file in DATA_DIR

library(tidyverse)
library(vroom)
library(furrr)

FIPS='01 02 04 05 06 08 09 10 12 13 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 44 45 46 47 48 49 50 51 53 54 55 56'


cl <- makeClusterPSOCK(availableCores()-4)
plan(cluster, workers = cl)


###############################################################
# Calculate diagnostics
  
# Evaluate Rhat for the parameter bo_raw
# using all chains
diagnostic_df <- tibble(State = strsplit(FIPS, split=' ')[[1]]) %>%
  mutate(files = map(.x=State, ~get_sample_paths(.x, DATA_DIR))) %>%
  mutate(diag = map(.x=files, ~read_and_summarize(.x, WARMUP, 
                                                  par_select=c(starts_with('b0_raw')))))
# Repeat, leaving out chain k
for(k in 1:NCHAINS){
  varname <- paste0('diag', 'k')
  diagnostic_df <- diagnostic_df %>%
    mutate({{varname}} := future_map(.x=files, ~read_and_summarize(.x, 'chain_csv', WARMUP, 'b0_raw',k=k)))
}

###################################################################
# Determine which chains to use here....
# 1. if diag is good, if so, use all chains
# 2. else   Check which diagk has the lowest max Rhat 
#           if min_diagk is ok
#              remove k sample from
#           else break

# REPLACE THIS NEXT LINE!!!!!!!!!!!!
diagnostic_df <- mutate(diagnostic_df, good_files=files)


save(diagnostic_df, file=file.path(DATA_DIR, 'diagnostic.Rdata'))

stopCluster(cl)
