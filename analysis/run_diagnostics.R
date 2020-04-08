ROOT <- '/home/nnagle/Dropbox/students/Piburn/covid-model'
DATA_DIR <- file.path(ROOT, 'tmp','2020-04-04')
WARMUP = 500

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
# FUNCTIONS

# Function to take a state and rood directory, and return a data_frame
#  with one column containing all sample csv paths
get_sample_paths <- function(STATE, ROOT){
  SAMP_DIR <- file.path(ROOT, STATE)
  csv.files <- list.files(SAMP_DIR,
                          pattern='samples.*[0-9]\\.csv$',
                          full.names=TRUE)
  return(tibble(chain_csv = csv.files))
}

# Function to create tidy object from list of files:
read_stan_draws <- function(files, csv_files_col_name='chain_csv', par='*', warmup=500) {
  #Input: files: a dataframe
  #csv_files_col_name: character: names of column with paths to csvs
  # Output: a dataframe with samles
  samples <- files %>%
    mutate(.chain=1:n()) %>% 
    mutate(samples = map(.x=!!sym(csv_files_col_name),
                         ~read_csv(.x, comment='#'))) %>%
    mutate(samples = map(.x=samples, ~mutate(.x,.iteration=1:n() ))) %>%
    unnest(cols=samples) %>%
    mutate(.draw = 1:n()) %>%
    select(.draw, .chain, .iteration, everything()) %>%
    ungroup() %>%
    filter(.iteration > warmup) %>%
    select(.chain, .iteration, starts_with(par)) %>%
    pivot_longer(cols=starts_with(par),names_to='.variable', values_to='.value')
  return(samples)
}

# Function to take a dataframe of samples, and a paramater string and 
#  calculate rhat for all matching pars (start_with(par))
calculate_rhat <- function(samples, warmup=0,par){
  temp <- samples %>%
    filter(.iteration > warmup) %>%
    select(.chain, .iteration, starts_with(par)) %>%
    pivot_longer(cols=starts_with(par),names_to='.variable', values_to='.value')
  NCHAINS = max(temp$.chain)
  NITER = max(temp$.iteration)
  grand_mean <- temp %>% group_by(.variable) %>% summarize(gmean=mean(.value))
  chain_summary <- temp %>% 
    group_by(.chain, .variable) %>% 
    summarize(wmean = mean(.value), wvar = var(.value)) %>%
    ungroup() %>%
    left_join(grand_mean, by='.variable') %>%
    group_by(.variable)   %>%
    summarize(W=mean(wvar), B = NITER /(NCHAINS-1) * sum((wmean-gmean)^2)) %>%
    mutate(V = (1-1/NITER) * W + 1/NITER * B,
           rhat = sqrt(V/W)) %>%
    summarize(max_rhat = max(rhat),
              which_max = .variable[which.max(rhat)],
              q99_rhat = quantile(rhat, .99),
              frac_gt_101 = mean(rhat>1.01))
  return(chain_summary)
}

# Function to take file names and read csv and calculate rhat summary
read_and_summarize <- function(files,  warmup=0, par_select, k=NULL){
  par_expr <- rlang::expr(par_select)
  par_enquo <- rlang::enquo(par_select)
  if(nrow(files)<3) {
    return(tibble(max_rhat = NA*0,
                  which_max = as(NA,'character'),
                  q99_rhat = NA*0,
                  frac_gt_101 = NA*0))
  }
  samples <- files %>%
      mutate(.chain=1:n()) %>%
      #mutate(samples = map(.x=!!sym(csv_files_col_name),
      #                     ~read_csv(.x, comment='#'))) %>%
      mutate(samples = map(.x=chain_csv,
                           ~vroom_stan(.x, col_select=!!par_enquo))) %>%
    mutate(samples = map(.x=samples, ~mutate(.x,.iteration=1:n() ))) %>%
      unnest(cols=samples) %>%
      mutate(.draw = 1:n()) %>%
      select(.draw, .chain, .iteration, everything()) %>%
      ungroup()
  temp <- samples %>%
    filter(.iteration > warmup) %>%
    select(.chain, .iteration, par_select) %>%
    pivot_longer(cols=!!par_enquo,names_to='.variable', values_to='.value')
  NCHAINS = max(temp$.chain)
  NITER = max(temp$.iteration)
  if(!is.null(k)){
    # Delete a chain
    temp <- temp %>% filter(.chain != k)
    NCHAINS=NCHAINS-1
  }
  grand_mean <- temp %>% group_by(.variable) %>% summarize(gmean=mean(.value))
  chain_summary <- temp %>% 
    group_by(.chain, .variable) %>% 
    summarize(wmean = mean(.value), wvar = var(.value)) %>%
    ungroup() %>%
    left_join(grand_mean, by='.variable') %>%
    group_by(.variable)   %>%
    summarize(W=mean(wvar), B = NITER /(NCHAINS-1) * sum((wmean-gmean)^2)) %>%
    mutate(V = (1-1/NITER) * W + 1/NITER * B,
           rhat = sqrt(V/W)) %>%
    summarize(max_rhat = max(rhat),
              which_max = .variable[which.max(rhat)],
              q99_rhat = quantile(rhat, .99),
              frac_gt_101 = mean(rhat>1.01))
  return(chain_summary)
}
  
vroom_stan <- function(file, ...){
  # Stan sample files have commented out lines in the start, middle and end of the file
  # vroom can figure out the start, but not the middle and end
  # Use grep to delete those
  tfile <- paste0(file, '.tmp')
  grepcmd <- paste0("grep -vh '^#' ", file, " > ", tfile)
  system(grepcmd)
  out <- vroom::vroom(tfile, delim=',', num_threads=1, ...)
  unlink(tfile)
  return(out)
}

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


save(diagnostic_df, file=file.path(DATA_DIR, 'diagnostic.Rdata'))

stopCluster(cl)
