library(sf)
library(tidyverse)
library(rstan)
library(tidybayes)
library(foreach)
library(doParallel)
library(covidmodeldata)
################################################################################

source('analysis-reformat/00-PARAMS.R')

if(CLEAN_DIR) unlink(DATA_DIR, recursive=TRUE)

if (is.null(NYT_FILE)) {
  
  nyt_data <- get_nyt() %>%
    format_nyt(skip_assignment = c("44")) # don't assign Rhode Island cases
  
} else nyt_data <- read_csv(NYT_FILE)


if (is.null(ACS_FILE)) {
  acs_data <- acs_data # from the covidmodeldata package 
  } else acs_data <- readr::read_rds(ACS_FILE)






######################################################
# COMPILE MODEL
# Negative Binomial with Gaussian Random Walk on slope

stan_mod <- stan_model(file='analysis-reformat/nb1_spline.stan', model_name='nb1_spline')


######################################################
# GET DATA


# We need a list of states.
# NO LONGER NECESSARY Check the coviddf to make sure that there are enough cases statewide to model.
# Except for DC (11), which we'll process with Maryland (24)

# Step 1. get a list of states from the NYT file, 
# and order from most counties to least.
state_list <- nyt_data %>%
  group_by(state_fips) %>%
  distinct(county_name) %>%
  summarize(n=n()) %>%
  filter(!(state_fips %in%  c('11'))) %>% # DC 
  filter(!is.na(state_fips)) %>%
  arrange(desc(n)) %>%
  pull(state_fips)

#state_list <- state_list[c(20:25)]


#####################################################
# DATE PREP
# CREATE A LIST WITH ONE ENTRY PER CHAIN (50 * NITER)
# 
stan_fit_list <- vector(mode='list', length=length(state_list)*NCHAINS)

slot = 1 # A counter
for(i in 1:length(state_list)){

  STATE_FIPS = state_list[i]
  # Create sample directory
  STATE_SAMPLES_DIR <- file.path(SAMPLES_ROOT, DATE, STATE_FIPS)
  if(!dir.exists(STATE_SAMPLES_DIR)) dir.create(STATE_SAMPLES_DIR, recursive=TRUE)
  
  geodf <- acs_data %>%
    filter(state_fips == STATE_FIPS)
  
  coviddf <- nyt_data %>% 
    filter(
      state_fips == STATE_FIPS,
      new_cases_mdl >= 0, # using the modeled cases
      date >= as.Date(DATE_0)
      )
  
  # Add DC(11) to Maryland (24) by changing its state_fips
  if(STATE_FIPS == '24') {
    
    coviddf.11 <-  nyt_data %>% 
      filter(
        state_fips == '11',
        new_cases_mdl >= 0, # using the modeled cases
        date >= as.Date(DATE_0)
      )
    
    coviddf <- rbind(coviddf, coviddf.11)
    rm(coviddf.11)
    
    geodf.11 <- acs_data %>%
      filter(state_fips == '11')
    
    geodf <- rbind(geodf, geodf.11)
    rm(geodf.11)
  }
  
  # extract unique counties from coviddf, create id, and add to geo df:
  covid_cty <- coviddf %>% 
    select(geoid) %>% 
    arrange() %>% 
    unique() %>%
    mutate(
      i = as.integer(as.factor(geoid))
      ) %>%
    left_join(
      geodf %>% st_drop_geometry, by = 'geoid') %>%
    arrange(i)
  
  # Create the county-level regression data matrix
  Xdf <- covid_cty %>%
    select(
      geoid, 
      i, 
      county_name, 
      csa_code, 
      csa_title, 
      pop = acs_total_pop_e, 
      inc = acs_median_income_e
      ) %>%
    mutate(
      pop_s = scale(pop),
      inc_s = scale(inc),
      lpop_s = scale(log(pop))
      )
  
  # Create metro code
  metro_df <- Xdf %>% 
    arrange(csa_code) %>%
    group_by(csa_code) %>%
    summarize(
      ncounties_j = n()
      ) %>%
    mutate(
      msa_j = 1:n()
      )
  
  # If there is a metro with two counties, we can create a level for it.
  # Otherwise, we have to eliminate the metro-level.
  level_1_bool <- metro_df %>% 
    filter(
      !is.na(csa_code)
      ) %>% 
    summarize(
      bool = max(ncounties_j) > 1
      ) %>% 
    pull(bool)
  
  # Insert that id back into the main dataframe, and create a code for NA:
  if(level_1_bool){
    Xdf <- Xdf %>% left_join(metro_df %>% select(csa_code, msa_j), by='csa_code') 
  } else {
    Xdf$msa_j = NA
  }

  # # Give the i and t indexes to coviddata
  # coviddf <- coviddf %>%
  #   left_join(
  #     Xdf %>% select(geoid, i), by='geoid'
  #     ) %>%
  #   mutate(
  #     t = as.numeric(date) - as.numeric(as.Date(DATE_0))+1
  #     ) %>%
  #   select(geoid, state_fips, state_name, county_name, date, t, everything())
  
  # I want to add zeros to the data.
  # The data begin with the first observation for each county.
  # This introduces a positive bias to the beginning.
  date_df <- tibble(date = seq.Date(from = as.Date(DATE_0),
                                    to = max(coviddf$date),
                                    by = 1))
  coviddf_template <- crossing(date_df, 
                              coviddf %>% 
                                select(geoid,
                                       state_fips, 
                                       state_name,
                                       county_name) %>% 
                                unique)
  coviddf <- coviddf_template %>%
    left_join(coviddf) %>%
    mutate(new_cases=ifelse(is.na(new_cases),0, new_cases),
           new_cases_mdl = ifelse(is.na(new_cases_mdl), 0, new_cases_mdl))
  
  
  # Add a time integer and county integer
  coviddf <- coviddf %>%
    mutate(t = as.numeric(date)-as.numeric(as.Date(DATE_0))+1) %>%
    left_join(
              Xdf %>% select(geoid, i), by='geoid'
    ) %>%
    select(geoid, i, state_fips, state_name, county_name, date, t, everything())
  
  # That's too many zeros. So only add up to 7 days before the first 1 observation:
  coviddf <- coviddf %>%
    group_by(geoid,state_fips, state_name, county_name) %>%
    mutate(first_t = t[min(which(new_cases_mdl>0))]) %>%
    filter( t>= first_t-7) %>%
    select(-first_t)
  
  
  ######################################
  # Create spline basis for the time series
  # With one knot pert date, the spline is Ba = Ia
  # 1. Set up a random walk smoothing matrix S=DD'
  # i.e. a ~ N(0,DD')
  # 2. Reparameterize so that Ba = Za', a' ~ N(0,I)
  T <- max(coviddf$t) + TPRED
  #spl_K <- T-2 # max: T-2 to use full rank | min: 1 to use parabolic

  ######################################
  # Create day of week design matrix
  X_dow <- tibble(t = seq(1:(T+TPRED))) %>%
    mutate(date = as.Date(DATE_0)+t-1) %>%
    mutate(dow = factor(weekdays(date), 
                        levels= c("Monday",  "Tuesday", "Wednesday", 
                                  "Thursday", "Friday", "Saturday","Sunday" ))) %>%
    model.matrix(~dow-1, data=.)
  
  
    
  ####################################
  # RE method 1
  # Create the intrinsic RW penalty on coef of a Identity matrix spline
  # This will require specifying slope separately
  D <- diff(diag(T+TPRED),differences = 2)
  #D[T-2,] <- D[T-2,]*100 # really penalize curves at the end of time
  D.svd <- svd(t(D))
  Z <- D.svd$u %*% diag(1/D.svd$d)
  Z <- Z[,(ncol(Z)+1-SPL_K):ncol(Z)]
  ####################################
  ## RW method 2
  ## intercept + tp*beta, beta ~ N(0,sigma)
  ## truncated polynomial basis
  ## This requires that there by NO separate slope
  #tp <- spam::toeplitz.spam(y=seq(1:(T+TPRED))-1, x=rep(0,(T+TPRED)))
  #tp <- tp[,1:(T+TPRED-1)]
  ## qr decomposition
  #tp.qr <- qr(tp, LAPACK=FALSE)
  #Q <- qr.Q(tp.qr)
  #R <- qr.R(tp.qr)
  ## svd decomp of R
  #R.svd <- svd(R)
  #Z <- Q %*% R.svd$u %*% diag(R.svd$d)
  #Z <- Z[1:spl_K]
  ## Now, we can do intercept + Za, a~N(0,sigma)
  ## compared to tp, Z has orthogonal columns
  ## PRIOR on variance os spline
  ## Expect the value to grow by 3 orders of magnitude over 30 days
  ## The variance seems to follow a quadratic equation with time, but I
  ## don't know what the parameters are.
  ## Trial and error suggests .02 is a good prior for the variance.
  

  
  stan_data <- list(Pop    = Xdf %>% pull(pop),
                    Y      = coviddf %>% filter(new_cases_mdl>0) %>% pull(new_cases_mdl), # using modeled cases
                    Y_i    = coviddf %>% filter(new_cases_mdl>0) %>% pull(i),
                    Y_t    = coviddf %>% filter(new_cases_mdl>0) %>% pull(t),
                    Y0_i   = coviddf %>% filter(new_cases_mdl==0) %>% pull(i),
                    Y0_t   = coviddf %>% filter(new_cases_mdl==0) %>% pull(t),
                    I      = max(coviddf$i),
                    T      = T+TPRED,
                    N      = nrow(coviddf %>% filter(new_cases_mdl>0) ),
                    N0     = nrow(coviddf %>% filter(new_cases_mdl==0) ),
                    K      = 1,
                    X      = Xdf %>% select(lpop_s),
                    X_dow  = X_dow,
                    spl_K  = SPL_K,
                    Z_spl  = Z,
                    J1     = if(all(is.finite(Xdf$msa_j))) max(Xdf$msa_j) else 0,
                    group1 = if(all(is.finite(Xdf$msa_j))) Xdf %>% pull(msa_j) else  rep(0,max(coviddf$i)),
                    taub_scale = .05,
                    sample_flag = TRUE
  )
  
  for(j in 1:NCHAINS){
    stan_fit_list[[slot]] <- list(stan_data=stan_data,
                                  sample_file = file.path(STATE_SAMPLES_DIR, paste0('samples_grw_',j,'.csv')),
                                  diagnostic_file = file.path(STATE_SAMPLES_DIR, paste0('diagnostic_grw_',j,'.csv')),
                                  Xdf = Xdf,
                                  coviddf = coviddf)
    slot = slot + 1
  }
  stan_dat <- stan_fit_list[[slot-1]]
  save(stan_dat, 
       file = file.path(dirname(stan_fit_list[[slot-1]]$sample_file), 'standata.RData'))
}



#################################################################
# Fit model here
######################################################
# Set up cluster here
cl <- makeCluster(NNODES)
registerDoParallel(cl)

##############################################
# Initialize function
init_fun <- function(chain_id) list( a = runif(1,-13,-9),
                                     b = runif(1,0,.4),
                                     tau_a0 = runif(1,0,2),
                                     tau_a1 = runif(1,0,2),
                                     raw_tau_b0 = runif(1,0,2), 
                                     raw_tau_splb = runif(1,0,2),
                                     phi = runif(1,.1,.9))


r <- foreach(i = 1:length(stan_fit_list),
             .verbose=TRUE, .packages='rstan') %dopar% {
  stan_fit <- try(sampling(object =          stan_mod,
                           data =            stan_fit_list[[i]]$stan_data,
                           iter =            NITER, 
                           thin =            NTHIN,
                           chains =          1,
                           cores =           1,
                           warmup =          WARMUP*NTHIN,
                           sample_file =     stan_fit_list[[i]]$sample_file,
                           diagnostic_file = stan_fit_list[[i]]$diagnostic_file,
                           append_samples =  FALSE,
                           init =            init_fun,
                           control =         list(adapt_delta = 0.99, max_treedepth=10)) ) 
  #return(1)
}
stopCluster(cl)

stan_opt <- try(optimizing(object =          stan_mod,
                         data =            stan_fit_list[[i]]$stan_data,
                         init =            init_fun(),
                         init_alpha =      1e-8,
                         tol_rel_grad =    1e4,
                         algorithm =       "BFGS",
                         verbose=          TRUE)) 

stan_pars <- tibble(value = stan_opt$par, variable = names(stan_opt$par))

