#library(sf)
library(tidyverse)
#library(lubridate)
library(rstan)
#library(tidybayes)
library(foreach)
library(doParallel)
library(covidmodeldata)
################################################################################

source('analysis-reformat/00-PARAMS.R')
source('analysis-reformat/00-proc_covid.R')
source('analysis-reformat/00-proc_county.R')

#if(CLEAN_DIR) unlink(DATA_DIR, recursive=TRUE)

if (is.null(NYT_FILE)) {
  
  nyt_data <- get_nyt() %>%
    format_nyt(skip_assignment = c("09","44")) # don't assign Rhode Island cases 
  
} else nyt_data <- read_csv(NYT_FILE)
nyt_data <- mutate(nyt_data, date=as_date(date))


covid_df <- proc_covid(
  raw_data = nyt_data %>% select(geoid, date, total_cases=total_cases_mdl),
  DATE_0 = DATE_0,
  ZERO_PAD = ZERO_PAD
)


if (is.null(ACS_FILE)) {
  acs_data <- acs_data # from the covidmodeldata package 
} else {acs_data <- readr::read_rds(ACS_FILE)}

  
county_df <- proc_county(
  raw_data = acs_data,
  geoid.list = unique(covid_df$geoid)
)

covid_df <- remove_prisons(covid_df, county_df)

date_df <- tibble(
  date = seq.Date(from=as.Date(DATE_0),
                  t = max(covid_df$date)+TPRED,
                  by=1)
  ) %>% 
  mutate(t = as.integer(date-as.Date(DATE_0)+1))


###################################################
# Create the time splines:
X_dow <- date_df %>%
  mutate(dow = factor(weekdays(date), 
                      levels= c("Monday",  "Tuesday", "Wednesday", 
                                "Thursday", "Friday", "Saturday","Sunday" ))) %>%
  arrange(t) %>%
  model.matrix(~dow-1, data=.)

######################################
# Create spline basis for the time series
# With one knot per date, the spline is B * alpha = I * alpha
# The coefficients alpha have a random walk
# 1. Set up a random walk smoothing matrix S=DD'
# i.e. alpha ~ N(0,DD')
# 2. Reparameterize so that Ba = Za', a' ~ N(0,I)

####################################
# RE method 1
# Create the intrinsic RW penalty on coef of a Identity matrix spline
# This will require specifying slope separately
NTIME <- max(date_df$t)
D <- pracma::tril(toeplitz(c(1,-2,1, rep(0, NTIME-4 ))))
svd <- svd((D))
Z <- svd$v %*% diag(1/svd$d)

# Flip order so most important is first:
Z <- Z[,rev(1:ncol(Z))]

ZZ <- Z
# Now, take a low-rank-approximation
Z <- Z[,1:SPL_K]
# Now add the first time to it
Z <- rbind(0, Z)


## PRIOR on variance os spline
## Expect the value to grow by 3 orders of magnitude, over 90 days
## The variance seems to follow a quadratic equation with time, but I
## don't know what the parameters are.
## Trial and error suggests .02 is a good prior for the variance.

# variance on last day is:
max_var <- max(diag(chol2inv(chol(crossprod(D)))))
# sd that should give us log(10^4) change  
sd_scale <- log(10^4) / sqrt(max_var)
  
  
if(!dir.exists(DATA_DIR)) dir.create(DATA_DIR, recursive=TRUE)
save(covid_df, county_df, X_dow, Z, date_df, 
     file=file.path(DATA_DIR, 'data_frames.Rdata'))




# We need a list of states to loop through.
# NO LONGER NECESSARY Check the coviddf to make sure that there are enough cases statewide to model.
# Except for DC (11), which we'll process with Maryland (24)

# Step 1. get a list of states from the NYT file, 
# and order from most counties to least.
state_list <- county_df %>%
  group_by(mystate) %>%
  distinct(county_name) %>%
  summarize(n=n())  %>%
  arrange(desc(n)) %>%
  pull(mystate)

#state_list <- state_list[c(20:25)]

#state_list <- c('15','44','02','56','30')

if(TN_ONLY) state_list <- '47'


#####################################################
# CREATE A LIST WITH ONE ENTRY PER CHAIN (50 * NITER)
# 
stan_fit_list <- vector(mode='list', length=length(state_list)*NCHAINS)

slot = 1 # A counter
for(i in 1:length(state_list)){

  STATE_FIPS = state_list[i]
  
  # Create sample directory
  STATE_SAMPLES_DIR <- file.path(SAMPLES_ROOT, DATE, STATE_FIPS)
  if(!dir.exists(STATE_SAMPLES_DIR)) dir.create(STATE_SAMPLES_DIR, recursive=TRUE)

  
  
    
  ####################
  # Pull out geo data for this state
  this_county_df <- county_df %>%
    filter(mystate == STATE_FIPS) %>%
    arrange(i)
  
  # Pull out covid data for this state and attribute with i, j, and t
  this_covid_df <- covid_df %>%
    filter(geoid %in% this_county_df$geoid) %>%
    left_join(this_county_df %>% select(geoid, i, group1, group2),
              by = 'geoid')
  
  
  # Create the county-level regression data matrix
  Xdf <- this_county_df %>%
    arrange(i) %>%
    select(
      i, 
      group1,
      group2,
      pop
      ) %>%
    mutate(
      pop_s = as.vector(scale(pop)),
      lpop_s = as.vector(scale(log(pop)))
      )
  
  
  stan_data <- list(Pop    = Xdf %>% pull(pop),
                    Y      = this_covid_df %>% filter(Y>0) %>% pull(Y), # using modeled cases
                    Y_i    = this_covid_df %>% filter(Y>0) %>% pull(i),
                    Y_t    = this_covid_df %>% filter(Y>0) %>% pull(t),
                    Y0_i   = this_covid_df %>% filter(Y==0) %>% pull(i),
                    Y0_t   = this_covid_df %>% filter(Y==0) %>% pull(t),
                    I      = max(this_covid_df$i),
                    T      = max(date_df$t),
                    N      = nrow(this_covid_df %>% filter(Y>0) ),
                    N0     = nrow(this_covid_df %>% filter(Y==0) ),
                    K      = 1,
                    X      = Xdf %>% select(lpop_s),
                    X_dow  = X_dow,
                    spl_K  = SPL_K,
                    Z_spl  = Z,
                    J1     = if(max(this_county_df$group1)>1) max(Xdf$group1) else 0,
                    J2     = if(max(this_county_df$group2)>1) max(Xdf$group2) else 0,
                    group1 = if(max(this_county_df$group1)>1) Xdf %>% pull(group1) else  rep(0, nrow(Xdf)),
                    group2 = if(max(this_county_df$group2)>1) Xdf %>% pull(group2) else  rep(0, nrow(Xdf)),
                    taub0_scale = .5*sd_scale,
                    taub1_scale = .5*sd_scale,
                    taub2_scale = sd_scale,
                    sample_flag = TRUE,
                    lppd_flag = FALSE,
                    post_pred = FALSE
  )
  
  for(j in 1:NCHAINS){
    stan_fit_list[[slot]] <- list(stan_data=stan_data,
                                  sample_file = file.path(STATE_SAMPLES_DIR, paste0('samples_grw_',j,'.csv')),
                                  diagnostic_file = file.path(STATE_SAMPLES_DIR, paste0('diagnostic_grw_',j,'.csv')),
                                  Xdf = Xdf,
                                  covid_df = this_covid_df)
    slot = slot + 1
  }
  stan_dat <- stan_fit_list[[slot-1]]
  save(stan_dat, 
       file = file.path(dirname(stan_fit_list[[slot-1]]$sample_file), 'standata.RData'))
}


######################################################
# COMPILE MODEL
# Negative Binomial with Gaussian Random Walk on slope

stan_mod <- stan_model(file='analysis-reformat/nb1_spline.stan', model_name='nb1_spline')





#################################################################
# Fit model here
######################################################
# Set up cluster here
cl <- makeCluster(NNODES)
registerDoParallel(cl)

##############################################
# Initialize function
init_fun <- function(chain_id) list( a = runif(1,-13,-9),
                                     tau_a0 = runif(1,0,.02),
                                     tau_a1 = runif(1,0,.02),
                                     tau_a2 = runif(1,0,2),
                                     raw_tau_splb0 = runif(1,0,.02),
                                     raw_tau_splb1 = runif(1,0,.02),
                                     raw_tau_splb2 = runif(1,0,2),
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
                           append_samples =  FALSE,
                           init =            init_fun,
                           control =         list(adapt_delta = 0.99, max_treedepth=12)) ) 
  #return(1)
}
stopCluster(cl)

#stan_opt <- try(optimizing(object =          stan_mod,
#                         data =            stan_fit_list[[i]]$stan_data,
#                         init =            init_fun(),
#                         init_alpha =      1e-8,
#                         tol_rel_grad =    1e4,
#                         algorithm =       "BFGS",
#                         verbose=          TRUE)) 

#stan_pars <- tibble(value = stan_opt$par, variable = names(stan_opt$par))

