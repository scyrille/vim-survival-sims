
library(tidyverse)

source(here::here("R", "generate_data.R"))
source(here::here("R", "pred_generators.R"))
source(here::here("R", "mc_truth_helpers.R"))

n <- 2e7 
seed <- 2026
set.seed(seed)

# Scenario 1 --------------------------------------------------------------

## Settings 

p <- 4
sigma <- diag(x = 1, nrow = p, ncol = p)
lambda0 <- 0.05
prev <- 0.4
alpha <- c(0.1,    # strong
           0.05,   # moderate
           0.005,  # weak
           0       # null
)       
target_cens <- 0.2

## Calibrate uniform censoring

c_max <- calibrate_c_max_aalen(lambda0     = lambda0,
                               alpha       = alpha,
                               prev        = rep(prev, p), 
                               sigma       = sigma,
                               target_cens = target_cens,
                               seed        = seed)

## Generate data 

data_scen1 <- simulate_aalen_data(n       = n , 
                                  prev    = rep(prev, p), 
                                  sigma   = sigma, 
                                  lambda0 = lambda0, 
                                  alpha   = alpha, 
                                  c_max   = c_max) 


tau <- estimate_tau75_aalen(
  n       = n, 
  lambda0 = lambda0,
  alpha   = alpha,
  prev    = 0.4,
  data    = data_scen1,
  seed    = seed
)


# Time horizon predictions 


## Predictions on time grid 

approx_times <- seq(0, tau, by = 0.05)

# Conditional survival function 
S0_mat <- sapply(
  approx_times,
  true_S0_aalen,
  data = data_scen1,
  lambda0 = lambda0,
  alpha = alpha
)

# Conditional censoring function 
G0_mat <- sapply(
  approx_times,
  true_G0_uniform,
  data = data_scen1,
  c_max = c_max
)

# Full oracle prediction function at tau
f0_tau <- 1 - S0_mat[,tau]

# Residual oracle prediction functions at tau
f0_X1_tau <- 1-true_S0_aalen(
  t       = tau,
  data    = data_scen1[,-3],
  lambda0 = lambda0,
  alpha   = alpha[2:4]
)

## Variable importance

### Brier 
V0_full_brier <- survML:::estimate_brier(
  time         = data_scen1$time,
  event        = data_scen1$event,
  approx_times = approx_times,
  tau          = tau,
  preds        = f0_tau,
  S_hat        = S0_mat,
  G_hat        = G0_mat
)

V0_reduced_X1_brier <- survML:::estimate_brier(
  time         = data_scen1$time,
  event        = data_scen1$event,
  approx_times = approx_times,
  tau          = tau,
  preds        = f0_X1_tau,
  S_hat        = S0_mat,
  G_hat        = G0_mat
)

