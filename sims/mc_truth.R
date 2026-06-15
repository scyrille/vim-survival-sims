
library(tidyverse)

source(here::here("R", "generate_data.R"))
source(here::here("R", "pred_generators.R"))

n <- 2e7 
set.seed(2026)


# Scenario 1 --------------------------------------------------------------

p <- 4
sigma <- diag(x = 1, nrow = p, ncol = p)
prev <- rep(0.4, p)
alpha <- c(0.1, 0.05, 0.005, 0)
lambda0 <- 0.05
target_cens <- 0.2
c_max <- calibrate_c_max_aalen(lambda0     = lambda0,
                               alpha       = alpha,
                               prev        = prev, 
                               sigma       = sigma,
                               target_cens = target_cens,
                               seed        = 2026)
c_max
data_scen1 <- simulate_aalen_data(n       = n , 
                                  prev    = prev, 
                                  sigma   = sigma, 
                                  lambda0 = lambda0, 
                                  alpha   = alpha, 
                                  c_max   = c_max) 

time <- data_scen1$time
event <- data_scen1$event
X <- data_scen1 %>% dplyr::select(starts_with("X")) 
surv_aalen_fit <- survival::survfit(survival::Surv(time, event)~1)
landmark_time <- quantile(surv_aalen_fit, probs = 0.75)$quantile
approx_times <- seq(0, landmark_time, by = 0.05)
