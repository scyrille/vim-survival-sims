
#--------------------------#
#  Calibration parameters  #
#--------------------------#

library(here)
library(tidyverse)
library(MASS)

source(here::here("R", "generate_data.R"))

# Large Monte Carlo sample 
n <- 2e7

# Overall seed 
set.seed(2026)

# Scenario 1 ----
param_1 <- calibrate_parameters(n = n, scenario = "1", 
                                max_upper = 1e6, probs = 0.75)

# Scenario 2 ----
param_2 <- calibrate_parameters(n = n, scenario = "2", 
                                max_upper = 1e6, probs = 0.75)


calibration_param <- list(
  params_1 = param_1,
  params_2 = param_2
)

saveRDS(calibration_param, here::here("outputs","results",
                                      "mc_calibration_param.rds"))
