
#--------------------------#
#  Calibration parameters  #
#--------------------------#

library(tidyverse)

source(here::here("R", "generate_data.R"))

# Large Monte Carlo sample 
n <- 2e7

# Overall seed 
set.seed(2026)

# Scenario 1 ----
params_1 <- calibrate_parameters(n = n, scenario = "1", 
                                 max_upper = 1e6, probs = 0.75)

# Scenario 2 ----
params_2 <- calibrate_parameters(n = n, scenario = "2", 
                                 max_upper = 1e6, probs = 0.75)


calibration_params <- list(
  params_1 = params_1,
  params_2 = params_2
)

saveRDS(calibration_params, here::here("outputs","results",
                                       "mc_calibration_param.rds"))
