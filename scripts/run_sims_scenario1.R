
library(here)
library(tidyverse)
library(MASS)
library(survML)
library(SuperLearner)
library(ranger)
library(xgboost)
library(timereg)
library(survivalSL)

source(here::here("R", "generate_data.R"))
source(here::here("R", "pred_generators.R"))
source(here::here("R", "do_one.R"))

# Calibration parameters 
calibration_param <- readRDS(here::here("outputs","results",
                                        "mc_calibration_param.rds"))

# Number of Monte Carlo replications 
# R <- 500
R <- 2

# Sample size 
# n <- c(500,1000,1500)
n <- 100

fs::dir_create(here::here("outputs", "results", "scenario1"))

c_max <- calibration_param$params_1$c_max
tau <- calibration_param$params_1$tau

nuisances <- c("stackG", "aalen")

param_grid <- expand.grid(
  n = n,
  nuisance = nuisances
)

seeds <- sample.int(1e9, nrow(param_grid))

output <- purrr::pmap(
  param_grid,
  function(n, nuisance) {
    
    res <- replicate(
      R,
      do_one(
        n = n,
        scenario = "1",
        c_max = c_max,
        tau = tau,
        nuisance = nuisance
      ),
      simplify = FALSE
    )
    
    saveRDS(
      res,
      here::here("outputs","results","scenario1",
                 paste0("sims_n", n, "_", nuisance,".rds"))
    )
    
    res
  }
)

# names(output) <- paste0(
#   "n", param_grid$n,
#   "_", param_grid$nuisance
# )
