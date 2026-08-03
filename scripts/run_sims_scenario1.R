
library(here)
library(tidyverse)
library(MASS)
library(survML)
library(SuperLearner)
library(ranger)
library(xgboost)
library(timereg)
library(survivalSL)

fs::dir_create(here::here("outputs", "results", "scenario1"))

source(here::here("R", "generate_data.R"))
source(here::here("R", "pred_generators.R"))
source(here::here("R", "do_one.R"))

# Calibration parameters 
calibration_param <- readRDS(here::here("outputs","results",
                                        "mc_calibration_param.rds"))
c_max <- calibration_param$params_1$c_max
tau <- calibration_param$params_1$tau

# Number of Monte Carlo replications 
R <- 500

# Sample size 
n <- c(500,1000)

nuisances <- c("aalen","stackG","survivalSL")

param_grid <- expand.grid(
  n = n,
  nuisance = nuisances,
  stringsAsFactors = FALSE
)

param_grid$seed <- sample.int(1e9, nrow(param_grid))

purrr::pmap(
  
  param_grid,
  
  function(n, nuisance, seed) {
    
    set.seed(seed)
    
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
  }
)
