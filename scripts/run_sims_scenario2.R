
library(here)
library(tidyverse)
library(MASS)
library(survML)
library(SuperLearner)
library(ranger)
library(xgboost)
library(timereg)
library(survivalSL)

fs::dir_create(here::here("outputs", "results", "scenario2"))

source(here::here("R", "generate_data.R"))
source(here::here("R", "pred_generators.R"))
source(here::here("R", "do_one.R"))

# Calibration parameters 
calibration_param <- readRDS(here::here("outputs","results",
                                        "mc_calibration_param.rds"))
c_max <- calibration_param$params_2$c_max
tau <- calibration_param$params_2$tau

# Number of Monte Carlo replications 
R <- 500

# Sample size 
n <- c(500,1000)

nuisances <- c("cox.aalen", "stackG", "survivalSL")

param_grid <- expand.grid(
  n = n,
  stringsAsFactors = FALSE
)

param_grid$seed <- sample.int(1e9, nrow(param_grid))

purrr::pmap(
  
  param_grid,
  
  function(n, seed) {
    
    set.seed(seed)
    
    res <- replicate(
      R,
      do_one(
        n = n,
        scenario = "2",
        c_max = c_max,
        tau = tau,
        nuisance = nuisances
      ),
      simplify = FALSE
    )
    
    saveRDS(
      res,
      here::here("outputs","results","scenario2",
                 paste0("sims_n", n, ".rds"))
    )
  }
)
