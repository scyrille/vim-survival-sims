
#-------------------#
#  TRUE VIM VALUES  #
#-------------------#

library(tidyverse)
library(here)

source(here::here("R", "generate_data.R"))
source(here::here("R", "pred_generators.R"))
source(here::here("R", "mc_truth_helpers.R"))

# Large Monte Carlo sample 
n <- 2e7

# Overall seed 
set.seed(2026)

# Calibration parameters
calibration_param <- readRDS(here::here("outputs","results",
                                        "mc_calibration_param.rds"))

# Scenario 1 --------------------------------------------------------------

## Generate data 
data_1 <- generate_data(n = n, scenario = "1",
                            c_max = calibration_param$params_1$c_max)

## True VIM values 
true_param_scen1 <- purrr::map_dfr(calibration_param$params_1$tau, function(t) {
  
  y <- as.integer(data_1$time > t)
  xvars <- names(dplyr::select(data_1, starts_with("X")))
  
  # Full model
  f_full <- true_S0_aalen(
    t       = t,
    data    = data_1,
    lambda0 = lambda0,
    alpha   = alpha
  )
  
  V_full_brier <- -vimp::measure_mse(f_full, y)$point_est
  V_full_auc   <- cvAUC::AUC(f_full, y)
  
  # Reduced models
  reduced_res <- purrr::map_dfr(seq_along(xvars), function(j) {
    
    xj <- xvars[j]
    
    f_red <- true_S0_aalen(
      t       = t,
      data    = dplyr::select(data_1, -dplyr::all_of(xj)),
      lambda0 = lambda0,
      alpha   = alpha[-j]
    )
    
    tibble::tibble(
      variable = xj,
      V_brier  = -vimp::measure_mse(f_red, y)$point_est,
      V_auc    = cvAUC::AUC(f_red, y)
    )
  })
  
  dplyr::bind_rows(
    reduced_res %>%
      dplyr::transmute(
        vim      = "BS(t)",
        tau      = t,
        variable,
        V_full   = V_full_brier,
        V_reduced = V_brier,
        vim_value = V_full_brier - V_brier
      ),
    reduced_res %>%
      dplyr::transmute(
        vim      = "AUC(t)",
        tau      = t,
        variable,
        V_full   = V_full_auc,
        V_reduced = V_auc,
        vim_value = pmax(V_full_auc - V_auc, 0)
      )
  )
}) %>%
  dplyr::mutate(n_mc = n,
                scenario_id = "1")%>%
  relocate(scenario_id)%>%
  arrange(variable, tau)


# Scenario 2 --------------------------------------------------------------

## Generate data 
data_2 <- generate_data(n = n, scenario = "2",
                        c_max = calibration_param$params_2$c_max)

## True VIM values
true_param_scen2 <- purrr::map_dfr(calibration_param$params_2$tau, function(t) {
  
  y <- as.integer(data_2$time > t)
  
  add_vars  <- names(dplyr::select(data_2, starts_with("X")))
  mult_vars <- names(dplyr::select(data_2, starts_with("Z")))
  all_vars  <- c(add_vars, mult_vars)
  
  f_full <- true_S0_cox_aalen(
    t         = t,
    data      = data_2,
    lambda0   = lambda0,
    alpha     = alpha,
    beta      = beta,
    add_vars  = add_vars,
    mult_vars = mult_vars
  )
  
  V_full_brier <- -vimp::measure_mse(f_full, y)$point_est
  V_full_auc   <- cvAUC::AUC(f_full, y)
  
  reduced_res <- purrr::map_dfr(all_vars, function(vj) {
    
    is_add  <- vj %in% add_vars
    is_mult <- vj %in% mult_vars
    
    add_vars_red  <- setdiff(add_vars, vj)
    mult_vars_red <- setdiff(mult_vars, vj)
    
    alpha_red <- if (is_add) {
      alpha[add_vars != vj]
    } else {
      alpha
    }
    
    beta_red <- if (is_mult) {
      beta[mult_vars != vj]
    } else {
      beta
    }
    
    f_red <- true_S0_cox_aalen(
      t         = t,
      data      = dplyr::select(data_2, -dplyr::all_of(vj)),
      lambda0   = lambda0,
      alpha     = alpha_red,
      beta      = beta_red,
      add_vars  = add_vars_red,
      mult_vars = mult_vars_red
    )
    
    tibble::tibble(
      variable = vj,
      # type     = dplyr::case_when(
      #   is_add  ~ "additive",
      #   is_mult ~ "multiplicative"
      # ),
      V_brier = -vimp::measure_mse(f_red, y)$point_est,
      V_auc   = cvAUC::AUC(f_red, y)
    )
  })
  
  dplyr::bind_rows(
    reduced_res %>%
      dplyr::transmute(
        vim       = "BS(t)",
        tau       = t,
        variable,
        # type,
        V_full    = V_full_brier,
        V_reduced = V_brier,
        vim_value = pmax(V_full_brier - V_brier, 0)
      ),
    
    reduced_res %>%
      dplyr::transmute(
        vim       = "AUC(t)",
        tau       = t,
        variable,
        # type,
        V_full    = V_full_auc,
        V_reduced = V_auc,
        vim_value = pmax(V_full_auc - V_auc, 0)
      )
  )
}) %>%
  dplyr::mutate(
    n_mc = n,
    scenario_id = "2"
  ) %>%
  dplyr::relocate(scenario_id) %>%
  dplyr::arrange(variable, tau)


#--------- Save all files ----------#

true_param <- bind_rows(true_param_scen1, true_param_scen2) 
saveRDS(true_param, here::here("outputs","results","mc_truth_param.rds"))
