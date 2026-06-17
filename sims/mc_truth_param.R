
#-------------------#
#  TRUE VIM VALUES  #
#-------------------#

library(tidyverse)

source(here::here("R", "generate_data.R"))
source(here::here("R", "pred_generators.R"))
source(here::here("R", "mc_truth_helpers.R"))

n <- 2e7
# n <- 1000
seed <- 2026
set.seed(seed)

# Scenario 1 --------------------------------------------------------------

## Generate data
data_scen1 <- generate_data(n        = n, 
                            scenario = "1",
                            seed     = seed
)

## Fixed time horizons -----
tau <- round(estimate_tau_aalen(
  data    = data_scen1, 
  lambda0 = lambda0, 
  alpha   = alpha,
  seed    = seed, 
  probs   = 0.75
),2)

## True VIM values ----
true_param_scen1 <- purrr::map_dfr(tau, function(t) {
  
  y <- as.integer(data_scen1$time > t)
  xvars <- names(dplyr::select(data_scen1, starts_with("X")))
  
  # Full model
  f_full <- true_S0_aalen(
    t       = t,
    data    = data_scen1,
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
      data    = dplyr::select(data_scen1, -dplyr::all_of(xj)),
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
        vim      = "Brier score",
        tau      = t,
        variable,
        V_full   = V_full_brier,
        V_reduced = V_brier,
        vim_value = V_full_brier - V_brier
      ),
    reduced_res %>%
      dplyr::transmute(
        vim      = "AUC",
        tau      = t,
        variable,
        V_full   = V_full_auc,
        V_reduced = V_auc,
        vim_value = pmax(V_full_auc - V_auc, 0)
      )
  )
}) %>%
  dplyr::mutate(n_mc = n,
                scenario_id = "S1")%>%
  relocate(scenario_id)%>%
  arrange(variable, tau)


# Scenario 2 --------------------------------------------------------------

## Generate data
data_scen2 <- generate_data(n        = n, 
                            scenario = "2",
                            seed     = seed
)

## Fixed time horizons -----
tau <- round(estimate_tau_cox_aalen(
  data    = data_scen2, 
  lambda0 = lambda0, 
  alpha   = alpha,
  beta    = beta, 
  seed    = seed,
  probs   = 0.75
),2)

## True VIM values ----
true_param_scen2 <- purrr::map_dfr(tau, function(t) {
  
  y <- as.integer(data_scen2$time > t)
  
  add_vars  <- names(dplyr::select(data_scen2, starts_with("X")))
  mult_vars <- names(dplyr::select(data_scen2, starts_with("Z")))
  all_vars  <- c(add_vars, mult_vars)
  
  f_full <- true_S0_cox_aalen(
    t         = t,
    data      = data_scen2,
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
      data      = dplyr::select(data_scen2, -dplyr::all_of(vj)),
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
        vim       = "Brier score",
        tau       = t,
        variable,
        # type,
        V_full    = V_full_brier,
        V_reduced = V_brier,
        vim_value = pmax(V_full_brier - V_brier, 0)
      ),
    
    reduced_res %>%
      dplyr::transmute(
        vim       = "AUC",
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
    scenario_id = "S2"
  ) %>%
  dplyr::relocate(scenario_id) %>%
  dplyr::arrange(variable, tau)


#--------- Save all files ----------#

true_param <- bind_rows(true_param_scen1, true_param_scen2) 
saveRDS(true_param, "outputs/tables/mc_truth_param.rds")
