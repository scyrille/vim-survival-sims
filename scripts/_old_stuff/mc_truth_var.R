
#-------------------------------#
#  MC TRUE ASYMPTOTIC VARIANCE  #
#-------------------------------#

library(tidyverse)

source(here::here("R", "generate_data.R"))
source(here::here("R", "pred_generators.R"))
source(here::here("R", "mc_truth_helpers.R"))

n <- 1000
seed <- 2026
set.seed(seed)

# Scenario 1 --------------------------------------------------------------

## Generate data with sample splitting  ----
data1_scen1 <- generate_data(n        = n, 
                             scenario = "1",
                             seed     = seed
)
data2_scen1 <- generate_data(n        = n, 
                             scenario = "1",
                             seed     = seed
)
x_vars <- colnames(data1_scen1 %>% dplyr::select(-c(time,event)))

## Fixed time horizons -----
tau <- round(estimate_tau_aalen(
  data    = data_scen1, 
  lambda0 = lambda0, 
  alpha   = alpha,
  seed    = seed, 
  probs   = 0.75
),2)

## Time grid ----
approx_times <- seq(0, 20, by = 0.01)

## Conditional survival function ----
S0_1 <- sapply(approx_times, true_S0_aalen, data1_scen1, lambda0, alpha)
colnames(S0_1) <- approx_times
S0_2 <- sapply(approx_times, true_S0_aalen, data2_scen1, lambda0, alpha)
colnames(S0_2) <- approx_times

## Independent uniform function ----
G0_1 <- sapply(approx_times, true_G0_uniform, data1_scen1, c_max)
colnames(G0_1) <- approx_times
G0_2 <- sapply(approx_times, true_G0_uniform, data2_scen1, c_max)
colnames(G0_2) <- approx_times

## True asymptotic variance ----
#' MC estimate of the true asymptotic variance of the estimator using the EIF 
#' with known nuisances plugged 

### Brier score ----
brier_X1 <- rep(NA, length(tau))
brier_X1_split <- rep(NA, length(tau))

for (t in tau){
  f0_1 <- 1 - S0_1[,t]
  f0_2 <- 1 - S0_2[,t]
  f0_X1_1 <-  1 - true_S0_aalen(
    t       = t,
    data    = data1_scen1 %>% dplyr::select(-X1),
    lambda0 = lambda0,
    alpha   = alpha[2:4]
  )
  f0_X1_2 <-  1 - true_S0_aalen(
    t       = t,
    data    = data2_scen1 %>% dplyr::select(-X1),
    lambda0 = lambda0,
    alpha   = alpha[2:4]
  )
  
  brier_V0_1 <- survML:::estimate_brier(time = data1_scen1$time,
                                        event = data1_scen1$event,
                                        approx_times = approx_times,
                                        tau = t,
                                        preds = f0_1,
                                        S_hat = S0_1,
                                        G_hat = G0_1)
  
  brier_V0_X1_1 <- survML:::estimate_brier(time = data1_scen1$time,
                                           event = data1_scen1$event,
                                           approx_times = approx_times,
                                           tau = t,
                                           preds = f0_X1_1,
                                           S_hat = S0_1,
                                           G_hat = G0_1)
  
  brier_V0_2 <- survML:::estimate_brier(time = data2_scen1$time,
                                        event = data2_scen1$event,
                                        approx_times = approx_times,
                                        tau = t,
                                        preds = f0_2,
                                        S_hat = S0_2,
                                        G_hat = G0_2)
  
  brier_V0_X1_2 <- survML:::estimate_brier(time = data2_scen1$time,
                                           event = data2_scen1$event,
                                           approx_times = approx_times,
                                           tau = t,
                                           preds = f0_X1_2,
                                           S_hat = S0_2,
                                           G_hat = G0_2)
  
  var_X1 <- mean((brier_V0_1$EIF - brier_V0_X1_1$EIF)^2)
  brier_X1[which(tau == t)] <- var_X1
  
  var_X1_split <- mean((brier_V0_2$EIF)^2) + mean((brier_V0_X1_2$EIF)^2)
  brier_X1_split[which(tau == t)] <- var_X1_split
}

brier_scen1 <- data.frame(vim = rep("brier", length(tau)),
                          tau = tau,
                          vim_X1 = brier_X1,
                          vim_X1_split = brier_X1_split)
                     

### AUC ----

auc_X1 <- rep(NA, length(tau))
auc_X1_split <- rep(NA, length(tau))

for (t in tau){
  f0_1 <- 1 - S0_1[,t]
  f0_2 <- 1 - S0_2[,t]
  f0_X1_1 <-  1 - true_S0_aalen(
    t       = t,
    data    = data1_scen1 %>% dplyr::select(-X1),
    lambda0 = lambda0,
    alpha   = alpha[2:4]
  )
  f0_X1_2 <-  1 - true_S0_aalen(
    t       = t,
    data    = data2_scen1 %>% dplyr::select(-X1),
    lambda0 = lambda0,
    alpha   = alpha[2:4]
  )
  
  auc_V0_1 <- survML:::estimate_AUC(time = data1_scen1$time,
                                    event = data1_scen1$event,
                                    approx_times = approx_times,
                                    tau = t,
                                    preds = f0_1,
                                    S_hat = S0_1,
                                    G_hat = G0_1,
                                    robust = TRUE)
  
  auc_V0_X1_1 <- survML:::estimate_AUC(time = data1_scen1$time,
                                       event = data1_scen1$event,
                                       approx_times = approx_times,
                                       tau = t,
                                       preds = f0_X1_1,
                                       S_hat = S0_1,
                                       G_hat = G0_1,
                                       robust = TRUE)
  
  auc_V0_2 <- survML:::estimate_AUC(time = data2_scen1$time,
                                    event = data2_scen1$event,
                                    approx_times = approx_times,
                                    tau = t,
                                    preds = f0_2,
                                    S_hat = S0_2,
                                    G_hat = G0_2,
                                    robust = TRUE)
  
  auc_V0_X1_2 <- survML:::estimate_AUC(time = data2_scen1$time,
                                       event = data2_scen1$event,
                                       approx_times = approx_times,
                                       tau = t,
                                       preds = f0_X1_2,
                                       S_hat = S0_2,
                                       G_hat = G0_2,
                                       robust = TRUE)
  
  var_X1 <- mean((auc_V0_1$EIF - auc_V0_X1_1$EIF)^2)
  auc_X1[which(tau == t)] <- var_X1
  
  var_X1_split <- mean((auc_V0_2$EIF)^2) + mean((auc_V0_X1_2$EIF)^2)
  auc_X1_split[which(tau == t)] <- var_X1_split
}

auc_scen1 <- data.frame(vim = rep("auc", length(tau)),
                        tau = tau,
                        vim_X1 = auc_X1,
                        vim_X1_split = auc_X1_split)
