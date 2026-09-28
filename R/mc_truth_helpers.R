
# "Oracle" conditional survival function 
true_S0_aalen <- function(t, data, lambda0, alpha) {
  
  X <- data %>% dplyr::select(starts_with("X")) %>% as.matrix()
  haz <- as.numeric(lambda0 + X %*% alpha)
  # exp(-t * haz)
  if (any(haz < 0))
    warning(sum(haz < 0), " subjet with negative hazard.")
  S0 <- exp(-t * pmax(haz, 0))
  
}

true_S0_cox_aalen <- function(t, data, lambda0,
                              alpha, beta,
                              add_vars, mult_vars) {
  
  X_add <- data %>%
    dplyr::select(dplyr::all_of(add_vars)) %>%
    as.matrix()
  
  X_mult <- data %>%
    dplyr::select(dplyr::all_of(mult_vars)) %>%
    as.matrix()
  
  haz_add <- as.numeric(lambda0 + X_add %*% alpha)
  
  if (any(haz_add < 0))
    warning(sum(haz_add < 0),
            " subjects with negative additive hazard.")
  
  haz <- pmax(haz_add, 0) * exp(as.numeric(X_mult %*% beta))
  
  exp(-t * haz)
}
