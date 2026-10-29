
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


# Reduced oracle for the independent covariates in S2_base and S2_prev.
# Marginalize an omitted binary covariate over its two possible values.
# In S2_base, the continuous covariate Z3 has coefficient zero.
true_S0_reduced_cox_aalen <- function(t, data, lambda0,
                                      alpha, beta,
                                      add_vars, mult_vars,
                                      variable, prev_add, prev_mult) {
  
  is_add <- variable %in% add_vars
  is_mult <- variable %in% mult_vars
  
  if (length(variable) != 1L || !xor(is_add, is_mult)) {
    stop("variable must identify one additive or multiplicative covariate.")
  }
  
  j <- if (is_add) match(variable, add_vars) else match(variable, mult_vars)
  coefficient <- if (is_add) alpha[j] else beta[j]
  
  # A null covariate leaves the full oracle prediction unchanged.
  if (coefficient == 0) {
    return(true_S0_cox_aalen(
      t = t, data = data, lambda0 = lambda0,
      alpha = alpha, beta = beta,
      add_vars = add_vars, mult_vars = mult_vars
    ))
  }
  
  if (is_add) {
    p_j <- prev_add[j]
    
    # Compute the rate without the omitted covariate directly, so that
    # predictions cannot depend on its observed value, even through rounding.
    add_vars_red <- add_vars[-j]
    haz_add_0 <- as.numeric(
      lambda0 + as.matrix(data[, add_vars_red, drop = FALSE]) %*% alpha[-j]
    )
    mult_rate <- exp(as.numeric(
      as.matrix(data[, mult_vars, drop = FALSE]) %*% beta
    ))
    
    S_0 <- exp(-t * pmax(haz_add_0, 0) * mult_rate)
    S_1 <- exp(-t * pmax(haz_add_0 + alpha[j], 0) * mult_rate)
    
  } else {
    if (j > length(prev_mult)) {
      stop(
        "A continuous omitted covariate with a nonzero coefficient requires ",
        "integration over its conditional distribution."
      )
    }
    
    p_j <- prev_mult[j]
    mult_vars_red <- mult_vars[-j]
    lp_mult_0 <- as.numeric(
      as.matrix(data[, mult_vars_red, drop = FALSE]) %*% beta[-j]
    )
    haz_add <- pmax(as.numeric(
      lambda0 + as.matrix(data[, add_vars, drop = FALSE]) %*% alpha
    ), 0)
    
    S_0 <- exp(-t * haz_add * exp(lp_mult_0))
    S_1 <- exp(-t * haz_add * exp(lp_mult_0 + beta[j]))
  }
  
  if (length(p_j) != 1L || !is.finite(p_j) || p_j < 0 || p_j > 1) {
    stop("The omitted binary covariate must have a prevalence in [0, 1].")
  }
  
  (1 - p_j) * S_0 + p_j * S_1
}
