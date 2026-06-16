
# Fixed time horizon 
estimate_tau75_aalen <- function(
    n,
    lambda0,
    alpha,
    prev,
    data, 
    seed
) {
  set.seed(seed)
  X <- data %>% dplyr::select(starts_with("X")) %>% as.matrix()
  haz <- as.numeric(lambda0 + X %*% alpha)
  T <- rexp(n, rate = haz)
  as.numeric(quantile(T, probs = c(0.25, 0.75)))
}

# Survival function 
true_S0_aalen <- function(t, data, lambda0, alpha) {
  X <- data %>% dplyr::select(starts_with("X")) %>% as.matrix()
  haz <- as.numeric(lambda0 + X %*% alpha)
  # exp(-t * haz)
  if (any(haz < 0))
    warning(sum(haz < 0), " subjet with negative hazard.")
  exp(-t * pmax(haz, 0))
}

# Cumulative distribution function (CDF) : F(x)=P(X ≤ x)
true_f0_aalen <- function(t, data, lambda0, alpha) {
  1 - true_S0_aalen(t, data, lambda0, alpha)
}

# Independent uniform censoring
true_G0_uniform <- function(t, data, c_max) {
  rep(pmax(1 - t / c_max, 0), nrow(data))
}
