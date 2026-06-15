

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
  as.numeric(quantile(T, probs = 0.75))
}


true_S0_aalen <- function(t, data, lambda0, alpha) {
  X <- data %>% dplyr::select(starts_with("X")) %>% as.matrix()
  haz <- as.numeric(lambda0 + X %*% alpha)
  exp(-t * haz)
}

true_f0_aalen <- function(t, data, lambda0, alpha) {
  1 - true_S0_aalen(t, data, lambda0, alpha)
}

true_G0_uniform <- function(t, data, c_max) {
  rep(pmax(1 - t / c_max, 0), nrow(data))
}
