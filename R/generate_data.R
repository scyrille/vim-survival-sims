#' Generate a block-structured correlation matrix
#'
#' @param p Total number of variables.
#' @param block_size Size of each block (assumed equal size).
#' @param rho_within Numeric vector of within-block correlations.
#'        Length must equal the number of blocks.
#' @param rho_between Correlation between blocks.
#'
#' @return A p x p correlation matrix.
#'
#' @export
generate_block_corr <- function(p = 20, 
                                block_size = 5, 
                                rho_within = c(0.6, 0.5, 0.4, 0.3),
                                rho_between = 0.1){
  
  # Number of blocks
  n_blocks <- ceiling(p / block_size)
  
  if(length(rho_within) != n_blocks){
    stop("Length of rho_within must equal number of blocks.")
  }
  
  # Initialize matrix with between-block correlation
  Sigma <- matrix(rho_between, p, p)
  diag(Sigma) <- 1
  
  # Define blocks
  blocks <- split(1:p, ceiling((1:p)/block_size))
  
  # Fill each block with its own rho
  for(i in seq_along(blocks)){
    b <- blocks[[i]]
    Sigma[b, b] <- rho_within[i]
    diag(Sigma)[b] <- 1
  }
  
  return(Sigma)
}


#' Simulate correlated binary variables via a latent Gaussian model.
#'
#' @param n Integer. Sample size.
#' @param prev Numeric vector of length `p` with marginal prevalences of DNA alterations (values in (0,1)).
#' @param Sigma Numeric `p x p` correlation (or covariance) matrix for the latent Gaussian variables.
#'
#' @return A data.frame with `p` binary variables named `X1, ..., Xp`.
#'
#' @export
simulate_dna_pathways <- function(n, 
                                  prev, 
                                  Sigma){
  
  p <- length(prev)
  
  stopifnot(is.numeric(prev), all(prev > 0 & prev < 1))
  stopifnot(is.matrix(Sigma), nrow(Sigma) == p, ncol(Sigma) == p)
  
  # Latent Gaussian generation
  Z <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  # Thresholds matching marginal prevalences
  thresholds <- qnorm(prev)
  
  # Dichotomization
  X <- sweep(Z, 2, thresholds, "<") * 1L
  
  colnames(X) <- paste0("D", seq_len(p))
  as.matrix(X)
}

#' Simulate GSVA-like RNA pathway scores
#'
#' @param n Number of samples
#' @param p Number of RNA pathways (default = 20)
#' @param mu Mean
#' @param Sigma Numeric `p x p` correlation (or covariance) matrix for the latent Gaussian variables.
#'
#' @return A standardized matrix (n × p) of RNA pathway scores
#' 
#' @export
simulate_rna_pathways <- function(n,
                                  p = 20,
                                  mu = 0, 
                                  Sigma) {
  
  # --- Multivariate normal simulation ---
  raw_scores <- MASS::mvrnorm(n = n,
                              mu = rep(mu, p),
                              Sigma = Sigma)
  
  colnames(raw_scores) <- paste0("R", 1:p)
  
  # --- Standardization (GSVA-like scaling) ---
  scaled_scores <- scale(raw_scores)
  
  return(as.data.frame(scaled_scores))
}

#' 
#' @param X matrix of covariates.
#' @param lambda0 Constant baseline hazard in the additive hazards model.
#' @param beta Additive regression coefficients.
#' @param censor_max Maximum censoring. Independent censoring times are 
#' generated as C \sim Uniform(0, censor_max).
#' 
#' @return A data.frame containing observed time, event indicator and 
#' covariates.
#' 
#' @export
simulate_aalen_time <- function(X, lambda0, beta, censor_max){
  
  n <- nrow(X)
  
  # Events times from exponential with rate lambdaX
  lambdaX <- as.vector(lambda0 + X %*% beta) 
  U <- runif(n)
  T <- -log(U)/lambdaX
  
  # Generate final censoring times
  C <- runif(n, 0, censor_max)
  
  time <- pmin(T, C)
  event <- as.integer(T <= C)
  
  data.frame(time, event, X)
}


#' Calibrate
#' @param n_cal Number of individuals used in the Monte Carlo calibration step.
#' @param lambda0 Constant baseline hazard.
#' @param beta Additive regression coefficients.
#' @param prev Marginal prevalences of the binary covariates.
#' @param target_cens Desired censoring proportion.
#' @param seed Random seed used during the calibration procedure. 
#' @param max_upper Upper bound used during the root-finding procedure 
#' (`uniroot`) to locate the value of cmax. 
#' 
#' @return Single numerical value which corresponds to the calibrated parameter 
#' cmax such that the expected censoring proportion is approximately equal to 
#' target_cens. 
#'
#' @export
calibrate_censor_max_aalen <- function(
    n_cal = 200000,
    lambda0 = 0.05, 
    beta = c(0.06, 0.02, 0.005, 0),
    prev = c(0.005, 0.02, 0.05, 0.20),
    target_cens = 0.20, 
    seed = 2026, 
    max_upper = 1e6
    ){
  set.seed(seed)
  stopifnot(length(beta)==4, length(prev)==4)
  
  # Simulate binary covariates 
  X <- sapply(prev, function(prob) rbinom(n_cal, 1, prob))
  
  # Hazards 
  lambdaX <- as.vector(lambda0 + X %*% beta)
  
  # Event times 
  U <- runif(n_cal)
  T <- -log(U)/lambdaX
  
  # Deterministic expected censoring rate for Uniform(0, cmax)
  cens_rate_det <- function(cmax) mean(pmin(1, T / cmax))
  f <- function(cmax) cens_rate_det(cmax) - target_cens
  
  lower <- 1e-8
  upper <- 1
  while(f(upper) > 0 && upper < max_upper) upper <- upper * 2
  if (f(lower)* f(upper) > 0) stop("Could not bracket root; increase max_upper or check target_cens.")
  
  censor_max <- uniroot(f, interval = c(lower, upper))$root
  censor_max
}
  
  