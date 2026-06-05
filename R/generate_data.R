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
simulate_binary <- function(n, 
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
  
  colnames(X) <- paste0("X", seq_len(p))
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
simulate_continuous <- function(n,
                                p = 20,
                                mu = 0, 
                                Sigma) {
  
  # --- Multivariate normal simulation ---
  raw_scores <- MASS::mvrnorm(n = n,
                              mu = rep(mu, p),
                              Sigma = Sigma)
  
  colnames(raw_scores) <- paste0("X", 1:p)
  
  # --- Standardization (GSVA-like scaling) ---
  scaled_scores <- scale(raw_scores)
  
  return(as.data.frame(scaled_scores))
}


#' Simulate data under a Cox proportional hazards model
#'
#' @param n Integer. Sample size.
#' @param prev Numeric vector of length `p` with marginal prevalences of binary
#'   covariates.
#' @param sigma Numeric `p x p` correlation/covariance matrix for latent Gaussian variables.
#' @param lambda0 Numeric. Constant baseline hazard.
#' @param beta Numeric vector of log-hazard ratios.
#' @param c_max Numeric. Maximum censoring time. C ~ Uniform(0, c_max).
#'
#' @return A data.frame containing observed time, event indicator and covariates.
#'
#' @export
simulate_cox_data <- function(n, prev, sigma, lambda0, beta, c_max) {
  
  p <- length(prev)
  
  stopifnot(length(beta) == p)
  stopifnot(is.numeric(prev), all(prev > 0 & prev < 1))
  stopifnot(is.matrix(sigma), all(dim(sigma) == c(p, p)))
  stopifnot(lambda0 > 0)
  stopifnot(c_max > 0)
  
  # Latent Gaussian generation
  Z <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = sigma)
  
  # Thresholds matching marginal prevalences
  thresholds <- qnorm(prev)
  
  # Dichotomization: P(X_j = 1) = prev_j
  X <- sweep(Z, 2, thresholds, "<") * 1L
  X <- as.matrix(X)
  colnames(X) <- paste0("X", seq_len(p))
  
  # Individual hazard rate under Cox PH model
  # lambda(t | X) = lambda0 * exp(X beta)
  linpred <- as.vector(X %*% beta)
  lambdaX <- lambda0 * exp(linpred)
  
  # Event times from exponential distribution
  U <- runif(n)
  T <- -log(U) / lambdaX
  
  # Independent censoring times
  C <- runif(n, min = 0, max = c_max)
  
  # Observed data
  time <- pmin(T, C)
  event <- as.integer(T <= C)
  
  data.frame(time = time, event = event, X, check.names = FALSE)
}

#' Simulate data under an Aalen additive model
#'
#' @param n Integer. Sample size.
#' @param prev Numeric vector of length `p` with marginal prevalences of binary
#'   covariates (values in (0,1)).
#' @param sigma Numeric `p x p` correlation (or covariance) matrix for the
#'   latent Gaussian variables.
#' @param lambda0 Numeric. Constant baseline hazard in the additive hazards model.
#' @param alpha Numeric vector of additive regression coefficients.
#' @param c_max Numeric. Maximum censoring time. Independent censoring
#'   times are generated as C ~ Uniform(0, c_max).
#'
#' @return A data.frame containing observed time, event indicator and covariates.
#'
#' @export
simulate_aalen_data <- function(n, prev, sigma, lambda0, alpha, c_max) {
  
  p <- length(prev)
  
  stopifnot(length(alpha) == p)
  stopifnot(is.numeric(prev), all(prev > 0 & prev < 1))
  stopifnot(is.matrix(sigma), all(dim(sigma) == c(p, p)))
  stopifnot(c_max > 0)
  
  # Theoretical positivity check for binary covariates
  min_lambda <- lambda0 + sum(alpha[alpha < 0])
  if (min_lambda <= 0) {
    stop(
      "The additive hazard can be non-positive for some covariate patterns. ",
      "Increase lambda0 or modify alpha."
    )
  }
  
  # Latent Gaussian generation
  Z <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = sigma)
  
  # Thresholds matching marginal prevalences
  thresholds <- qnorm(prev)
  
  # Dichotomization
  X <- sweep(Z, 2, thresholds, "<") * 1L
  X <- as.matrix(X)
  colnames(X) <- paste0("X", seq_len(p))
  
  # Event times from exponential with individual rate lambdaX
  lambdaX <- as.vector(lambda0 + X %*% alpha)
  
  if (any(lambdaX <= 0)) {
    stop(
      "Some sampled values of the additive hazard are <= 0. ",
      "Increase lambda0 or modify alpha."
    )
  }
  
  U <- runif(n)
  T <- -log(U) / lambdaX
  
  # Generate censoring times
  C <- runif(n, min = 0, max = c_max)
  
  # Observed data
  time <- pmin(T, C)
  event <- as.integer(T <= C)
  
  data.frame(time = time, event = event, X, check.names = FALSE)
}


#' Simulate data under a Cox-Aalen model
#'
#' @param n Integer. Sample size.
#' @param lambda0 Constant baseline hazard in the additive part.
#' @param alpha Numeric vector of additive regression coefficients.
#' @param beta Numeric vector of multiplicative regression coefficients.
#' @param prev_add Numeric vector with marginal prevalences of additive binary covariates.
#' @param sigma_add Numeric correlation/covariance matrix for the latent Gaussian
#'   variables used to generate additive binary covariates.
#' @param prev_mult Numeric vector with marginal prevalences of multiplicative binary covariates.
#' @param mean_mult Numeric vector of means for multiplicative continuous covariates.
#' @param sigma_mult_bin Numeric correlation/covariance matrix for the latent Gaussian
#'   variables used to generate multiplicative binary covariates.
#' @param sigma_mult_cont Numeric covariance matrix for multiplicative continuous covariates.
#' @param c_max Maximum censoring time. Independent censoring times are
#'   generated as C ~ Uniform(0, censor_max).
#'
#' @return A data.frame containing observed time, event indicator, and covariates.
#'
#' @export
simulate_cox_aalen_data <- function(n,
                                    lambda0,
                                    alpha,
                                    beta,
                                    prev_add,
                                    sigma_add,
                                    prev_mult,
                                    mean_mult,
                                    sigma_mult_bin,
                                    sigma_mult_cont,
                                    c_max) {
  
  # ---------- Dimensions ---------- #
  p_add <- length(prev_add)
  p_mult_bin <- length(prev_mult)
  p_mult_cont <- length(mean_mult)
  p_mult <- p_mult_bin + p_mult_cont
  
  stopifnot(length(alpha) == p_add)
  stopifnot(length(beta) == p_mult)
  stopifnot(is.matrix(sigma_add), all(dim(sigma_add) == c(p_add, p_add)))
  stopifnot(is.matrix(sigma_mult_bin), all(dim(sigma_mult_bin) == c(p_mult_bin, p_mult_bin)))
  stopifnot(is.matrix(sigma_mult_cont), all(dim(sigma_mult_cont) == c(p_mult_cont, p_mult_cont)))
  stopifnot(c_max > 0)
  
  # --------- Additive covariates (binary) --------- #
  Z_add <- MASS::mvrnorm(n = n, mu = rep(0, p_add), Sigma = sigma_add)
  thresholds_add <- qnorm(prev_add)
  X_add <- sweep(Z_add, 2, thresholds_add, "<") * 1L
  X_add <- as.matrix(X_add)
  colnames(X_add) <- paste0("X", seq_len(p_add))
  
  # --------- Multiplicative covariates --------- #
  ## Binary
  Z_mult_bin <- MASS::mvrnorm(
    n = n,
    mu = rep(0, p_mult_bin),
    Sigma = sigma_mult_bin
  )
  thresholds_mult <- qnorm(prev_mult)
  X_mult_bin <- sweep(Z_mult_bin, 2, thresholds_mult, "<") * 1L
  X_mult_bin <- as.matrix(X_mult_bin)
  
  ## Continuous
  X_mult_cont <- scale(
    MASS::mvrnorm(n = n, mu = mean_mult, Sigma = sigma_mult_cont)
  )
  X_mult_cont <- as.matrix(X_mult_cont)
  
  ## Combine multiplicative covariates
  X_mult <- cbind(X_mult_bin, X_mult_cont)
  X_mult <- as.matrix(X_mult)
  colnames(X_mult) <- paste0("Z", seq_len(p_mult))
  
  # --------- Hazard components --------- #
  add_part <- as.vector(lambda0 + X_add %*% alpha)
  
  if (any(add_part <= 0)) {
    stop("Some values of the additive part are <= 0. Increase lambda0 or modify alpha.")
  }
  
  lambda_i <- as.vector(exp(X_mult %*% beta) * add_part)
  
  # --------- Event times --------- #
  U <- runif(n)
  T <- -log(U) / lambda_i
  
  # --------- Censoring times --------- #
  C <- runif(n, min = 0, max = c_max)
  
  # --------- Observed data --------- #
  time <- pmin(T, C)
  event <- as.integer(T <= C)
  
  data.frame(time = time, event = event, X_add, X_mult, check.names = FALSE)
}

#' Calibrate the upper bound of a Uniform censoring distribution
#' to achieve a target censoring proportion under a Cox proportional hazards model.
#'
#' @param n_cal Integer. Number of individuals used in the Monte Carlo calibration step.
#' @param lambda0 Numeric. Constant baseline hazard.
#' @param beta Numeric vector of log-hazard ratios.
#' @param prev Numeric vector of marginal prevalences of the binary covariates.
#' @param sigma Numeric correlation/covariance matrix for the latent Gaussian variables.
#' @param target_cens Numeric. Desired censoring proportion.
#' @param seed Integer. Random seed used during the calibration procedure.
#' @param max_upper Numeric. Upper bound used during the root-finding procedure.
#'
#' @return A single numeric value corresponding to the calibrated parameter
#'   `c_max` such that the expected censoring proportion is approximately
#'   equal to `target_cens`.
#'
#' @export
calibrate_c_max_cox <- function(
    n_cal = 200000,
    lambda0 = 0.1,
    beta,
    prev,
    sigma,
    target_cens = 0.20,
    seed = 2026,
    max_upper = 1e6
) {
  set.seed(seed)
  
  p <- length(prev)
  
  stopifnot(length(beta) == p)
  stopifnot(is.numeric(prev), all(prev > 0 & prev < 1))
  stopifnot(is.matrix(sigma), all(dim(sigma) == c(p, p)))
  stopifnot(lambda0 > 0)
  stopifnot(target_cens > 0, target_cens < 1)
  
  # ------ Simulate binary covariates ----- #
  Z <- MASS::mvrnorm(n = n_cal, mu = rep(0, p), Sigma = sigma)
  
  # Thresholds matching marginal prevalences
  thresholds <- qnorm(prev)
  
  # Dichotomization: P(X_j = 1) = prev_j
  X <- sweep(Z, 2, thresholds, "<") * 1L
  X <- as.matrix(X)
  colnames(X) <- paste0("X", seq_len(p))
  
  # ------ Individual hazards under Cox PH ----- #
  # lambda(t | X) = lambda0 * exp(X beta)
  linpred <- as.vector(X %*% beta)
  lambdaX <- lambda0 * exp(linpred)
  
  # ------ Event times ----- #
  U <- runif(n_cal)
  T <- -log(U) / lambdaX
  
  # ------ Expected censoring proportion under C ~ Uniform(0, c_max) ----- #
  cens_rate_det <- function(cmax) {
    mean(pmin(1, T / cmax))
  }
  
  f <- function(cmax) cens_rate_det(cmax) - target_cens
  
  lower <- 1e-8
  upper <- 1
  
  while (f(upper) > 0 && upper < max_upper) {
    upper <- upper * 2
  }
  
  if (f(lower) * f(upper) > 0) {
    stop("Could not bracket root; increase max_upper or check target_cens.")
  }
  
  censor_max <- uniroot(f, interval = c(lower, upper))$root
  
  censor_max
}


#' Calibrate the upper bound of a Uniform censoring distribution
#' to achieve a target censoring proportion under an Aalen additive model.
#'
#' @param n_cal Integer. Number of individuals used in the Monte Carlo calibration step.
#' @param lambda0 Numeric. Constant baseline hazard.
#' @param alpha Numeric vector of additive regression coefficients.
#' @param prev Numeric vector of marginal prevalences of the binary covariates.
#' @param Sigma Numeric correlation/covariance matrix for the latent Gaussian variables.
#' @param target_cens Numeric. Desired censoring proportion.
#' @param seed Integer. Random seed used during the calibration procedure.
#' @param max_upper Numeric. Upper bound used during the root-finding procedure
#'   (`uniroot`) to locate the value of `censor_max`.
#'
#' @return A single numeric value corresponding to the calibrated parameter
#'   `c_max` such that the expected censoring proportion is approximately
#'   equal to `target_cens`.
#'
#' @export
calibrate_c_max_aalen <- function(
    n_cal = 200000,
    lambda0 = 0.05,
    alpha = c(0.06, 0.02, 0.005, 0),
    prev = c(0.005, 0.02, 0.05, 0.20),
    sigma,
    target_cens = 0.20,
    seed = 2026,
    max_upper = 1e6
) {
  set.seed(seed)
  
  p <- length(prev)
  
  stopifnot(length(alpha) == p)
  stopifnot(is.numeric(prev), all(prev > 0 & prev < 1))
  stopifnot(is.matrix(sigma), all(dim(sigma) == c(p, p)))
  stopifnot(target_cens > 0, target_cens < 1)
  
  # Theoretical positivity check for binary covariates
  min_lambda <- lambda0 + sum(alpha[alpha < 0])
  if (min_lambda <= 0) {
    stop(
      "The additive hazard can be non-positive for some covariate patterns. ",
      "Increase lambda0 or modify alpha."
    )
  }
  
  # ------ Simulate binary covariates ----- #
  Z <- MASS::mvrnorm(n = n_cal, mu = rep(0, p), Sigma = sigma)
  
  # Thresholds matching marginal prevalences
  thresholds <- qnorm(prev)
  
  # Dichotomization
  X <- sweep(Z, 2, thresholds, "<") * 1L
  X <- as.matrix(X)
  colnames(X) <- paste0("X", seq_len(p))
  
  # ------ Individual hazards ----- #
  lambdaX <- as.vector(lambda0 + X %*% alpha)
  
  if (any(lambdaX <= 0)) {
    stop(
      "Some sampled values of the additive hazard are <= 0. ",
      "Increase lambda0 or modify alpha."
    )
  }
  
  # ------ Event times ----- #
  U <- runif(n_cal)
  T <- -log(U) / lambdaX
  
  # ------ Expected censoring proportion under C ~ Uniform(0, cmax) ----- #
  cens_rate_det <- function(cmax) {
    mean(pmin(1, T / cmax))
  }
  
  f <- function(cmax) cens_rate_det(cmax) - target_cens
  
  lower <- 1e-8
  upper <- 1
  
  while (f(upper) > 0 && upper < max_upper) {
    upper <- upper * 2
  }
  
  if (f(lower) * f(upper) > 0) {
    stop("Could not bracket root; increase max_upper or check target_cens.")
  }
  
  censor_max <- uniroot(f, interval = c(lower, upper))$root
  censor_max
}

#' Calibrate the upper bound of a Uniform censoring distribution
#' to achieve a target censoring proportion under a Cox-Aalen model.
#'
#' @param n_cal Integer. Number of individuals used in the Monte Carlo calibration step.
#' @param lambda0 Constant baseline hazard in the additive part.
#' @param alpha Numeric vector of additive regression coefficients.
#' @param beta Numeric vector of multiplicative regression coefficients.
#' @param prev_add Numeric vector with marginal prevalences of additive binary covariates.
#' @param sigma_add Numeric correlation/covariance matrix for the latent Gaussian
#'   variables used to generate additive binary covariates.
#' @param prev_mult Numeric vector with marginal prevalences of multiplicative binary covariates.
#' @param mean_mult Numeric vector of means for multiplicative continuous covariates.
#' @param sigma_mult_bin Numeric correlation/covariance matrix for the latent Gaussian
#'   variables used to generate multiplicative binary covariates.
#' @param sigma_mult_cont Numeric covariance matrix for multiplicative continuous covariates.
#' @param target_cens Target censoring proportion.
#' @param seed Integer. Random seed.
#' @param max_upper Numeric. Maximum upper bound used in the root-finding procedure.
#'
#' @return A single numeric value corresponding to the calibrated c_max.
#'
#' @export
calibrate_c_max_cox_aalen <- function(
    n_cal = 200000,
    lambda0,
    alpha,
    beta,
    prev_add,
    sigma_add,
    prev_mult,
    mean_mult,
    sigma_mult_bin,
    sigma_mult_cont,
    target_cens,
    seed,
    max_upper = 1e6
) {
  set.seed(seed)
  
  # ---------- Dimensions ---------- #
  p_add <- length(prev_add)
  p_mult_bin <- length(prev_mult)
  p_mult_cont <- length(mean_mult)
  p_mult <- p_mult_bin + p_mult_cont
  
  stopifnot(length(alpha) == p_add)
  stopifnot(length(beta) == p_mult)
  stopifnot(is.matrix(sigma_add), all(dim(sigma_add) == c(p_add, p_add)))
  stopifnot(is.matrix(sigma_mult_bin), all(dim(sigma_mult_bin) == c(p_mult_bin, p_mult_bin)))
  stopifnot(is.matrix(sigma_mult_cont), all(dim(sigma_mult_cont) == c(p_mult_cont, p_mult_cont)))
  stopifnot(target_cens > 0, target_cens < 1)
  
  # ---------- Quick theoretical positivity check ---------- #
  # Since X_add is binary, the minimum possible additive part is obtained
  # when all covariates with negative coefficients are equal to 1.
  min_add_theoretical <- lambda0 + sum(alpha[alpha < 0])
  
  if (min_add_theoretical <= 0) {
    stop(
      "The additive part can be non-positive for some binary covariate patterns. ",
      "Increase lambda0 or modify alpha."
    )
  }
  
  # ---------- Additive covariates (binary) ---------- #
  Z_add <- MASS::mvrnorm(
    n = n_cal,
    mu = rep(0, p_add),
    Sigma = sigma_add
  )
  
  thresholds_add <- qnorm(prev_add)
  X_add <- sweep(Z_add, 2, thresholds_add, "<") * 1L
  X_add <- as.matrix(X_add)
  
  # ---------- Multiplicative covariates ---------- #
  ## Binary
  Z_mult_bin <- MASS::mvrnorm(
    n = n_cal,
    mu = rep(0, p_mult_bin),
    Sigma = sigma_mult_bin
  )
  
  thresholds_mult <- qnorm(prev_mult)
  X_mult_bin <- sweep(Z_mult_bin, 2, thresholds_mult, "<") * 1L
  X_mult_bin <- as.matrix(X_mult_bin)
  
  ## Continuous
  X_mult_cont <- scale(
    MASS::mvrnorm(
      n = n_cal,
      mu = mean_mult,
      Sigma = sigma_mult_cont
    )
  )
  X_mult_cont <- as.matrix(X_mult_cont)
  
  ## All multiplicative covariates
  X_mult <- cbind(X_mult_bin, X_mult_cont)
  X_mult <- as.matrix(X_mult)
  
  # ---------- Individual hazards ---------- #
  add_part <- as.vector(lambda0 + X_add %*% alpha)
  
  if (any(add_part <= 0)) {
    stop(
      "Some sampled values of the additive part are <= 0. ",
      "Increase lambda0 or modify alpha."
    )
  }
  
  lambda_i <- as.vector(exp(X_mult %*% beta) * add_part)
  
  # ---------- Event times ---------- #
  U <- runif(n_cal)
  T <- -log(U) / lambda_i
  
  # ---------- Expected censoring under C ~ Uniform(0, cmax) ---------- #
  cens_rate_det <- function(cmax) {
    mean(pmin(1, T / cmax))
  }
  
  f <- function(cmax) cens_rate_det(cmax) - target_cens
  
  lower <- 1e-8
  upper <- 1
  
  while (f(upper) > 0 && upper < max_upper) {
    upper <- upper * 2
  }
  
  if (f(lower) * f(upper) > 0) {
    stop("Could not bracket root; increase max_upper or check target_cens.")
  }
  
  uniroot(f, interval = c(lower, upper))$root
}

