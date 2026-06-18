
get_scenario_param <- function(scenario){
  
  if (scenario == "1"){
    
    p           = 4
    lambda0     = 0.05 
    alpha       = c(0.1,    # strong
                    0.05,   # moderate
                    0.005,  # weak
                    0)      # null, 
    sigma       = diag(x = 1, nrow = p, ncol = p)
    prev        = rep(0.4, p)
    target_cens = 0.2
    
    params = list(p           = p,
                  lambda0     = lambda0, 
                  alpha       = alpha,
                  sigma       = sigma, 
                  prev        = prev,
                  target_cens = target_cens)
    
    list2env(params, .GlobalEnv)
    
  } else if (scenario == "2"){
    
    p <- 6
    prev_add <- c(0.4, 0.4, 0.4)
    sigma_add <- diag(x = 1, nrow = 3, ncol = 3)
    alpha <- c(0.1,  # strong
               0.05, # intermediate
               0)    # null
    prev_mult <- c(0.4, 0.4)
    sigma_mult_bin <- diag(x = 1, nrow = 2, ncol = 2)
    mean_mult <- 0
    sigma_mult_cont <- diag(x = 1, nrow = 1, ncol = 1)
    beta <- c(1,   # strong
              0.5, # intermediate
              0)   # null
    lambda0 <- 0.05
    target_cens <- 0.2
    
    params = list(p               = p,
                  lambda0         = lambda0,
                  alpha           = alpha,
                  beta            = beta,
                  prev_add        = prev_add,
                  sigma_add       = sigma_add,
                  prev_mult       = prev_mult,
                  mean_mult       = mean_mult,
                  sigma_mult_bin  = sigma_mult_bin,
                  sigma_mult_cont = sigma_mult_cont,
                  target_cens     = target_cens)
    
    list2env(params, .GlobalEnv)
  }
}


generate_data <- function(n, scenario, c_max){
  
  if (scenario == "1"){
    
    get_scenario_param("1")
    
    p <- length(prev)
    
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
    
    # ------ Event times ----- #
    haz <- as.numeric(lambda0 + X %*% alpha)
    T <- rexp(n, rate = haz)
    
    # --------- Censoring times --------- #
    C <- runif(n, min = 0, max = c_max)
    
    # --------- Observed data --------- #
    time <- pmin(T, C)
    event <- as.integer(T <= C)
    
    data <- data.frame(time = time, event = event, X, 
                       check.names = FALSE)
  
  } else if (scenario == "2"){
    
    get_scenario_param("2")
    
    # ---------- Dimensions ---------- #
    p_add <- length(prev_add)
    p_mult_bin <- length(prev_mult)
    p_mult_cont <- length(mean_mult)
    p_mult <- p_mult_bin + p_mult_cont
    
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
    
    # --------- Event times --------- #
    haz_add <- as.numeric(lambda0 + X_add %*% alpha)
    
    if (any(haz_add < 0))
      warning(sum(haz_add < 0),
              " subjects with negative additive hazard.")
    
    haz <- pmax(haz_add, 0) * exp(as.numeric(X_mult %*% beta))
    T <- rexp(n, rate = haz)
    
    # --------- Censoring times --------- #
    C <- runif(n, min = 0, max = c_max)
    
    # --------- Observed data --------- #
    time <- pmin(T, C)
    event <- as.integer(T <= C)
    
    data <- data.frame(time = time, event = event, X_add, X_mult, 
                       check.names = FALSE)
    
  }
  return(data)
}


calibrate_parameters <- function(n, scenario, max_upper, probs){
  
  if (scenario == 1){
    
    get_scenario_param("1")
    
    # ------ Simulate binary covariates ----- #
    Z <- MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = sigma)
    
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
    haz <- as.numeric(lambda0 + X %*% alpha)
    T <- rexp(n, rate = haz)
    
    # ------ Expected censoring proportion under C ~ Uniform(0, cmax) ----- #
    cens_rate_det_1 <- function(cmax) {
      mean(pmin(1, T / cmax))
    }
    
    f <- function(cmax) cens_rate_det_1(cmax) - target_cens
    
    lower <- 1e-8
    upper <- 1
    
    while (f(upper) > 0 && upper < max_upper) {
      upper <- upper * 2
    }
    
    if (f(lower) * f(upper) > 0) {
      stop("Could not bracket root; increase max_upper or check target_cens.")
    }
  
    c_max <- uniroot(f, interval = c(lower, upper))$root
    
    #------------- Time horizon -------------#
    # 75th percentile of the true event-time distribution 
    tau <- as.numeric(quantile(T, probs = probs))
    
  } else if (scenario == "2"){
    
    get_scenario_param("2")
    
    # ---------- Dimensions ---------- #
    p_add <- length(prev_add)
    p_mult_bin <- length(prev_mult)
    p_mult_cont <- length(mean_mult)
    p_mult <- p_mult_bin + p_mult_cont
    
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
      n = n,
      mu = rep(0, p_add),
      Sigma = sigma_add
    )
    
    thresholds_add <- qnorm(prev_add)
    X_add <- sweep(Z_add, 2, thresholds_add, "<") * 1L
    X_add <- as.matrix(X_add)
    
    # ---------- Multiplicative covariates ---------- #
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
      MASS::mvrnorm(
        n = n,
        mu = mean_mult,
        Sigma = sigma_mult_cont
      )
    )
    X_mult_cont <- as.matrix(X_mult_cont)
    
    ## All multiplicative covariates
    X_mult <- cbind(X_mult_bin, X_mult_cont)
    X_mult <- as.matrix(X_mult)
    
    # ---------- Individual hazards ---------- #
    haz_add <- as.vector(lambda0 + X_add %*% alpha)
    
    if (any(haz_add <= 0)) {
      stop(
        "Some sampled values of the additive part are <= 0. ",
        "Increase lambda0 or modify alpha."
      )
    }
    
    # --------- Event times --------- #
    haz <- pmax(haz_add, 0) * exp(as.numeric(X_mult %*% beta))
    T <- rexp(n, rate = haz)
    
    # ---------- Expected censoring under C ~ Uniform(0, c_max) ---------- #
    cens_rate_det_2 <- function(cmax) {
      mean(pmin(1, T / cmax))
    }
    
    f <- function(cmax) cens_rate_det_2(cmax) - target_cens
    
    lower <- 1e-8
    upper <- 1
    
    while (f(upper) > 0 && upper < max_upper) {
      upper <- upper * 2
    }
    
    if (f(lower) * f(upper) > 0) {
      stop("Could not bracket root; increase max_upper or check target_cens.")
    }
    
    c_max <- uniroot(f, interval = c(lower, upper))$root
    
    #------------- Time horizon -------------#
    # 75th percentile of the true event-time distribution 
    tau <- as.numeric(quantile(T, probs = probs))
  }
  
  return(list(c_max = c_max,
              tau   = tau))
  
}
