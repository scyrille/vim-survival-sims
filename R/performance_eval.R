
#' Estimate empirical bias
#'
#' @param est Numeric vector of estimates
#' @param true True parameter value
#'
#' @return Numeric scalar (bias)
estimate_bias <- function(est, true){
  mean(est - true)
}

#' Standardized bias is typically defined as Bias / SD(true quantity).
#' In simulation papers for scalar parameters, a common choice is:
#'   std_bias = mean(est - true) / sd(est)
#' (or sd of the true influence / true parameter across reps if it varies).
#'
#' @param est Numeric vector of estimates
#' @param true True parameter value (scalar)
#' @param denom Character: denominator for standardization.
#'   - "sd_est" (default): sd(est)
#'   - "sd_error": sd(est - true)
#'
#' @return Numeric scalar (standardized bias)
estimate_std_bias <- function(est, true, denom = c("sd_est", "sd_error")){
  denom <- match.arg(denom)
  err <- est - true
  d <- switch(
    denom,
    sd_est   = stats::sd(est, na.rm = TRUE),
    sd_error = stats::sd(err, na.rm = TRUE)
  )
  if (is.na(d) || d == 0) return(NA_real_)
  mean(err, na.rm = TRUE) / d
}

#' Estimate root mean squared error (RMSE)
#'
#' @param est Numeric vector of estimates
#' @param true True parameter value
#'
#' @return Numeric scalar (RMSE)
estimate_rmse <- function(est, true){
  sqrt(mean((est - true)^2, na.rm = TRUE))
}


#' Estimate coverage probability
#'
#' @param lower Vector of lower CI bounds
#' @param upper Vector of upper CI bounds
#' @param true True parameter value
#'
#' @return Numeric scalar (coverage)
estimate_coverage <- function(lower, upper, true){
  mean(lower <= true & upper >= true, na.rm = TRUE)
}

#' Estimate mean confidence interval width
#'
#' @param lower Vector of lower CI bounds
#' @param upper Vector of upper CI bounds
#'
#' @return Numeric scalar (mean CI width)
estimate_ci_width <- function(lower, upper){
  mean(upper - lower, na.rm = TRUE)
}

#' Estimate type I error
#'
#' @param pval Vector of p-values
#' @param alpha Significance level (default 0.05)
#'
#' @return Numeric scalar (type I error rate)
estimate_type1_error <- function(pval, alpha = 0.05){
  mean(pval < alpha, na.rm = TRUE)
}

#' Estimate rank correlation (Spearman by default)
#'
#' Useful to assess whether an estimator preserves ordering (e.g., risk scores,
#' pathway rankings, predicted effects).
#'
#' @param x Numeric vector (e.g., predicted/ranked values)
#' @param y Numeric vector (e.g., true values)
#' @param method Correlation method: "spearman" (default) 
#'
#' @return Numeric scalar (rank correlation)
estimate_rank_correlation <- function(x, y, method = "spearman"){
  method <- match.arg(method)
  stats::cor(x, y, method = method, use = "complete.obs")
}

simulation_metrics <- function(est, lower, upper, true,
                               pval = NULL, alpha = 0.05,
                               std_bias_denom = c("sd_est", "sd_error"),
                               rank_x = NULL, rank_y = NULL,
                               rank_method = "spearman"){
  
  std_bias_denom <- match.arg(std_bias_denom)
  rank_method <- match.arg(rank_method)
  
  out <- list(
    bias      = estimate_bias(est, true),
    std_bias  = estimate_std_bias(est, true, denom = std_bias_denom),
    rmse      = estimate_rmse(est, true),
    coverage  = estimate_coverage(lower, upper, true),
    ci_width  = estimate_ci_width(lower, upper)
  )
  
  if(!is.null(pval)){
    out$type1_or_power <- mean(pval < alpha, na.rm = TRUE)
  }
  
  if(!is.null(rank_x) && !is.null(rank_y)){
    out$rank_cor <- estimate_rank_correlation(rank_x, rank_y, method = rank_method)
  }
  
  out
}