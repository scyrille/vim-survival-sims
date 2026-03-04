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
  as.data.frame(X)
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

