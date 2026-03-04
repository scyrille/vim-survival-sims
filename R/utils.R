#' Pairwise Fisher's exact tests 
#'
#' @param data data.frame containing binary variables
#' @param vars character vector of variable names
#'
#' @return A list with:
#'   - p.value.matrix : symmetric matrix of raw p-values
#'   - tidy.results   : long-format tibble of pairwise results
#'   
#' @export
pairwise_fisher <- function(data, vars) {
  
  stopifnot(all(vars %in% names(data)))
  
  combn_pairs <- combn(vars, 2, simplify = FALSE)
  
  results <- purrr::map_dfr(combn_pairs, function(pair) {
    
    x <- data[[pair[1]]]
    y <- data[[pair[2]]]
    
    tab <- table(x, y)
    
    test <- fisher.test(tab)
    
    tibble::tibble(
      var1 = pair[1],
      var2 = pair[2],
      p.value = test$p.value,
      odds.ratio = unname(test$estimate)
    )
  })
  
  # Create symmetric matrix of raw p-values
  mat <- matrix(NA, length(vars), length(vars),
                dimnames = list(vars, vars))
  
  for(i in seq_len(nrow(results))) {
    mat[results$var1[i], results$var2[i]] <- results$p.value[i]
    mat[results$var2[i], results$var1[i]] <- results$p.value[i]
  }
  
  diag(mat) <- 0
  
  list(
    p.value.matrix = mat,
    tidy.results = results
  )
}