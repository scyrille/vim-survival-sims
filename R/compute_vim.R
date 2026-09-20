
compute_vim <- function(time, 
                        event, 
                        X, 
                        indxs, 
                        tau, 
                        approx_times, 
                        cf_folds, 
                        ss_folds, 
                        nuisance, 
                        scale_est = FALSE){
  
  start <- Sys.time()
  
  vims <- c("BS(t)","AUC(t)")
  
  message(paste0("Estimating conditional survival functions by ",
                 nuisance, 
                 "..."))
  
  V0_preds <- CV_generate_full_predictions_landmark(
    time = time,
    event = event,
    X = X,
    tau = tau,
    approx_times = approx_times,
    nuisance = nuisance,
    cf_folds = cf_folds
  )
  
  CV_full_preds <- V0_preds$CV_full_preds
  CV_full_preds_train <- V0_preds$CV_full_preds_train
  CV_S_preds <- V0_preds$CV_S_preds
  CV_G_preds <- V0_preds$CV_G_preds

  output <- purrr::map_dfr(indxs, function(indx_i) {
    char_indx <- as.character(indx_i)
    indx <- as.numeric(strsplit(char_indx, split = ",")[[1]])
    variable <- names(X)[indx]
    
    message(paste0("Estimating the importance of ",
                   variable,
                   " with nuisance functions estimated by ",
                   nuisance, 
                   "..."))
    
    CV_reduced_preds <- CV_generate_reduced_predictions_landmark(
      time = time,
      event = event,
      X = X,
      tau = tau,
      cf_folds = cf_folds,
      indx = indx,
      full_preds_train = CV_full_preds_train
    )
    
    purrr::map_dfr(vims, function(vim) {
      output <- switch(
        vim,
        `BS(t)` = survML::vim_brier(
          time = time,
          event = event,
          approx_times = approx_times,
          landmark_times = tau,
          f_hat = CV_full_preds,
          fs_hat = CV_reduced_preds,
          S_hat = CV_S_preds,
          G_hat = CV_G_preds,
          sample_split = TRUE,
          cf_folds = cf_folds,
          ss_folds = ss_folds,
          scale_est = scale_est
        ), 
        `AUC(t)` = survML::vim_AUC(
          time = time,
          event = event,
          approx_times = approx_times,
          landmark_times = tau,
          f_hat = purrr::map(CV_full_preds, ~ 1 - .x),
          fs_hat = purrr::map(CV_reduced_preds, ~ 1 - .x),
          S_hat = CV_S_preds,
          G_hat = CV_G_preds,
          sample_split = TRUE,
          cf_folds = cf_folds,
          ss_folds = ss_folds,
          scale_est = scale_est
        )
      )
      output %>%
        dplyr::mutate(
          vim = vim,
          variable = variable
        )
    })
    
  })
  
  end <- Sys.time()
  runtime <- as.numeric(difftime(end, start, units = "mins"))
  
  output %>%
    dplyr::rename(., tau = landmark_time)%>%
    dplyr::mutate(runtime = runtime, 
                  nuisance = nuisance) 
}
