
do_one <- function(n, scenario, c_max, tau, nuisance){
  
  start <- Sys.time()
  
  vims <- c("AUC","brier")
  
  if (scenario == 1){
    
    get_scenario_param("1")
    data <- generate_data(n = n, scenario = "1", c_max = c_max)
    
  } else if (scenario == "2"){
    
    get_scenario_param("2")
    data <- generate_data(n = n, scenario = "2", c_max = c_max)
    
  }
  
  sample_split <- TRUE
  V <- 5 
  
  time <- data$time
  event <- data$event
  X <- data %>% dplyr::select(-c(time,event))
  indxs <- paste0(1:ncol(X))
  
  approx_times <- sort(unique(c(0, time[event == 1], tau)))
  approx_times <- approx_times[approx_times <= max(tau)]
  
  folds <- survML:::generate_folds(n, V, sample_split)
  cf_folds <- folds$cf_folds
  ss_folds <- folds$ss_folds
  
  V0_preds <- CV_generate_full_predictions_landmark(
    time = time,
    event = event,
    X = X,
    landmark_times = tau,
    approx_times = approx_times,
    nuisance = nuisance,
    cf_folds = cf_folds)
  
  CV_full_preds <- V0_preds$CV_full_preds
  CV_full_preds_train <- V0_preds$CV_full_preds_train
  CV_S_preds <- V0_preds$CV_S_preds
  CV_G_preds <- V0_preds$CV_G_preds
  
  shared_settings <- expand.grid(indx = indxs, vim = vims) 
  
  output <- purrr::map_dfr(indxs, function(indx_i) {
    char_indx <- as.character(indx_i)
    indx <- as.numeric(strsplit(char_indx, split = ",")[[1]])
    variable <- names(X)[indx]
    
    CV_reduced_preds <- CV_generate_reduced_predictions_landmark(
      time = time,
      event = event,
      X = X,
      landmark_times = tau,
      cf_folds = cf_folds,
      indx = indx,
      full_preds_train = CV_full_preds_train
    )
    
    purrr::map_dfr(vims, function(vim) {
      output <- switch(
        vim,
        brier = survML::vim_brier(
          time = time,
          event = event,
          approx_times = approx_times,
          landmark_times = tau,
          f_hat = CV_full_preds,
          fs_hat = CV_reduced_preds,
          S_hat = CV_S_preds,
          G_hat = CV_G_preds,
          cf_folds = cf_folds,
          sample_split = sample_split,
          ss_folds = ss_folds#,
          # scale_est = TRUE
        ), 
        AUC = survML::vim_AUC(
          time = time,
          event = event,
          approx_times = approx_times,
          landmark_times = tau,
          f_hat = purrr::map(CV_full_preds, ~ 1 - .x),
          fs_hat = purrr::map(CV_reduced_preds, ~ 1 - .x),
          S_hat = CV_S_preds,
          G_hat = CV_G_preds,
          cf_folds = cf_folds,
          sample_split = sample_split,
          ss_folds = ss_folds#,
          # scale_est = TRUE
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
  output <- output %>%
    dplyr::mutate(scenario = scenario, 
                  runtime = runtime, 
                  n = n, 
                  nuisance = nuisance)%>%
    dplyr::relocate(vim, variable)
  return(output)
}
