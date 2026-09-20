
do_one <- function(n, scenario, c_max, tau, nuisances, scale_est = FALSE){
  
  start <- Sys.time()
  
  vims <- c("BS(t)","AUC(t)")
  
  message("Generating data...")
  
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
  
  output <- nuisances %>%
    purrr::map(
      ~compute_vim(
        time = time, 
        event = event, 
        X = X, 
        indxs = indxs, 
        tau = tau, 
        approx_times = approx_times, 
        cf_folds = cf_folds, 
        ss_folds = ss_folds, 
        nuisance = .x, 
        scale_est = scale_est
      )
    )%>% 
    bind_rows() %>%
    dplyr::mutate(scenario = scenario, 
                  n = n)
  
  output 
}
