
generate_full_predictions <- function(time,
                                      event,
                                      X,
                                      X_holdout, 
                                      landmark_times,
                                      approx_times,
                                      nuisance){
  
  if (nuisance == "stackG"){
    
    SL.library <- c("SL.mean", "SL.gam","SL.ranger")
    bin_size <- 0.5
    surv_out <- survML::stackG(time = time,
                               event = event,
                               X = X,
                               newX = rbind(X_holdout, X),
                               newtimes = approx_times,
                               time_grid_approx = approx_times,
                               bin_size = bin_size,
                               time_basis = "continuous",
                               surv_form = "PI",
                               SL_control = list(SL.library = SL.library,
                                                 V = 5))
    S_hat <- surv_out$S_T_preds[1:nrow(X_holdout),]
    G_hat <- surv_out$S_C_preds[1:nrow(X_holdout),]
    f_hat <- S_hat[,which(approx_times %in% landmark_times),drop=FALSE]
    S_hat_train <- surv_out$S_T_preds[(nrow(X_holdout)+1):(nrow(X_holdout)+nrow(X)),]
    G_hat_train <- surv_out$S_C_preds[(nrow(X_holdout)+1):(nrow(X_holdout)+nrow(X)),]
    f_hat_train <- S_hat_train[,which(approx_times %in% landmark_times),drop=FALSE]
  } 
  
  else if (nuisance == "aalen"){
    
    # Event model: S(t|X) 
    datS <- data.frame(time = time, event = event, X)
    S_fit <- timereg::aalen(survival::Surv(time, event) ~ ., data = datS)
    
    # Censoring model: G(t|X) 
    cens_event <- 1 - event
    datG <- data.frame(time = time, cens_event = cens_event, X)
    G_fit <- timereg::aalen(survival::Surv(time, cens_event) ~ ., data = datG)
    
    S_hat       <- pec::predictSurvProb(S_fit, newdata = X_holdout, times = approx_times)
    G_hat       <- pec::predictSurvProb(G_fit, newdata = X_holdout, times = approx_times)
    f_hat       <- S_hat[,which(approx_times %in% landmark_times),drop=FALSE]
    S_hat_train <- pec::predictSurvProb(S_fit, newdata = X,         times = approx_times)
    G_hat_train <- pec::predictSurvProb(G_fit, newdata = X,         times = approx_times)
    f_hat_train <- S_hat_train[,which(approx_times %in% landmark_times),drop=FALSE]
  
  }
  
  else if (nuisance == "cox.aalen"){
    
    x_vars <- colnames(X)
    prop_vars <- x_vars[grepl("Z",x_vars)]
    add_vars  <- x_vars[grepl("X",x_vars)]
    
    rhs <- c(
      if (length(prop_vars) > 0) paste0("prop(", prop_vars, ")"),
      add_vars
    )
    
    formS <- as.formula(
      paste("survival::Surv(time, event) ~", paste(rhs, collapse = " + "))
    )
    
    formG <- as.formula(
      paste("survival::Surv(time, cens_event) ~", paste(rhs, collapse = " + "))
    )
    
    # Event model: S(t|X) 
    datS <- data.frame(time = time, event = event, X)
    S_fit <- timereg::cox.aalen(formS, data = datS)
    
    # Censoring model: G(t|X) 
    cens_event <- 1 - event
    datG <- data.frame(time = time, cens_event = cens_event, X)
    G_fit       <- timereg::cox.aalen(formG, data = datG)
  
    S_hat       <- pec::predictSurvProb(S_fit, newdata = X_holdout, times = approx_times)
    G_hat       <- pec::predictSurvProb(G_fit, newdata = X_holdout, times = approx_times)
    f_hat       <- S_hat[,which(approx_times %in% landmark_times),drop=FALSE]
    S_hat_train <- pec::predictSurvProb(S_fit, newdata = X, times = approx_times)
    G_hat_train <- pec::predictSurvProb(G_fit, newdata = X, times = approx_times)
    f_hat_train <- S_hat_train[,which(approx_times %in% landmark_times),drop=FALSE]

  }
  # 
  # else if (nuisance == "survivalSL"){
  #   
  # }
  
  return(list(S_hat = S_hat,
              G_hat = G_hat,
              f_hat = f_hat,
              f_hat_train = f_hat_train,
              S_hat_train = S_hat_train,
              G_hat_train = G_hat_train))
}




generate_reduced_predictions <- function(f_hat,
                                         X_reduced,
                                         X_reduced_holdout){

  SL.library <- c("SL.mean", "SL.gam","SL.ranger")
  long_dat <- data.frame(f_hat = f_hat, X_reduced)
  long_new_dat <- data.frame(X_reduced_holdout)
  reduced_fit <- SuperLearner::SuperLearner(Y = long_dat$f_hat,
                                            X = long_dat[,2:ncol(long_dat),drop=FALSE],
                                            family = stats::gaussian(),
                                            SL.library = SL.library,
                                            method = "method.NNLS",
                                            verbose = FALSE)
  fs_hat <- matrix(predict(reduced_fit, newdata = long_new_dat)$pred,
                   nrow = nrow(X_reduced_holdout),
                   ncol = 1)
  
  return(list(fs_hat = fs_hat))
}

CV_generate_full_predictions_landmark <- function(time,
                                                  event,
                                                  X,
                                                  landmark_times,
                                                  approx_times,
                                                  nuisance,
                                                  cf_folds) {
  
  V <- length(unique(cf_folds))
  
  res <- purrr::map(seq_len(V), function(j) {
    
    train_id <- cf_folds != j
    test_id  <- cf_folds == j
    
    full_preds <- generate_full_predictions(
      time = time[train_id],
      event = event[train_id],
      X = X[train_id, , drop = FALSE],
      X_holdout = X[test_id, , drop = FALSE],
      
      landmark_times = landmark_times,
      approx_times = approx_times,
      nuisance = nuisance
    )
    
    list(
      CV_full_preds_train = full_preds$f_hat_train,
      CV_full_preds = full_preds$f_hat,
      CV_S_preds = full_preds$S_hat,
      CV_S_preds_train = full_preds$S_hat_train,
      CV_G_preds = full_preds$G_hat,
      CV_G_preds_train = full_preds$G_hat_train
    )
  })
  
  purrr::transpose(res)
}

CV_generate_reduced_predictions_landmark <- function(time,
                                                     event,
                                                     X,
                                                     landmark_times,
                                                     cf_folds,
                                                     indx,
                                                     full_preds_train) {
  
  V <- length(unique(cf_folds))
  
  purrr::map(seq_len(V), function(j) {
    
    train_id <- cf_folds != j
    test_id  <- cf_folds == j
    
    X_reduced_train <- X[train_id, -indx, drop = FALSE]
    X_reduced_holdout <- X[test_id, -indx, drop = FALSE]
    
    purrr::map_dfc(seq_along(landmark_times), function(k) {
      
      reduced_preds <- generate_reduced_predictions(
        f_hat = full_preds_train[[j]][, k],
        X_reduced = X_reduced_train,
        X_reduced_holdout = X_reduced_holdout
      )
      
      tibble::tibble(!!paste0("t", landmark_times[k]) := reduced_preds$fs_hat)
      
    }) %>% as.matrix()
  })
}
