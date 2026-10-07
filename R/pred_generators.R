
library(tidyr)
library(purrr)
library(survML)
library(survivalSL)

generate_full_predictions <- function(time, 
                                      event, 
                                      X, 
                                      X_holdout, 
                                      tau, 
                                      approx_times,
                                      nuisance) {
  
  
  event <- as.integer(event)
  X <- as.data.frame(X)
  X_holdout <- as.data.frame(X_holdout)
  X_holdout <- X_holdout[, names(X), drop = FALSE]
  newX <- rbind(X_holdout, X)
  dat <- data.frame(time = time, event = event, X,
                    check.names = FALSE)
  datS <- dat
  datG <- dat %>%
    dplyr::mutate(event = 1L - event)
  
  if (nuisance == "stackG") {
    
    SL.library <- c(
      "SL.mean", 
      "SL.glm", 
      "SL.gam",
      "SL.ranger",
      "SL.xgboost"
    )
    
    fit <- survML::stackG(
      time = time, 
      event = event, 
      X = X, 
      newX = newX,
      newtimes = approx_times,
      time_grid_approx = approx_times,
      bin_size = 0.05,
      time_basis = "continuous",
      surv_form = "PI",
      SL_control = list(SL.library = SL.library, 
                        V = 5)
    )
    
    S_all <- fit$S_T_preds
    G_all <- fit$S_C_preds
    
  } else {
    
    if (nuisance == "aalen") {
      
      S_fit <- timereg::aalen(
        survival::Surv(time, event) ~ ., data = datS
      )
      # G_fit <- timereg::aalen(
      #   survival::Surv(time, event) ~ ., data = datG
      # )
      G_fit <- survival::survfit(
        survival::Surv(time, event) ~ 1, data = datG
      )
      
    } else if (nuisance == "cox.aalen") {
      
      x_vars <- names(X)
      prop_vars <- x_vars[grepl("Z",x_vars)]
      add_vars  <- x_vars[grepl("X",x_vars)]
      
      rhs <- c(
        if (length(prop_vars) > 0) paste0("prop(", prop_vars, ")"),
        add_vars
      )
      
      form <- stats::as.formula(
        paste("survival::Surv(time, event) ~",
              paste(rhs, collapse = " + ")),
        env = environment()
      )
      
      # Event model: S(t | X)
      S_fit <- timereg::cox.aalen(form, data = datS)
      
      # Censoring model: G(t | X)
      # G_fit <- timereg::cox.aalen(form, data = datG)
      G_fit <- survival::survfit(
        survival::Surv(time, event) ~ 1, data = datG
      )
      
    } else if (nuisance == "survivalSL") {  
      
      methods <- c(
        "LIB_COXall", 
        "LIB_COXridge", 
        "LIB_PHexponential",
        "LIB_AFTgamma", 
        "LIB_RSF"
      )
      
      # Event model: S(t | X)
      S_fit <- survivalSL::survivalSL(
        formula       = survival::Surv(time, event)~.,
        data          = datS,
        metric        = "auc",
        methods       = methods,
        cv            = 5,
        show_progress = FALSE
      )
      
      # Censoring model: G(t | X)
      G_fit <- survivalSL::survivalSL(
        formula       = survival::Surv(time, event)~.,
        data          = datG,
        metric        = "auc",
        methods       = methods,
        cv            = 5,
        show_progress = FALSE
      )
      
    }
    
    if (nuisance == "survivalSL") {
      
      predict_sl <- function(fit) {
        
        times_pos <- sort(unique(approx_times[approx_times > 0]))
        out <- matrix(1, nrow = nrow(newX), ncol = length(approx_times))
        
        if (length(times_pos) > 0L) {
          pred <- stats::predict(
            fit,
            newdata = as.data.frame(newX),
            newtimes = times_pos
          )
          mat <- pred$predictions[["sl"]]
          cols <- match(approx_times[approx_times > 0], pred$times)
          
          out[, approx_times > 0] <- mat[, cols, drop = FALSE]
        }
        
        out
      }
      
      S_all <- predict_sl(S_fit)
      G_all <- predict_sl(G_fit)
      
    } else {
      
      S_all <- pec::predictSurvProb(
        S_fit, newdata = newX, times = approx_times
      )
      # G_all <- pec::predictSurvProb(
      #   G_fit, newdata = newX, times = approx_times
      # )
      G_summary <- summary(
        G_fit,
        times = sort(unique(approx_times)),
        extend = TRUE
      )
      
      # Restaurer l'ordre initial des temps, y compris les doublons
      G_t <- G_summary$surv[match(approx_times, G_summary$time)]
      G_all <- matrix(
        G_t,
        nrow = nrow(newX),
        ncol = length(approx_times),
        byrow = TRUE
      )
    }
  }
  
  expected_dim <- c(nrow(newX), length(approx_times))
  stopifnot(
    is.matrix(S_all), identical(dim(S_all), expected_dim),
    is.matrix(G_all), identical(dim(G_all), expected_dim)
  )
  
  i_holdout <- seq_len(nrow(X_holdout))
  i_train <- nrow(X_holdout) + seq_len(nrow(X))
  j_tau <- match(tau, approx_times)
  
  S_hat <- S_all[i_holdout, , drop = FALSE]
  G_hat <- G_all[i_holdout, , drop = FALSE]
  S_hat_train <- S_all[i_train, , drop = FALSE]
  G_hat_train <- G_all[i_train, , drop = FALSE]
  
  list(
    S_hat = S_hat,
    G_hat = G_hat,
    f_hat = S_hat[, j_tau, drop = FALSE],
    f_hat_train = S_hat_train[, j_tau, drop = FALSE],
    S_hat_train = S_hat_train,
    G_hat_train = G_hat_train
  )
}


generate_reduced_predictions <- function(f_hat,
                                         X_reduced,
                                         X_reduced_holdout){

  SL.library <- c("SL.mean", 
                  "SL.glm", 
                  "SL.gam",
                  "SL.ranger",
                  "SL.xgboost")
  
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

CV_generate_full_predictions <- function(time,
                                         event,
                                         X,
                                         tau,
                                         approx_times,
                                         nuisance,
                                         cf_folds) {
  
  X <- as.data.frame(X)
  event <- as.integer(event)
  
  V <- length(unique(cf_folds))
  
  res <- purrr::map(seq_len(V), function(j) {
    
    train_id <- cf_folds != j
    test_id  <- cf_folds == j
    
    full_preds <- generate_full_predictions(
      time = time[train_id],
      event = event[train_id],
      X = X[train_id, , drop = FALSE],
      X_holdout = X[test_id, , drop = FALSE],
      tau = tau,
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

CV_generate_reduced_predictions <- function(time,
                                            event,
                                            X,
                                            tau,
                                            cf_folds,
                                            indx,
                                            full_preds_train) {
  
  V <- length(unique(cf_folds))
  
  purrr::map(seq_len(V), function(j) {
    
    train_id <- cf_folds != j
    test_id  <- cf_folds == j
    
    X_reduced_train <- X[train_id, -indx, drop = FALSE]
    X_reduced_holdout <- X[test_id, -indx, drop = FALSE]
    
    purrr::map_dfc(seq_along(tau), function(k) {
      
      reduced_preds <- generate_reduced_predictions(
        f_hat = full_preds_train[[j]][, k],
        X_reduced = X_reduced_train,
        X_reduced_holdout = X_reduced_holdout
      )
      
      reduced_preds
      tibble::tibble(!!paste0("t", tau[k]) := reduced_preds$fs_hat)
      
    }) %>% as.matrix()
  })
}
