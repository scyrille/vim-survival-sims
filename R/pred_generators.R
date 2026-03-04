#' Estimate conditional survival function nuisance parameters using the Aalen additive hazards model
#'
#' @param time \code{n x 1} numeric vector of observed
#' follow-up times. If there is censoring, these are the minimum of the
#' event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed.
#' @param X \code{n x p} data.frame of observed covariate values
#' @param X_holdout \code{m x p} data.frame of new observed covariate
#' values at which to obtain \code{m} predictions for the estimated algorithm.
#' Must have the same names and structure as \code{X}.
#' @param newtimes \code{k x 1} numeric vector of times at which to obtain \code{k}
#' predicted conditional survivals.
#'
#' @return A list containing elements \code{S_hat} (conditional event survival function, corresponding to \code{X_holdout} and \code{newtimes}),
#' \code{S_hat_train} (conditional event survival function, corresponding to \code{X} and \code{newtimes}),
#' \code{G_hat} (conditional censoring survival function, corresponding to \code{X_holdout} and \code{newtimes}),
#' and \code{G_hat_train} (conditional censoring survival function, corresponding to \code{X} and \code{newtimes})
#'
#' @export
conditional_surv_aalen <- function(time, event, X, X_holdout, newtimes){
  
  # ---- Event model: S(t|X) ----
  datS <- data.frame(time = time, event = event, X)
  S_fit <- timereg::aalen(survival::Surv(time, event) ~ ., data = datS)
  
  S_hat       <- pec::predictSurvProb(S_fit, newdata = X_holdout, times = newtimes)
  S_hat_train <- pec::predictSurvProb(S_fit, newdata = X,         times = newtimes)
  
  # ---- Censoring model: G(t|X) ----
  cens_event <- 1 - event
  datG <- data.frame(time = time, cens_event = cens_event, X)
  G_fit <- timereg::aalen(survival::Surv(time, cens_event) ~ ., data = datG)
  
  G_hat       <- pec::predictSurvProb(G_fit, newdata = X_holdout, times = newtimes)
  G_hat_train <- pec::predictSurvProb(G_fit, newdata = X,         times = newtimes)
  
  list(
    S_hat = S_hat,
    G_hat = G_hat,
    S_hat_train = S_hat_train,
    G_hat_train = G_hat_train
  )
}

#' Estimate conditional survival function nuisance parameters using survival stacking
#'
#' @param time \code{n x 1} numeric vector of observed
#' follow-up times. If there is censoring, these are the minimum of the
#' event and censoring times.
#' @param event \code{n x 1} numeric vector of status indicators of
#' whether an event was observed.
#' @param X \code{n x p} data.frame of observed covariate values
#' @param X_holdout \code{m x p} data.frame of new observed covariate
#' values at which to obtain \code{m} predictions for the estimated algorithm.
#' Must have the same names and structure as \code{X}.
#' @param approx_times Numeric vector of times at which to
#' approximate product integral or cumulative hazard interval. See [stackG] documentation.
#' @param SL.library Super Learner library
#' @param V Number of cross-validation folds, to be passed to \code{SuperLearner}
#' @param newtimes \code{k x 1} numeric vector of times at which to obtain \code{k}
#' predicted conditional survivals.
#' @param bin_size Size of time bin on which to discretize for estimation
#' of cumulative probability functions. Can be a number between 0 and 1,
#' indicating the size of quantile grid (e.g. \code{0.1} estimates
#' the cumulative probability functions on a grid based on deciles of
#' observed \code{time}s). If \code{NULL}, creates a grid of
#' all observed \code{time}s. See [stackG] documentation.
#'
#' @return A list containing elements \code{S_hat} (conditional event survival function, corresponding to \code{X_holdout} and \code{newtimes}),
#' \code{S_hat_train} (conditional event survival function, corresponding to \code{X} and \code{newtimes}),
#' \code{G_hat} (conditional censoring survival function, corresponding to \code{X_holdout} and \code{newtimes}),
#' and \code{G_hat_train} (conditional censoring survival function, corresponding to \code{X} and \code{newtimes})
#'
#' @seealso [stackG]
#'
#' @export
generate_nuisance_predictions_stackG <- function(time,
                                                 event,
                                                 X,
                                                 X_holdout,
                                                 newtimes,
                                                 SL.library,
                                                 V,
                                                 bin_size,
                                                 approx_times){
  
  surv_out <- survML::stackG(time = time,
                             event = event,
                             X = X,
                             newX = rbind(X_holdout, X),
                             newtimes = newtimes,
                             time_grid_approx = approx_times,
                             bin_size = bin_size,
                             time_basis = "continuous",
                             surv_form = "PI",
                             SL_control = list(SL.library = SL.library,
                                               V = V))
  S_hat <- surv_out$S_T_preds[1:nrow(X_holdout),]
  G_hat <- surv_out$S_C_preds[1:nrow(X_holdout),]
  S_hat_train <- surv_out$S_T_preds[(nrow(X_holdout)+1):(nrow(X_holdout)+nrow(X)),]
  G_hat_train <- surv_out$S_C_preds[(nrow(X_holdout)+1):(nrow(X_holdout)+nrow(X)),]
  return(list(S_hat = S_hat,
              G_hat = G_hat,
              S_hat_train = S_hat_train,
              G_hat_train = G_hat_train))
}

generate_oracle_predictions_DR <- function(time,
                                           event,
                                           X,
                                           X_holdout,
                                           nuisance_preds,
                                           outcome,
                                           landmark_times,
                                           restriction_time,
                                           approx_times,
                                           SL.library,
                                           V,
                                           indx){
  
  S_hat <- nuisance_preds$S_hat_train
  G_hat <- nuisance_preds$G_hat_train
  
  if (sum(indx) != 0){ # only remove column if there is a column to remove
    X <- X[,-indx,drop=FALSE]
    X_holdout <- X_holdout[,-indx,drop=FALSE]
  }
  
  if (outcome == "survival_probability"){
    newtimes <- landmark_times
  } else if (outcome == "restricted_survival_time"){
    newtimes <- restriction_time
  }
  DR_predictions_combined <- DR_pseudo_outcome_regression(time = time,
                                                          event = event,
                                                          X = X,
                                                          newX = rbind(X_holdout, X),
                                                          S_hat = S_hat,
                                                          G_hat = G_hat,
                                                          newtimes = newtimes,
                                                          outcome = outcome,
                                                          approx_times = approx_times,
                                                          SL.library = SL.library,
                                                          V = V)
  
  DR_predictions <- DR_predictions_combined[1:nrow(X_holdout),]
  DR_predictions_train <- DR_predictions_combined[(nrow(X_holdout) + 1):nrow(DR_predictions_combined),]
  
  if (outcome == "survival_probability"){
    f0_hat <- 1 - DR_predictions
    f0_hat_train <- 1 - DR_predictions_train
    if (length(landmark_times) == 1){
      f0_hat <- matrix(f0_hat, ncol = 1)
      f0_hat_train <- matrix(f0_hat_train, ncol = 1)
    }
  } else if (outcome == "restricted_survival_time"){
    f0_hat <- DR_predictions
    f0_hat_train <- DR_predictions_train
  }
  
  return(list(f0_hat = f0_hat,
              f0_hat_train = f0_hat_train))
  
}