
LIB_COXall_fixed <- function(formula, data) {
  
  data <- as.data.frame(data)
  
  variables_formula <- all.vars(formula)
  time_var <- variables_formula[[1L]]
  
  # model = TRUE empêche survfit() de devoir reconstruire
  # le model.frame à partir de l'appel data = data
  fit <- survival::coxph(
    formula = formula,
    data    = data,
    model   = TRUE
  )
  
  fit_surv <- survival::survfit(
    fit,
    newdata = data,
    se.fit  = FALSE
  )
  
  prediction_times <- sort(unique(data[[time_var]]))
  
  predictions <- t(
    summary(
      fit_surv,
      times = prediction_times
    )$surv
  )
  
  predictions <- cbind(
    rep(1, nrow(predictions)),
    predictions
  )
  
  result <- list(
    model       = fit,
    library     = "LIB_COXall",
    formula     = formula,
    data        = data,
    times       = c(0, prediction_times),
    predictions = predictions
  )
  
  class(result) <- "libsl"
  result
}

# Remplacement uniquement pour la session R courante
assignInNamespace(
  x     = "LIB_COXall",
  value = LIB_COXall_fixed,
  ns    = "survivalSL"
)
