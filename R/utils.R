

# ----------------------------- Documentation ------------------------------#

#' This file contains utils functions to:
#' 1. Make survival curves
#' 2. Compute Fisher pairwise test

# ----------------------------- Dependencies -------------------------------#

library(tidyverse)
library(survival)      
library(survminer) 
library(flextable)

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


#' Plot Kaplan–Meier survival curves (null or adjusted model)
#'
#' This function produce survival curves, either:
#'  - null model:   Surv(time, status) ~ 1       (no groups, no legend)
#'  - adjusted model: Surv(time, status) ~ group (groups as strata, legend shown)
#'
#' @param formula A survival formula, e.g. `Surv(time, status) ~ 1` or
#'   `Surv(time, status) ~ arm`.
#' @param data A data.frame containing the variables in `formula`.
#' @param title Plot title.
#' @param legend Position of the legend for adjusted models
#'   (e.g. `"top"`, `"bottom"`, `"right"`, `"left"`, `"none"`).
#' @param legend.title Legend title. 
#' @param xlab Label for x-axis (time).
#' @param ylab Label for y-axis (survival probability).
#' @param xlim Numeric length-2 vector for x-axis limits.
#' @param break.x.by Numeric: distance between x-axis breaks.
#' @param conf.int Logical: show confidence intervals.
#' @param risk.table Character or logical: type of risk table
#'   (e.g. `"nrisk_cumcensor"`, `TRUE`, or `"none"`).
#' @param risk.table.y.text Logical: show group labels in risk table y-axis.
#' @param risk.table.fontsize Numeric: font size for risk table.
#' @param tables.height Numeric: relative height of risk table vs plot.
#' @param surv.scale Character: `"percent"` or `"default"`.
#' @param axes.offset Logical: add offset to axes (ggsurvplot argument).
#' @param gg_theme A ggplot2 theme object for the main plot.
#' @param tables.theme A ggplot2 theme object for the risk table.
#' @param pval Logical: add log-rank p-value (for adjusted models).
#' @param ... Additional arguments passed to `ggsurvplot()`.
#'
#' @return A `ggsurvplot` object (list with ggplot & table components).

plot_surv <- function(formula,
                      data, 
                      title = "",
                      legend = "top", 
                      legend.title = "",
                      xlab = "Time in months",
                      ylab = "Survival probability",
                      xlim = c(0,40),  
                      break.x.by = 6,
                      conf.int = F, 
                      risk.table = "nrisk_cumcensor",
                      risk.table.y.text = F,
                      risk.table.fontsize = 4,
                      tables.height = 0.13,
                      surv.scale = "percent",
                      axes.offset = T,
                      gg_theme = theme_bw()+
                        theme(axis.title = element_text(size = 9),
                              axis.text  = element_text(size = 9)),
                      tables.theme = theme_cleantable()+
                        theme(plot.title = element_text(size = 10)),
                      pval = F,
                      ...){
  
  fit <- surv_fit(formula, data, match.fd = FALSE)
  plot <- vector(mode = "list", length(fit))
  
  for (i in seq_along(fit)){
    if (grepl("null_model", names(fit)[i])){ 
      
      plot[[i]] <- ggsurvplot(
        fit = fit[[i]], 
        legend = "none", 
        title = title, 
        xlab = xlab,
        ylab = ylab,
        xlim = xlim,  
        break.x.by = break.x.by,
        conf.int = conf.int, 
        risk.table = risk.table,
        risk.table.y.text = risk.table.y.text,
        risk.table.fontsize = risk.table.fontsize, 
        tables.height = tables.height,
        gg_theme = gg_theme, 
        tables.theme = tables.theme,
        surv.scale = surv.scale,
        axes.offset = axes.offset,
        ...
      ) } else {
        
        names(fit[[i]]$strata) <- sub(".*?=", "", names(fit[[i]]$strata)) 
        
        plot[[i]] <- 
          ggsurvplot(
            fit = fit[[i]], 
            xlab = xlab,
            ylab = ylab,
            title = title,
            legend = legend,
            legend.title = legend.title[i],
            xlim = xlim,  
            break.x.by = break.x.by,
            conf.int = conf.int, 
            risk.table = risk.table,
            risk.table.y.text = risk.table.y.text,
            risk.table.fontsize = risk.table.fontsize, 
            gg_theme = gg_theme, 
            tables.theme = tables.theme,
            surv.scale = surv.scale,
            axes.offset = axes.offset,
            pval = pval,
            ...
          )
      }
  }
  names(plot) <- gsub("::", "_", names(fit))
  plot
}

# Pour intégrer les formules aux flextables 
mixed_paragraph <- function(x) {
  
  # Repérage des expressions $...$
  positions <- stringr::str_locate_all(x, "\\$[^$]+\\$")[[1]]
  
  # Aucun code LaTeX dans la cellule
  if (nrow(positions) == 0L) {
    return(
      flextable::as_paragraph(
        flextable::as_chunk(x)
      )
    )
  }
  
  chunks <- list()
  start_text <- 1L
  
  for (k in seq_len(nrow(positions))) {
    
    eq_start <- positions[k, "start"]
    eq_end   <- positions[k, "end"]
    
    # Texte situé avant l’équation
    if (eq_start > start_text) {
      chunks[[length(chunks) + 1L]] <-
        flextable::as_chunk(
          substr(x, start_text, eq_start - 1L)
        )
    }
    
    # Équation sans les délimiteurs $
    equation <- substr(x, eq_start + 1L, eq_end - 1L)
    
    chunks[[length(chunks) + 1L]] <-
      flextable::as_equation(
        equation,
        width = 0.8,
        height = 0.22
      )
    
    start_text <- eq_end + 1L
  }
  
  # Texte situé après la dernière équation
  if (start_text <= nchar(x)) {
    chunks[[length(chunks) + 1L]] <-
      flextable::as_chunk(
        substr(x, start_text, nchar(x))
      )
  }
  
  do.call(flextable::as_paragraph, chunks)
}


tbl_vim <- function(output){
  
  full_model <- output %>%
    distinct(vim, large_predictiveness) %>%
    transmute(
      model = "Full model",
      vim,
      performance = large_predictiveness,
      importance = NA_real_,
      model_order = 0L
    )
  
  reduced_models <- output %>%
    transmute(
      model = paste0(
        "Without $X_{",
        stringr::str_remove(variable, "^X"),
        "}$"
      ),
      vim,
      performance = small_predictiveness,
      importance = est,
      model_order = readr::parse_number(variable)
    )
  
  tab <- bind_rows(full_model, reduced_models) %>%
    mutate(
      vim = recode(
        vim,
        "BS(t)"  = "BS",
        "AUC(t)" = "AUC"
      )
    ) %>%
    pivot_wider(
      id_cols = c(model, model_order),
      names_from = vim,
      values_from = c(performance, importance),
      names_glue = "{.value}_{vim}"
    ) %>%
    # arrange(model_order) %>%
    dplyr::select(
      model,
      performance_BS,
      importance_BS,
      performance_AUC,
      importance_AUC
    )%>%
    dplyr::mutate(
      across(2:5, ~ifelse(!is.na(.x), 
                          format(round(.x, 3), nsmall = 3), 
                          "—"))
    )%>%
    purrr::set_names(nm = c("Model",
                            "$\\widehat V_n^{\\mathrm{BS}}(\\tau)$",
                            "$\\widehat \\psi_{n,j}^{\\mathrm{BS}}(\\tau)$",
                            "$\\widehat V_n^{\\mathrm{AUC}}(\\tau)$",
                            "$\\widehat \\psi_{n,j}^{\\mathrm{AUC}}(\\tau)$"))
  
  ft_tab <- flextable::flextable(tab)
  
  for (i in seq_len(nrow(tab))) {
    ft_tab <- ft_tab %>%
      flextable::compose(
        i = i,
        j = "Model",
        value = mixed_paragraph(tab$Model[i]),
        part = "body"
      )
  }
  
  for (j in 2:ncol(tab)){
    ft_tab <- ft_tab %>%
      flextable::compose(
        j = j, 
        value = mixed_paragraph(names(tab)[j]),
        part = "header"
      )
  }
  
  notes <- "$\\widehat V_n^{\\mathrm{BS}}(\\tau)$: denotes the population-level predictive performance at horizon $\\tau$, defined as the negative Brier score; $\\widehat \\psi_{n,j}^{\\mathrm{BS}}(\\tau)$: denotes the loss in Brier-score-based predictive performance resulting from the exclusion of $X_j$; $\\widehat V_n^{\\mathrm{AUC}}(\\tau)$: denotes the population-level cumulative/dynamic AUC at horizon $\\tau$; $\\widehat \\psi_{n,j}^{\\mathrm{AUC}}(\\tau)$: denotes the loss in cumulative/dynamic AUC resulting from the exclusion of $X_j$."
  
  ft_tab %>%
    flextable::add_footer_lines(notes)%>%
    flextable::compose(
      value = mixed_paragraph(notes),
      part = "footer"
    )
}
