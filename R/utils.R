

# ----------------------------- Documentation ------------------------------#

#' This file contains utils functions to:
#' 1. Make survival curves
#' 2. Compute Fisher pairwise test

# ----------------------------- Dependencies -------------------------------#

library(tidyverse)
library(survival)      
library(survminer) 

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
