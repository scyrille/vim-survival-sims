
library(dplyr)
library(ggplot2)
library(cowplot)

simulation_metrics <- function(output){
  
  truth_param <- readRDS(here::here("outputs","results","mc_truth_param.RDS"))
  
  scenario <- output %>%
    dplyr::select(scenario)%>%
    dplyr::distinct()
 
  perf <- output %>%
    # map(~ .x %>% bind_rows())%>%
    # bind_rows()%>%
    left_join(truth_param %>% dplyr::filter(scenario == scenario), 
              by = c("scenario","vim","variable","tau"))%>%
    dplyr::mutate(err = est - true_vim,
                  cov = ifelse(cil <= true_vim & ciu >= true_vim, 1, 0),
                  reject = ifelse(p < 0.05, 1, 0),
                  width = (est - cil)*2)%>%
    dplyr::group_by(tau, n, nuisance, vim, variable)%>% 
    dplyr::summarize(
      bias = mean(err),
      nreps = n(),
      variance = var(est),
      bias_mc_se = sqrt( mean((err - bias) ^ 2) / (nreps-1)),
      coverage = mean(cov),
      power = mean(reject),
      cov_mc_se = sqrt(coverage * (1 - coverage)/nreps),
      power_mc_se = sqrt(power * (1 - power)/nreps),
      ci_width = mean(width),
      width_mc_se = sqrt( mean((width - ci_width) ^ 2) / (nreps-1)),
      var_mc_se = sqrt(2/(nreps-1)),
      runtime = mean(runtime),
      .groups = "drop"
    )
  return(perf)
}


plot_perf <- function(sims_perf){
  
  point_size <- 1
  axis_text_size <- 12
  legend_text_size_small_plot <- 13
  
  sims_perf <- sims_perf %>%
    dplyr::mutate(
      nuisance = factor(nuisance, 
                        levels = c("aalen",
                                   "cox.aalen",
                                   "stackG",
                                   "survivalSL"),
                        labels = c("Aalen",
                                   "Cox-Aalen",
                                   "Global survival stacking",
                                   "Survival Super Learner")))
  theme_plot <- 
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "grey100"),
      panel.grid.minor.y = element_blank(),
      legend.position = "none",
      axis.text.y = element_text(size = axis_text_size),
      # axis.text.x = element_blank(),
      # axis.title.x = element_blank(),
      plot.margin = unit(c(0.1, 0.3, 0.1, 0), "cm")
    )
  
  
  # Bias
  bias_plot <- sims_perf %>% 
    ggplot(aes(x = variable, y = bias, color = nuisance)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "gray30") +
    geom_point(size = point_size) +
    geom_errorbar(aes(ymin=bias-1.96*bias_mc_se, ymax=bias + 1.96*bias_mc_se), width=.1) +
    facet_wrap(~ vim, labeller = label_parsed, strip.position = "top") +
    labs(x = "Variable", y = "Empirical bias", color = "Nuisance")+
    theme_plot
  
  # Coverage
  cover_plot <- sims_perf %>% 
    ggplot(aes(x = variable, y = coverage, color = nuisance)) +
    geom_hline(yintercept = 0.95, linetype = "solid", color = "gray30") +
    geom_point(size = point_size) +
    geom_errorbar(aes(ymin=ifelse(coverage-1.96*cov_mc_se <0, 0, coverage-1.96*cov_mc_se),
                      ymax=ifelse(coverage+1.96*cov_mc_se >1, 1, coverage+1.96*cov_mc_se)),
                  width=.1)+
    facet_wrap(~ vim, labeller = label_parsed, strip.position = "top") +
    labs(x = "Variable", y = "Empirical coverage", color = "Nuisance")+
    theme_plot
  
  # Variance 
  var_plot <- sims_perf %>% 
    ggplot(aes(x = variable, y = variance, color = nuisance)) +
    geom_point(size = point_size) +
    geom_errorbar(aes(ymin=variance - 1.96*var_mc_se,
                      ymax=variance + 1.96*var_mc_se),
                  width=.1) +
    facet_wrap(~ vim, labeller = label_parsed, strip.position = "top") +
    labs(x = "Variable", y = "Empirical variance", color = "Nuisance")+
    theme_plot
  
  # CI width
  width_plot <- sims_perf %>% 
    ggplot(aes(x = variable, y = ci_width, color = nuisance)) +
    geom_point(size = point_size) +
    geom_errorbar(aes(ymin = ci_width - 1.96*width_mc_se,
                      ymax = ci_width + 1.96*width_mc_se),
                  width=.1) +
    facet_wrap(~ vim, labeller = label_parsed, strip.position = "top") +
    labs(x = "Variable", y = "Empirical confidence interval width", color = "Nuisance")+
    theme_plot
  
  # Power 
  power_plot <- sims_perf %>% 
    ggplot(aes(x = variable, y = power, color = nuisance)) +
    geom_point(size = point_size) +
    geom_errorbar(aes(ymin=ifelse(power-1.96*power_mc_se <0, 0, power-1.96*power_mc_se),
                      ymax=ifelse(power+1.96*power_mc_se >1, 1, power+1.96*power_mc_se)),
                  width=.1) +
    geom_hline(yintercept = 0.05, linetype = "solid", color = "gray30") +
    facet_wrap(~ vim, labeller = label_parsed, strip.position = "top") +
    labs(x = "Variable", y = "Empirical power", color = "Nuisance")+
    theme_plot
  
  legend_j <- get_legend(
    bias_plot +
      guides(color = guide_legend(nrow = 1, ncol = 3)) +
      theme(legend.direction = "horizontal",
            legend.position = "bottom",
            legend.title = element_text(size = legend_text_size_small_plot),
            legend.text = element_text(size = legend_text_size_small_plot))
  )
 
  panel_plots <- cowplot::plot_grid(
    bias_plot, 
    var_plot, 
    cover_plot, 
    width_plot,
    power_plot,
    nrow = 5, 
    ncol = 1
  )
  
  cowplot::plot_grid(
    panel_plots, 
    legend_j, 
    ncol = 1, 
    nrow = 2,
    rel_heights = c(1, .075)
  )
}
