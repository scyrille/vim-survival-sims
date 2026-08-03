
simulation_metrics <- function(scenario){
  
  truth_param <- readRDS(here::here("outputs","results","mc_truth_param.RDS"))
  
  output_files <- list.files(
    path = here::here("outputs","results", paste0("scenario", scenario)), 
    full.names = T
  )
 
  output <- lapply(output_files, readRDS)%>%
    map(~ .x %>%
        bind_rows())%>%
    bind_rows()%>%
    left_join(truth_param %>% dplyr::filter(scenario == scenario), 
              by = c("scenario","vim","variable","tau"))%>%
    dplyr::mutate(err = est - true_vim,
                  cov = ifelse(cil <= true_vim & ciu >= true_vim, 1, 0),
                  reject = ifelse(p < 0.05, 1, 0),
                  width = (est - cil)*2)%>%
    dplyr::group_by(tau, n, nuisance, vim, variable)%>% 
    dplyr::summarize(
      runtime = mean(runtime),
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
      .groups = "drop"
    )
  return(output)
}

