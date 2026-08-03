
library(here)
library(tidyverse)

# Load true parameters ----------------------------------------------------

truth_param <- readRDS(here::here("outputs","results","mc_truth_param.rds"))

# Load results ------------------------------------------------------------

output_files_1 <- list.files(path = here::here("outputs",
                                                  "results",
                                                  "scenario1"), 
                             full.names = T)
res_1 <- lapply(output_files_1, readRDS)%>%
  map(~ .x %>%
        bind_rows()%>%
        left_join(truth_param, by = c("scenario","vim"))%>%
        dplyr::mutate(err = (est - true_vim),
                      cov = ifelse(cil <= true_vim & ciu >= true_vim, 1, 0),
                      reject = ifelse(p < 0.05, 1, 0),
                      width = (est - cil)*2)%>%
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
          rank_cor = cor(est, true_vim, method = "pearson",
                         use = "complete.obs"), 
          .groups = "drop"))

