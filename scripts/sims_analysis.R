
# ---------------------------------------------------- #
#                S1_base simulations results           #
# ---------------------------------------------------- #

# Packages and project functions -----------------------------------------

packages <- c("here", "purrr", "ggplot2")

invisible(lapply(packages, library, character.only = TRUE))

source_files <- list.files(
  here::here("R"),
  pattern = "\\.[Rr]$",
  full.names = TRUE
)

invisible(lapply(source_files, source))

res_dir <- here::here("outputs", "results")
fig_dir <- here::here("outputs", "figures")

# S1_base -----------------------------------------------------------------

# True VIM values
true_param_S1_base <- readRDS(here::here(res_dir, "S1_base/mc_truth_param_S1_base.rds"))

# n = 500
sims_n500_S1_base <- list.files(
  path = here::here(res_dir, "S1_base/n500_seed837949287"), 
  full.names = T)%>%
  purrr::map(~readRDS(.x))%>%
  bind_rows()

perf_sims_n500_S1_base <- simulation_metrics(
  truth_param = true_param_S1_base,
  output      = sims_n500_S1_base)

plot_perf(perf_sims_n500_S1_base)%>%
  ggsave(filename = file.path(fig_dir, "perf_sims_n500_S1_base.png"),
         width = 13, height = 8, dpi = 600, units = "in", bg = "white")

# n = 1000 
sims_n1000_S1_base <- list.files(
  path = here::here(res_dir, "S1_base/n1000_seed837949287"),
  full.names = T)%>%
  purrr::map(~readRDS(.x))%>%
  bind_rows()

perf_sims_n1000_S1_base <- simulation_metrics(
  truth_param = true_param_S1_base,
  output      = sims_n1000_S1_base)

plot_perf(perf_sims_n1000_S1_base)%>%
  ggsave(filename = file.path(fig_dir, "perf_sims_n1000_S1_base.png"),
         width = 13, height = 8, dpi = 600, units = "in", bg = "white")
