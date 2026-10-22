# ---------------------------------------------------- #
#                 S1_prev simulations                  #
#        Additive model — varying prevalence           #
# ---------------------------------------------------- #

# Run this script in the project containing the functions in R/.
# Each replication is saved in a separate RDS file.
# Restarting the simulation section computes only replications with missing files.
# Reuse existing files only when the code and parameters are unchanged.
# Otherwise, choose a new out_dir.

# Packages and project functions -----------------------------------------

packages <- c(
  "here", "tidyverse", "MASS", "cvAUC", "vimp", "survML",
  "SuperLearner", "ranger", "xgboost", "timereg",
  "survivalSL", "parallel"
)

invisible(lapply(packages, library, character.only = TRUE))

source_files <- list.files(
  here::here("R"),
  pattern = "\\.[Rr]$",
  full.names = TRUE
)

invisible(lapply(source_files, source))

out_dir <- here::here("outputs", "results", "S1_prev")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir, mustWork = TRUE)

# Set the RNG explicitly to keep the Monte Carlo calculations reproducible
# when rerunning the entire script in the same R session.
RNGkind(
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)
set.seed(2026)

n_mc <- 2e7

# Calibration parameters -------------------------------------------------

param_S1_prev <- calibrate_parameters(
  n = n_mc,
  scenario = "S1_prev",
  max_upper = 1e6,
  probs = 0.75
)

c_max <- param_S1_prev$c_max
tau <- param_S1_prev$tau

saveRDS(
  param_S1_prev,
  file.path(out_dir, "mc_calibration_param_S1_prev.rds")
)

# True VIM values ---------------------------------------------------------

# Generate data
data_S1_prev <- generate_data(
  n = n_mc,
  scenario = "S1_prev",
  c_max = c_max
)

# True VIM values
true_param_S1_prev <- purrr::map_dfr(tau, function(t) {
  
  y <- as.integer(data_S1_prev$t > t)
  
  xvars <- names(
    dplyr::select(data_S1_prev, dplyr::starts_with("X"))
  )
  
  # Full model
  f_full <- true_S0_aalen(
    t = t,
    data = data_S1_prev,
    lambda0 = lambda0,
    alpha = alpha
  )
  
  V_full_brier <- -vimp::measure_mse(f_full, y)$point_est
  V_full_auc <- cvAUC::AUC(f_full, y)
  
  # Reduced models
  reduced_res <- purrr::map_dfr(seq_along(xvars), function(j) {
    
    xj <- xvars[j]
    
    f_red <- true_S0_aalen(
      t = t,
      data = dplyr::select(
        data_S1_prev,
        -dplyr::all_of(xj)
      ),
      lambda0 = lambda0,
      alpha = alpha[-j]
    )
    
    tibble::tibble(
      variable = xj,
      V_brier = -vimp::measure_mse(f_red, y)$point_est,
      V_auc = cvAUC::AUC(f_red, y)
    )
  })
  
  dplyr::bind_rows(
    reduced_res %>%
      dplyr::transmute(
        vim = "BS(t)",
        tau = t,
        variable,
        true_V_full = V_full_brier,
        true_V_reduced = V_brier,
        true_vim = V_full_brier - V_brier
      ),
    reduced_res %>%
      dplyr::transmute(
        vim = "AUC(t)",
        tau = t,
        variable,
        true_V_full = V_full_auc,
        true_V_reduced = V_auc,
        true_vim = pmax(V_full_auc - V_auc, 0)
      )
  )
}) %>%
  dplyr::mutate(
    n_mc = n_mc,
    scenario = "S1_prev"
  ) %>%
  dplyr::relocate(scenario) %>%
  dplyr::arrange(variable, tau)

saveRDS(
  true_param_S1_prev,
  file.path(out_dir, "mc_truth_param_S1_prev.rds")
)

# Release the Monte Carlo data set before starting the parallel workers.
rm(data_S1_prev)
invisible(gc())

# Run simulations ---------------------------------------------------------
#
# To resume only the simulations, first load the packages and source the
# project functions, then reload the saved calibration parameters:
#
# param_S1_base <- readRDS(
#   file.path(out_dir, "mc_calibration_param_S1_base.rds")
# )
# c_max <- param_S1_base$c_max
# tau <- param_S1_base$tau
#
# Then run this section. Rerunning the entire script also recomputes the
# calibration parameters and the Monte Carlo truth values.

# Number of Monte Carlo replications
R <- 500L

# Sample sizes
# n <- c(500,1000)
n <- 1000

nuisances <- c("aalen", "stackG","survivalSL")

# Add any other global constants used by the project functions, if needed.
# Do not export the large Monte Carlo data set.
extra_export_names <- character()

param_grid <- expand.grid(
  n = n,
  stringsAsFactors = FALSE
)

RNGkind("L'Ecuyer-CMRG")
set.seed(2026)
param_grid$seed <- sample.int(1e9, nrow(param_grid))

workers <- 8L

# Compute and save one replication 

run_replication <- function(job, n, c_max, tau, nuisances) {
  
  # The RNG stream is tied to the replication, regardless of which worker
  # receives the task from parLapplyLB.
  assign(".Random.seed", job$rng, envir = .GlobalEnv)
  
  result <- do_one(
    n = n,
    scenario = "S1_prev",
    c_max = c_max,
    tau = tau,
    nuisance = nuisances
  )
  
  # The final file is created only after writing has completed.
  # Files ending in .tmp are not treated as completed results.
  tmp <- paste0(job$file, ".tmp")
  on.exit(unlink(tmp), add = TRUE)
  
  saveRDS(result, tmp)
  
  if (!file.rename(tmp, job$file)) {
    stop("Failed to save: ", job$file)
  }
  
  # Return only the path to avoid sending large objects to the main process.
  job$file
}

# Export functions and global parameters

fun_names <- as.character(
  utils::lsf.str(envir = .GlobalEnv, all.names = TRUE)
)

export_names <- unique(c(
  fun_names,
  "c_max", "tau", "nuisances",
  extra_export_names
))

cl <- parallel::makePSOCKcluster(workers)

tryCatch({
  
  # Load packages on each worker
  parallel::clusterCall(
    cl,
    function(pkgs, lib_paths) {
      .libPaths(lib_paths)
      
      for (pkg in pkgs) {
        library(pkg, character.only = TRUE)
      }
      
      NULL
    },
    pkgs = packages,
    lib_paths = .libPaths()
  )
  
  parallel::clusterExport(
    cl,
    varlist = export_names,
    envir = .GlobalEnv
  )
  
  purrr::pwalk(param_grid, function(n, seed) {
    
    # Rebuild all RNG streams, including when resuming, to preserve the same
    # stream for each replication index.
    set.seed(seed)
    
    rng_seeds <- vector("list", R)
    rng_seeds[[1L]] <- .Random.seed
    
    for (r in seq_len(R)[-1L]) {
      rng_seeds[[r]] <- parallel::nextRNGStream(
        rng_seeds[[r - 1L]]
      )
    }
    
    # One directory per n / seed combination
    rep_dir <- file.path(
      out_dir,
      sprintf("n%s_seed%s", n, seed)
    )
    
    dir.create(
      rep_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
    
    rep_dir <- normalizePath(rep_dir, mustWork = TRUE)
    
    # Keep file names in the original replication order.
    files <- file.path(
      rep_dir,
      sprintf("rep_%04d.rds", seq_len(R))
    )
    
    # Schedule only missing replications
    pending <- which(!file.exists(files))
    
    message(sprintf(
      "n = %s: %d replications to compute out of %d",
      n, length(pending), R
    ))
    
    if (length(pending) > 0L) {
      
      jobs <- lapply(pending, function(r) {
        list(
          rng = rng_seeds[[r]],
          file = files[[r]]
        )
      })
      
      invisible(parallel::parLapplyLB(
        cl,
        X = jobs,
        fun = run_replication,
        n = n,
        c_max = c_max,
        tau = tau,
        nuisances = nuisances,
        chunk.size = 1L
      ))
    }
  })
  
}, finally = {
  
  # Stop the workers even if a replication or save operation fails.
  parallel::stopCluster(cl)
  
})
