
packages <- c("here","tidyverse","MASS","survML","SuperLearner",
              "ranger","xgboost","timereg","survivalSL","parallel")

invisible(lapply(packages, library, character.only = TRUE))

lapply(list.files(here::here("R"), full.names = T), source)

# Calibration parameters 
calibration_param <- readRDS(here::here("outputs","results",
                                        "mc_calibration_param.rds"))
c_max <- calibration_param$params_2$c_max
tau <- calibration_param$params_2$tau

# Number of Monte Carlo replications 
# R <- 500
R <- 5

# Sample size 
# n <- c(500,1000)
n <- 1000

nuisances <- c("cox.aalen", "stackG", "survivalSL")

param_grid <- expand.grid(
  n = n,
  stringsAsFactors = FALSE
)

RNGkind("L'Ecuyer-CMRG")
set.seed(2026)
param_grid$seed <- sample.int(1e9, nrow(param_grid))

workers <- min(R, 4L)  

out_dir <- here::here("outputs", "results", "scenario2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

run_replication <- function(rng, n, c_max, tau, nuisances) {
  assign(".Random.seed", rng, envir = .GlobalEnv)
  
  do_one(
    n = n,
    scenario = "2",
    c_max = c_max,
    tau = tau,
    nuisance = nuisances
  )
}

fun_names <- as.character(
  utils::lsf.str(envir = .GlobalEnv, all.names = TRUE)
)

export_names <- unique(c(
  fun_names,
  "c_max", "tau", "nuisances"
))

cl <- parallel::makePSOCKcluster(workers)

tryCatch({
  
  # Charger les packages sur chaque processus
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
    
    set.seed(seed)
    rng_seeds <- vector("list", R)
    rng_seeds[[1L]] <- .Random.seed
    
    for (r in seq_len(R)[-1L]) {
      rng_seeds[[r]] <- parallel::nextRNGStream(
        rng_seeds[[r - 1L]]
      )
    }
    
    res <- parallel::parLapplyLB(
      cl,
      X = rng_seeds,
      fun = run_replication,
      n = n,
      c_max = c_max,
      tau = tau,
      nuisances = nuisances
    )
    
    saveRDS(
      res,
      file.path(out_dir, paste0("sims_n", n, ".rds"))
    )
  })
  
}, finally = {
  parallel::stopCluster(cl)
})
