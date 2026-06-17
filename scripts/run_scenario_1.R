

#---------------------------#
#         SCENARIO 1        #
#---------------------------#

run_s1 <- function(n, nuisance, seed){
  
  start <- Sys.time()
  
  get_params_scenario(scenario = "1",
                      n_mc     = n, 
                      seed     = seed)
  tau <- 14.97 
  data <- generate_data(n = n, 
                        scenario = "1",
                        seed     = seed)
  vims <- 
}