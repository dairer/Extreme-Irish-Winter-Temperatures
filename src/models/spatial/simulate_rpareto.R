# simulate under true model on observational grid
read_csv("data/processed/obs_data_mintp_winter.csv") %>% dplyr::select(Long.projected, Lat.projected) %>% unique %>% write_csv("data/processed/obs_grid_simulated_on.csv")
rm(list=ls())

simulate = function(marg_mod, thresh_qnt){

  library(tidyverse)
  library(doParallel)
  library(foreach)
  
  file_name = paste0("output/simulations/simulations_on_obs_grid/true/thresh_qnt_", thresh_qnt,"_marg_mod_",marg_mod,"_run_")
  fit = readRDS(paste0("output/rparp_", marg_mod, "_thresh_qnt_", thresh_qnt))$par
  
  
  alpha_var <<- fit[1]
  beta_var <<- fit[2]
  
  
  # --- marginal model names
  locs_to_pred = as.data.frame(read_csv("data/processed/obs_grid_simulated_on.csv"))
  
  variogram_model = function(h){
    h = sqrt(norm(h,type = "2")^2)
    nu=0.1
    res = rep(NA, length(h))
    res[h == 0] = 0
    res[h != 0] = alpha_var*(1 - ((((h[h != 0] /beta_var)^nu) * besselK(x = (h[h != 0] / beta_var), nu = nu))/(gamma(nu)*2^(nu - 1))))
    res
  }

  for(i in seq(1, 100)){
    print(i)
    nCores <- 1
    cl <- parallel::makeCluster(nCores)
    clusterSetRNGStream(cl)
    
    # this returns samples in the unit frechet margin
    simulation_score <- mvPot::simulPareto(n = 1000,
                                           loc = locs_to_pred,
                                           vario = variogram_model,
                                           nCores = nCores,
                                           cl = cl)
    parallel::stopCluster(cl)
    
    simulation_score %>%
      saveRDS(paste0(file_name,i))
  }
}

job::job({simulate(marg_mod = 'E1', 0.96)})
job::job({simulate(marg_mod = 'E2', 0.96)})
job::job({simulate(marg_mod = 'E2', 0.95)})
job::job({simulate(marg_mod = 'E2', 0.97)})

