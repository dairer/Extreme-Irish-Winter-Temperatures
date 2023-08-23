# simulate under bootstrapped r pareto fits

rm(list=ls())
simulate_this_mod = function(bts_nums, marg_mod, thresh_qnt){
  library(tidyverse)
  library(doParallel)
  library(foreach)

  
  locs_to_pred = as.data.frame(read_csv("data/processed/obs_grid_simulated_on.csv"))
  
  variogram_model <- function(h){
    h = sqrt(norm(h,type = "2")^2)
    nu=0.1
    res = rep(NA, length(h))
    res[h == 0] = 0
    res[h != 0] = alpha_var*(1 - ((((h[h != 0] /beta_var)^nu) * besselK(x = (h[h != 0] / beta_var), nu = nu))/(gamma(nu)*2^(nu - 1))))
    res
  }
  
  my_simulate <- function(file_name){
    print(file_name)
    for(i in seq(1,100)){
      nCores <- 2
      cl <- parallel::makeCluster(nCores)
      clusterSetRNGStream(cl)
      
      simulation_score <- mvPot::simulPareto(n = 100,
                                             loc = locs_to_pred,
                                             vario = variogram_model,
                                             nCores = nCores,
                                             cl = cl)
      parallel::stopCluster(cl)
      
      
      file.remove(paste0(file_name,i))
      simulation_score %>%
        saveRDS(paste0(file_name,i))
    }
  }

  for(x in bts_nums){
    print(x)

    stationary_models = read_csv(paste0("output/thresh_qnt_",thresh_qnt,"_bts_rpareto_fits_model_", marg_mod, "_L-BFGS-B.csv"),
                                 col_names = c('bts', "alpha", "beta"))
    
    
    this_rparp_mod = stationary_models %>% filter(bts == x)
    
    if(!is.na(this_rparp_mod$alpha)){
      alpha_var <- this_rparp_mod$alpha[nrow(this_rparp_mod)]
      beta_var <- this_rparp_mod$beta[nrow(this_rparp_mod)]
      if(length(alpha_var)>0){
        my_simulate(paste0("output/simulations/simulations_on_obs_grid/bts/thresh_qnt_",thresh_qnt,"_marg_mod_",marg_mod,"_bootstrap_",this_rparp_mod$bts[nrow(this_rparp_mod)],"_run_"))
      }
    }
  }
}


job::job({simulate_this_mod(seq(1, 200),  'E2', 0.95)}) # 2.5 hours, 2 cores
job::job({simulate_this_mod(seq(1, 200),  'E2', 0.96)}) # 2.5 hours, 2 cores
job::job({simulate_this_mod(seq(1, 200),  'E2', 0.97)}) # 2.5 hours, 2 cores
job::job({simulate_this_mod(seq(1, 200),  'E1', 0.96)}) # 2.5 hours, 2 cores
