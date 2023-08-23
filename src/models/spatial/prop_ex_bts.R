# proportion of exceedence calculated on bootstrapped data sets
rm(list=ls())
library(tidyverse)

prop_ex_bts = function(marg_mod, temp_conditioned_on, var_qnt, thresh_qnt, bts_range){
  for(bts_num in bts_range){
    obs_sites = read_csv("data/processed/obs_data_mintp_winter.csv") %>%
      dplyr::select(Station, Long.projected, Lat.projected, dist_sea_logged) %>%
      unique()
    
    grid_simulated = as.data.frame(read_csv("data/processed/obs_grid_simulated_on.csv")) %>%
      left_join(obs_sites) %>% as.tibble()
    
    
    for(tmp in temp_conditioned_on){
      yr = 2022
      
      if(file.exists(paste0("output/simulations/rescaled_simulations_on_obs_grid/bts/thresh_qnt_", thresh_qnt,"_marg_mod_",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp,"_qnt_", var_qnt,"_bts_", bts_num))){
        my_simulations = readRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/bts/thresh_qnt_", thresh_qnt,"_marg_mod_",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp,"_qnt_", var_qnt,"_bts_", bts_num))
        
        prop_ex_2022 = c()
        
        
        frechet_val = grid_simulated %>%
          left_join(read_csv(paste0("output/bootstrapped_obs_sites_extreme_temps_frechet_scale_thresh_qnt_",thresh_qnt,"marg_mod",marg_mod,"_qnt_",var_qnt,".csv"),
                             col_names = c('bts', 'year', 'Station', 'temp','frechet_value')) %>% 
                      filter(bts == bts_num, temp == tmp, year == yr) %>%
                      unique()) %>% 
          pull(frechet_value)
        
        
        frechet_val[frechet_val == -Inf] = Inf
        
        num_exceed_tmp = my_simulations %>%
          map(~{ sum(.x > frechet_val)}) %>%
          unlist()
        
        prop_ex_2022 =  mean(num_exceed_tmp[num_exceed_tmp>0]/108)
        
        
        yr = 1950
        my_simulations = readRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/bts/thresh_qnt_", thresh_qnt,"_marg_mod_",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp,"_qnt_", var_qnt,"_bts_", bts_num))
        
        frechet_val = grid_simulated %>%
          left_join(read_csv(paste0("output/bootstrapped_obs_sites_extreme_temps_frechet_scale_thresh_qnt_",thresh_qnt,"marg_mod",marg_mod,"_qnt_",var_qnt,".csv"),
                             col_names = c('bts', 'year', 'Station', 'temp','frechet_value')) %>% 
                      filter(bts == bts_num, temp == tmp, year == yr) %>%
                      unique()) %>% 
          pull(frechet_value)
        
        frechet_val[frechet_val == -Inf] = Inf
        
        num_exceed_tmp = my_simulations %>% 
          map(~{ sum(.x > frechet_val)}) %>%
          unlist()
        
        prop_ex_1950 = mean(num_exceed_tmp[num_exceed_tmp>0]/108)
        
        print("---------")
        print(prop_ex_1950)
        
        tibble(bts_num, temp = tmp, prop_ex_1950, prop_ex_2022) %>%
          write_csv(paste0("output/simulations/simulation_summary/new_bootstrapped_prop_exceedance_model_true_thresh_qnt_",thresh_qnt,"_marg_mod_",marg_mod,"_qnt_", var_qnt,".csv"), append = T)
        
      }
    }
    
  }
}

job::job({prop_ex_bts('E2', seq(6, 12), 0.1, 0.96, seq(200))})
job::job({prop_ex_bts('E2', seq(6, 12), 0.5, 0.96, seq(200))})
job::job({prop_ex_bts('E2', seq(6, 12), 0.9, 0.96, seq(200))})


# 
job::job({prop_ex_bts('E2', seq(6, 12), 0.5, 0.95, seq(200))})
job::job({prop_ex_bts('E2', seq(6, 12), 0.5, 0.97, seq(200))})
job::job({prop_ex_bts('E1', seq(6, 12), 0.5, 0.96, seq(200))})

