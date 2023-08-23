# scale up simulations using extreme temps on frechet scale
rm(list=ls())
library(tidyverse)

scale_true_sim = function(marg_mods, yrs, tmps, var_qnt, thresh_qnt){
  for(marg_mod in marg_mods){
    for(yr in yrs){
      for(tmp in tmps){
        
        extreme_temp_frechet = read_csv(paste0("output/obs_sites_extreme_temps_frechet_scale_thresh_qnt_",thresh_qnt,"_model_",marg_mod,"_qnt_",var_qnt,".csv"))
      
        data_sets = list.files('output/simulations/simulations_on_obs_grid/true/')
        data_sets = data_sets[which(grepl(paste0("thresh_qnt_",thresh_qnt, "_"), data_sets) & grepl(marg_mod, data_sets))]
        
        
        obs_sites = read_csv("data/processed/obs_data_mintp_winter.csv") %>%
          dplyr::select(Station, 
                        scales_95_logged, threshold_95, 
                        scales_96_logged, threshold_96,
                        scales_97_logged, threshold_97,
                        Long.projected, Lat.projected, dist_sea_logged) %>%
          unique() 
        
        grid_simulated = as.data.frame(read_csv("data/processed/obs_grid_simulated_on.csv")) %>%
          left_join(obs_sites)
        
        
        
        if(thresh_qnt == 0.95){
          thresh_ex = read_csv(paste0("data/processed/thresh_exceedance_lambda_mintp_",var_qnt,".csv")) %>%
            dplyr::select(Station, year, thresh_exceedance_95_coast)  %>%
            filter(year == yr) %>%
            unique()

          this_grid = grid_simulated %>%
            left_join(thresh_ex) %>%
            left_join(extreme_temp_frechet %>% filter(temp == tmp, year == yr) %>% unique())
          
          this_grid$scales_logged = this_grid$scales_95_logged
          this_grid$thresh_exceedance = this_grid$thresh_exceedance_95_coast
          this_grid$threshold = this_grid$threshold_tn_l_o_95_w_coast
          
        }else if(thresh_qnt == 0.96){
          
          thresh_ex = read_csv(paste0("data/processed/thresh_exceedance_lambda_mintp_",var_qnt,".csv")) %>%
            dplyr::select(Station, year, thresh_exceedance_96_coast)  %>%
            filter(year == yr) %>%
            unique()
          
          this_grid = grid_simulated %>%
            left_join(thresh_ex) %>%
            left_join(extreme_temp_frechet %>% filter(temp == tmp, year == yr) %>% unique())
          
          this_grid$scales_logged = this_grid$scales_96_logged
          this_grid$thresh_exceedance = this_grid$thresh_exceedance_96_coast
          this_grid$threshold = this_grid$threshold_tn_l_o_96_w_coast
          
          
        }else if(thresh_qnt == 0.97){
          thresh_ex = read_csv(paste0("data/processed/thresh_exceedance_lambda_mintp_",var_qnt,".csv")) %>%
            dplyr::select(Station, year, thresh_exceedance_97_coast)  %>%
            filter(year == yr) %>%
            unique()
          
          this_grid = grid_simulated %>%
            left_join(thresh_ex) %>%
            left_join(extreme_temp_frechet %>% filter(temp == tmp, year == yr) %>% unique())
          
          this_grid$scales_logged = this_grid$scales_97_logged
          this_grid$thresh_exceedance = this_grid$thresh_exceedance_97_coast
          this_grid$threshold = this_grid$threshold_tn_l_o_97_w_coast
          
        }else{
          asfdgsafadafdgf()
        }
        
        

        my_simulations_extremes = c()
        for(i in seq(1, 100)){ 
          my_simulations_extremes = c(my_simulations_extremes, readRDS(paste0("output/simulations/simulations_on_obs_grid/true/thresh_qnt_", thresh_qnt,"_marg_mod_",marg_mod,"_run_",i)))
        }
        
        my_simulations_standardised = list()
        for (i in 1:length(my_simulations_extremes)) {
          this_cost = mean(my_simulations_extremes[[i]])
          my_simulations_standardised[[i]] = my_simulations_extremes[[i]]/this_cost
        }
        
        # --- maximum frechet margin at each of my simulation sites
        max_at_each_site = c()
        for(s in seq(length(my_simulations_standardised[[1]]))){
          max_at_each_site = c(max_at_each_site, lapply(my_simulations_standardised, "[[", s) %>% unlist %>% max)
        }
        
        this_grid$scaler = this_grid$frechet_value/max_at_each_site
        this_grid$scaler[this_grid$scaler == -Inf] = Inf
        br = min(this_grid$scaler)
        
        my_r = evd::rgpd(n=length(my_simulations_extremes), 1,1,1)
        # --- scale simulations back up
        my_simulations_rescaled = my_simulations_standardised
        for (i in 1:length(my_simulations_standardised)) {
          my_simulations_rescaled[[i]] = (my_r[i] * br)*my_simulations_standardised[[i]]
        }
        
        my_simulations_rescaled %>% saveRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/true/thresh_qnt_", thresh_qnt,"_marg_mod_",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp,"_qnt_", var_qnt))
        
      }
    }
  }
}

# 4 mins
# job::job({scale_true_sim(c("E2"), c(1950, 2020, 2022),  c(8,9,10), var_qnt = 0.5, thresh_qnt = 0.96)})
# job::job({scale_true_sim(c("E2"), c(1950, 2020, 2022),  c(8,9,10), var_qnt = 0.1, thresh_qnt = 0.96)})
# job::job({scale_true_sim(c("E2"), c(1950, 2020, 2022),  c(8,9,10), var_qnt = 0.9, thresh_qnt = 0.96)})

# job::job({scale_true_sim(c("E2"), c(1950, 2020, 2022),  c(8,9,10), var_qnt = 0.5, thresh_qnt = 0.97)})
# job::job({scale_true_sim(c("E2"), c(1950, 2020, 2022),  c(8,9,10), var_qnt = 0.5, thresh_qnt = 0.95)})
# job::job({scale_true_sim(c("E1"), c(1950, 2020, 2022),  c(8,9,10), var_qnt = 0.5, thresh_qnt = 0.96)})

# job::job({scale_true_sim(c("E2"), c(1950, 2020, 2022),  c(6,7,11,12), var_qnt = 0.5, thresh_qnt = 0.96)})
# job::job({scale_true_sim(c("E2"), c(1950, 2020, 2022),  c(6,7,11,12), var_qnt = 0.1, thresh_qnt = 0.96)})
# job::job({scale_true_sim(c("E2"), c(1950, 2020, 2022),  c(6,7,11,12), var_qnt = 0.9, thresh_qnt = 0.96)})
# job::job({scale_true_sim(c("E2"), c(1950, 2020, 2022),  c(6,7,11,12), var_qnt = 0.5, thresh_qnt = 0.97)})
# job::job({scale_true_sim(c("E2"), c(1950, 2020, 2022),  c(6,7,11,12), var_qnt = 0.5, thresh_qnt = 0.95)})
# job::job({scale_true_sim(c("E1"), c(1950, 2020, 2022),  c(6,7,11,12), var_qnt = 0.5, thresh_qnt = 0.96)})



scale_bts_sim = function(bts_range, marg_mods, yrs, tmps, var_qnt, thresh_qnt){  
  for(marg_mod in marg_mods){
    for(yr in yrs){
      for(tmp in tmps){
        for(bts_num in bts_range){

          if(!file.exists(paste0("output/simulations/rescaled_simulations_on_obs_grid/bts/",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp,"_qnt_", var_qnt,"_bts_", bts_num))){
            
            
            extreme_temp_frechet = vroom::vroom(paste0("output/bootstrapped_obs_sites_extreme_temps_frechet_scale_thresh_qnt_",thresh_qnt,"marg_mod",marg_mod,"_qnt_",var_qnt,".csv"),
                                                col_names = c('bts', 'year', 'Station', 'temp', 'frechet_value')) %>%
              filter(bts == bts_num) %>%
              dplyr::select(-bts) %>%
              filter(year == yr)
            
            data_sets = list.files('output/simulations/simulations_on_obs_grid/bts/')
            data_sets = data_sets[which(grepl(paste0("thresh_qnt_",thresh_qnt, "_"), data_sets) & grepl(paste0(marg_mod,"_bootstrap_",bts_num,"_"), data_sets))]

            obs_sites = read_csv("data/processed/obs_data_mintp_winter.csv") %>%
              dplyr::select(Station, 
                            scales_95_logged, threshold_95, 
                            scales_96_logged, threshold_96,
                            scales_97_logged, threshold_97,
                            Long.projected, Lat.projected, dist_sea_logged) %>%
              unique() 
            
            grid_simulated = as.data.frame(read_csv("data/processed/obs_grid_simulated_on.csv")) %>%
              left_join(obs_sites)
            
            
            
            if(thresh_qnt == 0.95){

              thresh_ex = read_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", marg_mod, "_bts_",bts_num,"_monthly_clim_covar_",var_qnt,".csv"),
                                   col_names = c("bts", "Station", "year", "thresh_exceedance_95_coast")) %>%
                dplyr::select(-bts) %>% 
                unique() %>%
                filter(year == yr) %>%
                unique()
              
              this_grid = grid_simulated %>%
                left_join(thresh_ex) %>%
                left_join(extreme_temp_frechet %>% filter(temp == tmp, year == yr) %>% unique())
              
              this_grid$scales_logged = this_grid$scales_95_logged
              this_grid$thresh_exceedance = this_grid$thresh_exceedance_95_coast
              this_grid$threshold = this_grid$threshold_tn_l_o_95_w_coast
              
            }else if(thresh_qnt == 0.96){
              
              thresh_ex = read_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", marg_mod, "_bts_",bts_num,"_monthly_clim_covar_",var_qnt,".csv"),
                                   col_names = c("bts", "Station", "year", "thresh_exceedance_96_coast")) %>%
                dplyr::select(-bts) %>% 
                unique() %>%
                filter(year == yr) %>%
                unique()
              
              this_grid = grid_simulated %>%
                left_join(thresh_ex) %>%
                left_join(extreme_temp_frechet %>% filter(temp == tmp, year == yr) %>% unique())
              
              this_grid$scales_logged = this_grid$scales_96_logged
              this_grid$thresh_exceedance = this_grid$thresh_exceedance_96_coast
              this_grid$threshold = this_grid$threshold_tn_l_o_96_w_coast
              
              
            }else if(thresh_qnt == 0.97){
              thresh_ex = read_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", marg_mod, "_bts_",bts_num,"_monthly_clim_covar_",var_qnt,".csv"),
                                   col_names = c("bts", "Station", "year", "thresh_exceedance_97_coast")) %>%
                dplyr::select(-bts) %>% 
                unique() %>%
                filter(year == yr) %>%
                unique()
              
              this_grid = grid_simulated %>%
                left_join(thresh_ex) %>%
                left_join(extreme_temp_frechet %>% filter(temp == tmp, year == yr) %>% unique())
              
              this_grid$scales_logged = this_grid$scales_97_logged
              this_grid$thresh_exceedance = this_grid$thresh_exceedance_97_coast
              this_grid$threshold = this_grid$threshold_tn_l_o_97_w_coast
              
            }else{
              asfdgsafadafdgf()
            }
            
            
            if(file.exists(paste0("output/simulations/simulations_on_obs_grid/bts/thresh_qnt_",thresh_qnt,"_marg_mod_",marg_mod,"_bootstrap_",bts_num,"_run_1"))){
              
              my_simulations_extremes = c()
              for(i in seq(1, 100)){ 
                my_simulations_extremes = c(my_simulations_extremes, readRDS(paste0("output/simulations/simulations_on_obs_grid/bts/thresh_qnt_",thresh_qnt,"_marg_mod_",marg_mod,"_bootstrap_",bts_num,"_run_",i)))
              }
              
              my_simulations_standardised = list()
              for (i in 1:length(my_simulations_extremes)) {
                this_cost = mean(my_simulations_extremes[[i]])
                my_simulations_standardised[[i]] = my_simulations_extremes[[i]]/this_cost
              }
              
              # --- maximum frechet margin at each of my simulation sites
              max_at_each_site = c()
              for(s in seq(length(my_simulations_standardised[[1]]))){
                max_at_each_site = c(max_at_each_site, lapply(my_simulations_standardised, "[[", s) %>% unlist %>% max)
              }
              
              this_grid$scaler = this_grid$frechet_value/max_at_each_site
              this_grid$scaler[this_grid$scaler == -Inf] = Inf
              br = min(this_grid$scaler)
              
              my_r = evd::rgpd(n=length(my_simulations_extremes), 1,1,1)
              
              # --- scale simulations back up
              my_simulations_rescaled = my_simulations_standardised
              for (i in 1:length(my_simulations_standardised)) {
                #my_simulations_rescaled[[i]] = (my_r[i] * this_grid$scaler)*my_simulations_standardised[[i]]
                my_simulations_rescaled[[i]] = (my_r[i] * br)*my_simulations_standardised[[i]]
              }
              
              my_simulations_rescaled %>% saveRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/bts/thresh_qnt_", thresh_qnt,"_marg_mod_",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp,"_qnt_", var_qnt,"_bts_", bts_num))
              
            }
          }
        }
      }
    }
  }
}


# job::job({scale_bts_sim(seq(200), c("E2"), c(1950, 2022),  c(8,9,10), var_qnt = 0.5, thresh_qnt = 0.96)})
# job::job({scale_bts_sim(seq(200),c("E2"), c(1950, 2022),  c(8,9,10), var_qnt = 0.1, thresh_qnt = 0.96)})
# job::job({scale_bts_sim(seq(200),c("E2"), c(1950, 2022),  c(8,9,10), var_qnt = 0.9, thresh_qnt = 0.96)})
# job::job({scale_bts_sim(seq(200),c("E2"), c(1950, 2022),  c(8,9,10), var_qnt = 0.5, thresh_qnt = 0.97)})
# job::job({scale_bts_sim(seq(200),c("E2"), c(1950, 2022),  c(8,9,10), var_qnt = 0.5, thresh_qnt = 0.95)})
# job::job({scale_bts_sim(seq(200),c("E1"), c(1950, 2022),  c(8,9,10), var_qnt = 0.5, thresh_qnt = 0.96)})
# 
# job::job({scale_bts_sim(seq(200),c("E2"), c(1950,  2022),  c(6,7,11,12), var_qnt = 0.5, thresh_qnt = 0.96)})
# job::job({scale_bts_sim(seq(200),c("E2"), c(1950,  2022),  c(6,7,11,12), var_qnt = 0.1, thresh_qnt = 0.96)})
# job::job({scale_bts_sim(seq(200),c("E2"), c(1950,  2022),  c(6,7,11,12), var_qnt = 0.9, thresh_qnt = 0.96)})
# job::job({scale_bts_sim(seq(200),c("E2"), c(1950,  2022),  c(6,7,11,12), var_qnt = 0.5, thresh_qnt = 0.97)})
# job::job({scale_bts_sim(seq(200),c("E2"), c(1950,  2022),  c(6,7,11,12), var_qnt = 0.5, thresh_qnt = 0.95)})
# job::job({scale_bts_sim(seq(200),c("E1"), c(1950,  2022),  c(6,7,11,12), var_qnt = 0.5, thresh_qnt = 0.96)})