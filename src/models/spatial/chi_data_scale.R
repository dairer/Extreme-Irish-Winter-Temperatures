
# this function calcualtes data scale chi on scaled simulations
library(tidyverse)

calc_chi_true = function(marg_mods, yrs, tmps, var_qnt, thresh_qnt){
  for(marg_mod in marg_mods){
    for(yr in yrs){
      for(tmp in tmps){
        # for(var_qnt in var_qnts){
          set.seed(123456)
          # ¸  ¸

          sites = readRDS("output/site_pairs")
        
          obs_sites = read_csv("data/processed/obs_data_mintp_winter.csv") %>%
            dplyr::select(Station, Long.projected, Lat.projected, dist_sea_logged) %>%
            unique()
          
          grid_simulated = as.data.frame(read_csv("data/processed/obs_grid_simulated_on.csv")) %>%
            left_join(obs_sites) %>% as.tibble()

          sites = sites %>% filter(site_1 %in% grid_simulated$Station)
          sites = sites %>% filter(site_2 %in% grid_simulated$Station)

          frechet_val = grid_simulated %>%
            left_join(read_csv(paste0("output/obs_sites_extreme_temps_frechet_scale_thresh_qnt_",thresh_qnt,"_model_",marg_mod,"_qnt_",var_qnt,".csv")) %>% 
                        filter(temp == tmp, year == yr) %>%
                        unique()) %>% 
            pull(frechet_value)
          
          if(file.exists(paste0("output/simulations/rescaled_simulations_on_obs_grid/true/thresh_qnt_", thresh_qnt,"_marg_mod_",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp,"_qnt_", var_qnt))){
            my_simulations = readRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/true/thresh_qnt_", thresh_qnt,"_marg_mod_",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp,"_qnt_", var_qnt))            
            
            for(s in seq(nrow(sites))){
              
              print(s)
              if((s %% 100) == 0){
                print(s)
              }

              id1 = sites[s,]$site_1
              id2 = sites[s,]$site_2
              
              # ---- to get all sims at loc with id v1 ---> unlist(lapply(my_simulations, "[[", sites[s,]$V1))
              id1_exceeds = (unlist(lapply(my_simulations, "[[", which(grid_simulated$Station == id1))) > frechet_val[which(grid_simulated$Station == id1)])
              id2_exceeds = (unlist(lapply(my_simulations, "[[", which(grid_simulated$Station == id2))) > frechet_val[which(grid_simulated$Station == id2)])
              
              tibble(sites[s,] %>% mutate(chi =sum(id1_exceeds & id2_exceeds)/sum(id1_exceeds) )) %>%
                write_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_thresh_qnt_", thresh_qnt,"model_",marg_mod,"_yr_",yr, "_conditioned_on_",tmp,"_qnt_", var_qnt, ".csv"),append = T)
            }
          }
        # }
      }
    }
  }
}


# # 15 mins 
job::job({calc_chi_true(marg_mods = "m22", yrs = 1950, tmps = 10, var_qnt = 0.1, thresh_qnt = 0.96)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 1950, tmps = 10, var_qnt = 0.5, thresh_qnt = 0.96)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 1950, tmps = 10, var_qnt = 0.9, thresh_qnt = 0.96)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 2020, tmps = 10, var_qnt = 0.1, thresh_qnt = 0.96)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 2020, tmps = 10, var_qnt = 0.5, thresh_qnt = 0.96)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 2020, tmps = 10, var_qnt = 0.9, thresh_qnt = 0.96)})

job::job({calc_chi_true(marg_mods = "m22", yrs = 1950, tmps = 9, var_qnt = 0.1, thresh_qnt = 0.96)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 1950, tmps = 9, var_qnt = 0.5, thresh_qnt = 0.96)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 1950, tmps = 9, var_qnt = 0.9, thresh_qnt = 0.96)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 2020, tmps = 9, var_qnt = 0.1, thresh_qnt = 0.96)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 2020, tmps = 9, var_qnt = 0.5, thresh_qnt = 0.96)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 2020, tmps = 9, var_qnt = 0.9, thresh_qnt = 0.96)})

job::job({calc_chi_true(marg_mods = "m22", yrs = 1950, tmps = 8, var_qnt = 0.1, thresh_qnt = 0.96)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 1950, tmps = 8, var_qnt = 0.5, thresh_qnt = 0.96)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 1950, tmps = 8, var_qnt = 0.9, thresh_qnt = 0.96)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 2020, tmps = 8, var_qnt = 0.1, thresh_qnt = 0.96)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 2020, tmps = 8, var_qnt = 0.5, thresh_qnt = 0.96)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 2020, tmps = 8, var_qnt = 0.9, thresh_qnt = 0.96)})

# 36 mins
job::job({calc_chi_true(marg_mods = "m21", yrs = 1950, tmps = c(8, 9, 10), var_qnt = 0.5, thresh_qnt = 0.96)})
job::job({calc_chi_true(marg_mods = "m21", yrs = 2020, tmps = c(8, 9, 10), var_qnt = 0.5, thresh_qnt = 0.96)})

job::job({calc_chi_true(marg_mods = "m22", yrs = 1950, tmps = c(8, 9, 10), var_qnt = 0.5, thresh_qnt = 0.95)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 2020, tmps = c(8, 9, 10), var_qnt = 0.5, thresh_qnt = 0.95)})

job::job({calc_chi_true(marg_mods = "m22", yrs = 1950, tmps = c(8, 9, 10), var_qnt = 0.5, thresh_qnt = 0.97)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 2020, tmps = c(8, 9, 10), var_qnt = 0.5, thresh_qnt = 0.97)})

# 36 mins
job::job({calc_chi_true(marg_mods = "m22", yrs = 2022, tmps = c(8, 9, 10), var_qnt = 0.1, thresh_qnt = 0.96)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 2022, tmps = c(8, 9, 10), var_qnt = 0.5, thresh_qnt = 0.96)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 2022, tmps = c(8, 9, 10), var_qnt = 0.9, thresh_qnt = 0.96)})

job::job({calc_chi_true(marg_mods = "m21", yrs = 2022, tmps = c(8, 9, 10), var_qnt = 0.5, thresh_qnt = 0.96)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 2022, tmps = c(8, 9, 10), var_qnt = 0.5, thresh_qnt = 0.95)})
job::job({calc_chi_true(marg_mods = "m22", yrs = 2022, tmps = c(8, 9, 10), var_qnt = 0.5, thresh_qnt = 0.97)})





calc_chi_bts =function(marg_mods, yrs, tmps, var_qnt, thresh_qnt, bts_range){
  for(marg_mod in marg_mods){
    for(yr in yrs){
      for(tmp in tmps){
        for(bts_num in bts_range){

        num_samples = 1000
        # set.seed(123456)
        
        sites = readRDS("output/site_pairs") %>% sample_n(num_samples)
        obs_sites = read_csv("data/processed/obs_data_mintp_winter.csv") %>%
          dplyr::select(Station, Long.projected, Lat.projected, dist_sea_logged) %>%
          unique()

        grid_simulated = as.data.frame(read_csv("data/processed/obs_grid_simulated_on.csv")) %>%
          left_join(obs_sites) %>% as.tibble()
        
        sites = sites %>% filter(site_1 %in% grid_simulated$Station)
        sites = sites %>% filter(site_2 %in% grid_simulated$Station)
    
        frechet_val = grid_simulated %>%
          left_join(read_csv(paste0("output/bootstrapped_obs_sites_extreme_temps_frechet_scale_thresh_qnt_",thresh_qnt,"marg_mod",marg_mod,"_qnt_",var_qnt,".csv"),
                             col_names = c('bts', 'year', 'Station', 'temp', 'frechet_value')) %>% 
                      filter(bts == bts_num, temp == tmp, year == yr) %>%
                      unique()) %>% 
          pull(frechet_value)

        if(file.exists(paste0("output/simulations/rescaled_simulations_on_obs_grid/bts/thresh_qnt_", thresh_qnt,"_marg_mod_",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp,"_qnt_", var_qnt,"_bts_", bts_num))){
          my_simulations = readRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/bts/thresh_qnt_", thresh_qnt,"_marg_mod_",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp,"_qnt_", var_qnt,"_bts_", bts_num))
          # my_simulations = my_simulations[sample(seq(100000), 50000, replace = F)]

          for(s in seq(nrow(sites))){

            id1 = sites[s,]$site_1
            id2 = sites[s,]$site_2
            
            # ---- to get all sims at loc with id v1 ---> unlist(lapply(my_simulations, "[[", sites[s,]$V1))
            id1_exceeds = (unlist(lapply(my_simulations, "[[", which(grid_simulated$Station == id1))) > frechet_val[which(grid_simulated$Station == id1)])
            id2_exceeds = (unlist(lapply(my_simulations, "[[", which(grid_simulated$Station == id2))) > frechet_val[which(grid_simulated$Station == id2)])
            
            tibble(sites[s,] %>% mutate(bts = bts_num, chi =sum(id1_exceeds & id2_exceeds)/sum(id1_exceeds) )) %>%
              write_csv(paste0("output/simulations/simulation_summary/new_bootstrap_chi_data_scale_model_thresh_qnt_", thresh_qnt,"model_",marg_mod,"_yr_",yr, "_conditioned_on_",tmp,"_qnt_", var_qnt, ".csv"),append = T)
          }
          }
        }
      }
    }
  }
}



# 6.5 hours
job::job({calc_chi_bts(marg_mods = "m22", yrs = c(2022, 1950), tmps = c(8, 9, 10), var_qnt = 0.5, thresh_qnt = c(0.96), seq(1, 200))})
job::job({calc_chi_bts(marg_mods = "m22", yrs = c(2022, 1950), tmps = c(8, 9, 10), var_qnt = 0.1, thresh_qnt = c(0.96), seq(1, 200))})
job::job({calc_chi_bts(marg_mods = "m22", yrs = c(2022, 1950), tmps = c(8, 9, 10), var_qnt = 0.9, thresh_qnt = c(0.96), seq(1, 200))})
job::job({calc_chi_bts(marg_mods = "m22", yrs = c(2022, 1950), tmps = c(8, 9, 10), var_qnt = 0.5, thresh_qnt = c(0.95), seq(1, 200))})
job::job({calc_chi_bts(marg_mods = "m22", yrs = c(2022, 1950), tmps = c(8, 9, 10), var_qnt = 0.5, thresh_qnt = c(0.97), seq(1, 200))})



# # ------ PLOT MODELS
marg_mod = 'm22'
thresh_qnt = 0.96

chi_true = rbind(
  read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_",2022, "_conditioned_on_",8,"_qnt_", 0.1, ".csv"),
           col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "-8°C", year = '2022', qnt = 0.1),
  read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_",2022, "_conditioned_on_",9,"_qnt_", 0.1, ".csv"),
           col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "-9°C", year = '2022', qnt = 0.1),
  read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_",2022, "_conditioned_on_",10,"_qnt_", 0.1, ".csv"),
           col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "-10°C", year = '2022', qnt = 0.1),
  read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_",1950, "_conditioned_on_",8,"_qnt_", 0.1, ".csv"),
           col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "-8°C", year = '1950', qnt = 0.1),
  read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_",1950, "_conditioned_on_",9,"_qnt_", 0.1, ".csv"),
           col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "-9°C", year = '1950', qnt = 0.1),
  read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_",1950, "_conditioned_on_",10,"_qnt_", 0.1, ".csv"),
           col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "-10°C", year = '1950', qnt = 0.1),
  read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_",2022, "_conditioned_on_",8,"_qnt_", 0.5, ".csv"),
           col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "-8°C", year = '2022', qnt = 0.5),
  read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_",2022, "_conditioned_on_",9,"_qnt_", 0.5, ".csv"),
           col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "-9°C", year = '2022', qnt = 0.5),
  read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_",2022, "_conditioned_on_",10,"_qnt_", 0.5, ".csv"),
           col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "-10°C", year = '2022', qnt = 0.5),
  read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_",1950, "_conditioned_on_",8,"_qnt_", 0.5, ".csv"),
           col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "-8°C", year = '1950', qnt = 0.5),
  read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_",1950, "_conditioned_on_",9,"_qnt_", 0.5, ".csv"),
           col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "-9°C", year = '1950', qnt = 0.5),
  read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_",1950, "_conditioned_on_",10,"_qnt_", 0.5, ".csv"),
           col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "-10°C", year = '1950', qnt = 0.5),
  read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_",2022, "_conditioned_on_",8,"_qnt_", 0.9, ".csv"),
           col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "-8°C", year = '2022', qnt = 0.9),
  read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_",2022, "_conditioned_on_",9,"_qnt_", 0.9, ".csv"),
           col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "-9°C", year = '2022', qnt = 0.9),
  read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_",2022, "_conditioned_on_",10,"_qnt_", 0.9, ".csv"),
           col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "-10°C", year = '2022', qnt = 0.9),
  read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_",1950, "_conditioned_on_",8,"_qnt_", 0.9, ".csv"),
           col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "-8°C", year = '1950', qnt = 0.9),
  read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_",1950, "_conditioned_on_",9,"_qnt_", 0.9, ".csv"),
           col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "-9°C", year = '1950', qnt = 0.9),
  read_csv(paste0("output/simulations/simulation_summary/chi_data_scale_clim_grid_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_",1950, "_conditioned_on_",10,"_qnt_", 0.9, ".csv"),
           col_names = c('s1', 's2', 'distance', 'chi')) %>% mutate(temp = "-10°C", year = '1950', qnt = 0.9)
  )






chi_bts = rbind(read_csv(paste0("output/simulations/simulation_summary/new_bootstrap_chi_data_scale_model_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_", 1950, "_conditioned_on_",8,"_qnt_", 0.1, ".csv"),
                         col_names = c('s1', 's2', 'distance', 'bts', 'chi')) %>% mutate(temp = "-8°C", year = '1950', qnt = 0.1),
                read_csv(paste0("output/simulations/simulation_summary/new_bootstrap_chi_data_scale_model_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_", 1950, "_conditioned_on_",9,"_qnt_", 0.1, ".csv"),
                         col_names = c('s1', 's2', 'distance', 'bts', 'chi')) %>% mutate(temp = "-9°C", year = '1950', qnt = 0.1),
                read_csv(paste0("output/simulations/simulation_summary/new_bootstrap_chi_data_scale_model_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_", 1950, "_conditioned_on_",10,"_qnt_", 0.1, ".csv"),
                         col_names = c('s1', 's2', 'distance', 'bts', 'chi')) %>% mutate(temp = "-10°C", year = '1950', qnt = 0.1),
                read_csv(paste0("output/simulations/simulation_summary/new_bootstrap_chi_data_scale_model_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_", 2022, "_conditioned_on_",8,"_qnt_", 0.1, ".csv"),
                         col_names = c('s1', 's2', 'distance', 'bts', 'chi')) %>% mutate(temp = "-8°C", year = '2022', qnt = 0.1),
                read_csv(paste0("output/simulations/simulation_summary/new_bootstrap_chi_data_scale_model_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_", 2022, "_conditioned_on_",9,"_qnt_", 0.1, ".csv"),
                         col_names = c('s1', 's2', 'distance', 'bts', 'chi')) %>% mutate(temp = "-9°C", year = '2022', qnt = 0.1),
                read_csv(paste0("output/simulations/simulation_summary/new_bootstrap_chi_data_scale_model_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_", 2022, "_conditioned_on_",10,"_qnt_", 0.1, ".csv"),
                         col_names = c('s1', 's2', 'distance', 'bts', 'chi')) %>% mutate(temp = "-10°C", year = '2022', qnt = 0.1),
                read_csv(paste0("output/simulations/simulation_summary/new_bootstrap_chi_data_scale_model_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_", 1950, "_conditioned_on_",8,"_qnt_", 0.5, ".csv"),
                         col_names = c('s1', 's2', 'distance', 'bts', 'chi')) %>% mutate(temp = "-8°C", year = '1950', qnt = 0.5),
                read_csv(paste0("output/simulations/simulation_summary/new_bootstrap_chi_data_scale_model_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_", 1950, "_conditioned_on_",9,"_qnt_", 0.5, ".csv"),
                         col_names = c('s1', 's2', 'distance', 'bts', 'chi')) %>% mutate(temp = "-9°C", year = '1950', qnt = 0.5),
                read_csv(paste0("output/simulations/simulation_summary/new_bootstrap_chi_data_scale_model_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_", 1950, "_conditioned_on_",10,"_qnt_", 0.5, ".csv"),
                         col_names = c('s1', 's2', 'distance', 'bts', 'chi')) %>% mutate(temp = "-10°C", year = '1950', qnt = 0.5),
                read_csv(paste0("output/simulations/simulation_summary/new_bootstrap_chi_data_scale_model_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_", 2022, "_conditioned_on_",8,"_qnt_", 0.5, ".csv"),
                         col_names = c('s1', 's2', 'distance', 'bts', 'chi')) %>% mutate(temp = "-8°C", year = '2022', qnt = 0.5),
                read_csv(paste0("output/simulations/simulation_summary/new_bootstrap_chi_data_scale_model_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_", 2022, "_conditioned_on_",9,"_qnt_", 0.5, ".csv"),
                         col_names = c('s1', 's2', 'distance', 'bts', 'chi')) %>% mutate(temp = "-9°C", year = '2022', qnt = 0.5),
                read_csv(paste0("output/simulations/simulation_summary/new_bootstrap_chi_data_scale_model_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_", 2022, "_conditioned_on_",10,"_qnt_", 0.5, ".csv"),
                         col_names = c('s1', 's2', 'distance', 'bts', 'chi')) %>% mutate(temp = "-10°C", year = '2022', qnt = 0.5),
                read_csv(paste0("output/simulations/simulation_summary/new_bootstrap_chi_data_scale_model_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_", 1950, "_conditioned_on_",8,"_qnt_", 0.9, ".csv"),
                         col_names = c('s1', 's2', 'distance', 'bts', 'chi')) %>% mutate(temp = "-8°C", year = '1950', qnt = 0.9),
                read_csv(paste0("output/simulations/simulation_summary/new_bootstrap_chi_data_scale_model_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_", 1950, "_conditioned_on_",9,"_qnt_", 0.9, ".csv"),
                         col_names = c('s1', 's2', 'distance', 'bts', 'chi')) %>% mutate(temp = "-9°C", year = '1950', qnt = 0.9),
                read_csv(paste0("output/simulations/simulation_summary/new_bootstrap_chi_data_scale_model_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_", 1950, "_conditioned_on_",10,"_qnt_", 0.9, ".csv"),
                         col_names = c('s1', 's2', 'distance', 'bts', 'chi')) %>% mutate(temp = "-10°C", year = '1950', qnt = 0.9),
                read_csv(paste0("output/simulations/simulation_summary/new_bootstrap_chi_data_scale_model_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_", 2022, "_conditioned_on_",8,"_qnt_", 0.9, ".csv"),
                         col_names = c('s1', 's2', 'distance', 'bts', 'chi')) %>% mutate(temp = "-8°C", year = '2022', qnt = 0.9),
                read_csv(paste0("output/simulations/simulation_summary/new_bootstrap_chi_data_scale_model_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_", 2022, "_conditioned_on_",9,"_qnt_", 0.9, ".csv"),
                         col_names = c('s1', 's2', 'distance', 'bts', 'chi')) %>% mutate(temp = "-9°C", year = '2022', qnt = 0.9),
                read_csv(paste0("output/simulations/simulation_summary/new_bootstrap_chi_data_scale_model_thresh_qnt_", thresh_qnt,"model_",'m22',"_yr_", 2022, "_conditioned_on_",10,"_qnt_", 0.9, ".csv"),
                         col_names = c('s1', 's2', 'distance', 'bts', 'chi')) %>% mutate(temp = "-10°C", year = '2022', qnt = 0.9)
                )
sum(chi_bts$chi == 0)
chi_bts = chi_bts %>%
  drop_na() 
chi_bts$bts %>% unique()
chi_bts$chi[is.na(chi_bts$chi)] = 0
chi_true$chi[is.na(chi_true$chi)] = 0


chi_bts_summary = chi_bts %>%
  mutate(dist_bin = cut(distance, breaks=seq(0, 4, length.out = 20))) %>%
  group_by(dist_bin) %>%
  summarise(bts, year, temp, chi, qnt, distance = mean(distance))  %>%
  group_by(distance, year, temp, bts, qnt) %>%
  drop_na() %>%
  summarise(mn = mean(chi)) %>%
  group_by(distance, year, temp, qnt)  %>%
  summarise(upper = quantile(mn, 0.975),
            lower = quantile(mn, 0.025))

chi_true_summary = chi_true %>%
  mutate(dist_bin = cut(distance, breaks=seq(0, 4, length.out = 20))) %>%
  group_by(dist_bin) %>%
  summarise(year, temp, chi,qnt,  distance = mean(distance))  %>%
  group_by(distance, year, qnt, temp) %>%
  drop_na() %>%
  summarise(mn = mean(chi))



prob_observing_T = rbind(read_csv(paste0("output/prob_extreme_temp_imp_samp_thresh_quant_",thresh_qnt,"_mar_mod_",marg_mod,"_var_qnt_",0.1,".csv"),
                                  col_names = c('temp', 'p_1950', 'p_2020', 'p_2022')) %>% mutate(qnt = 0.1),
                         read_csv(paste0("output/prob_extreme_temp_imp_samp_thresh_quant_",thresh_qnt,"_mar_mod_",marg_mod,"_var_qnt_",0.5,".csv"),
                                  col_names = c('temp', 'p_1950', 'p_2020', 'p_2022')) %>% mutate(qnt = 0.5),
                         read_csv(paste0("output/prob_extreme_temp_imp_samp_thresh_quant_",thresh_qnt,"_mar_mod_",marg_mod,"_var_qnt_",0.9,".csv"),
                                  col_names = c('temp', 'p_1950', 'p_2020', 'p_2022')) %>% mutate(qnt = 0.9))

prob_observing_T = prob_observing_T %>%
  filter(temp %in% c(8,9,10))


prob_observing_T = prob_observing_T%>%
  dplyr::select(temp,p_1950,p_2022, qnt) %>%
  pivot_longer(-c(temp,qnt)) %>%
  mutate(year = str_remove(name,"p_")) 

prob_observing_T = rbind(prob_observing_T %>% filter(name == 'p_1950') %>% 
                           mutate(typ = 'actual') %>% dplyr::select(-name),
                         prob_observing_T %>% filter(name == 'p_2022') %>% 
                           mutate(typ = 'actual') %>% dplyr::select(-name)) %>%
  rename(unconditional_factor = value)

prob_observing_T$temp[prob_observing_T$temp == 8]="-8°C"
prob_observing_T$temp[prob_observing_T$temp == 9]="-9°C"
prob_observing_T$temp[prob_observing_T$temp == 10]="-10°C"


chi_bts_summary = chi_bts_summary %>%
  left_join(prob_observing_T, by = c('qnt', 'temp', 'year'))
chi_true_summary = chi_true_summary %>%
  left_join(prob_observing_T, by = c('qnt', 'temp', 'year'))


bts_unconditioned = chi_bts_summary
true_unconditioned = chi_true_summary


bts_unconditioned$upper = bts_unconditioned$upper*bts_unconditioned$unconditional_factor
bts_unconditioned$lower = bts_unconditioned$lower*bts_unconditioned$unconditional_factor
true_unconditioned$mn = true_unconditioned$mn*true_unconditioned$unconditional_factor



bts_dat = rbind(chi_bts_summary %>%
                  mutate(lab = "conditioned"),
                bts_unconditioned %>%
                  mutate(lab = "unconditioned"))


true_data = rbind(chi_true_summary %>%
                    mutate(lab = "conditioned"),
                  true_unconditioned %>%
                    mutate(lab = "unconditioned"))


bts_dat$temp = factor(bts_dat$temp, levels = c("-8°C",  "-9°C", "-10°C"))
true_data$temp = factor(true_data$temp, levels = c("-8°C",  "-9°C", "-10°C"))



# ---- for text, prop of change
# qnt = 0.1
true_data %>%
  filter(distance > 0.9,
         distance <1) %>%
  filter(temp == '-8°C', qnt == '0.1', lab == 'unconditioned') %>%
  pull(mn)
#0.0069460072/0.0009389482
#7.397647

true_data %>%
  filter(distance > 0.9,
         distance <1) %>%
  filter(temp == '-9°C', qnt == '0.1', lab == 'unconditioned') %>%
  pull(mn)
#0.0027881367/0.0002030585
#13.73071

true_data %>%
  filter(distance > 0.9,
         distance <1) %>%
  filter(temp == '-10°C', qnt == '0.1', lab == 'unconditioned') %>%
  pull(mn)
#1.199187e-03/ 4.952731e-05
#24.21264

# qnt = 0.5
true_data %>%
  filter(distance > 0.9,
         distance <1) %>%
  filter(temp == '-8°C', qnt == '0.5', lab == 'unconditioned') %>%
  pull(mn)
#0.017198512/ 0.003680905
#4.67236

true_data %>%
  filter(distance > 0.9,
         distance <1) %>%
  filter(temp == '-9°C', qnt == '0.5', lab == 'unconditioned') %>%
  pull(mn)
#0.009189240/0.001329564
#6.911469

true_data %>%
  filter(distance > 0.9,
         distance <1) %>%
  filter(temp == '-10°C', qnt == '0.5', lab == 'unconditioned') %>%
  pull(mn)
#0.004583219/0.000383728
#11.94393


# qnt = 0.9
true_data %>%
  filter(distance > 0.9,
         distance <1) %>%
  filter(temp == '-8°C', qnt == '0.9', lab == 'unconditioned') %>%
  pull(mn)
#0.09147712/0.03248713
#2.815796

true_data %>%
  filter(distance > 0.9,
         distance <1) %>%
  filter(temp == '-9°C', qnt == '0.9', lab == 'unconditioned') %>%
  pull(mn)
#0.05778302/0.01743325
#3.314529

true_data %>%
  filter(distance > 0.9,
         distance <1) %>%
  filter(temp == '-10°C', qnt == '0.9', lab == 'unconditioned') %>%
  pull(mn)
#0.035550457/0.008983027
#3.957514

\
plt = gridExtra::grid.arrange(bts_dat %>%
  filter(lab == 'unconditioned', qnt == 0.1) %>%
  ggplot()+
  geom_ribbon(aes(x = distance*100, ymin = lower, ymax = upper, fill = year), alpha = 0.25)+
  geom_smooth(data = true_data %>% filter(lab == 'unconditioned', qnt == 0.1) , aes(x = distance*100, y = mn,col = year, linetype = year), se=F)+
  facet_wrap(~temp)+
  labs(x = "Distance (km)",
       y = expression(chi[o]),
       col = "Year",
       shape = "Year",
       linetype = "Year")+
  xlim(0, 375)+
  theme_minimal(12)+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.position = 'none',
        strip.text.y = element_text(size=0)),
bts_dat %>%
  filter(lab == 'unconditioned', qnt == 0.5) %>%
  ggplot()+
  geom_ribbon(aes(x = distance*100, ymin = lower, ymax = upper, fill = year), alpha = 0.25)+
  geom_smooth(data = true_data %>% filter(lab == 'unconditioned', qnt == 0.5) , aes(x = distance*100, y = mn,col = year, linetype = year), se=F)+
  facet_wrap(~temp)+
  labs(x = "Distance (km)",
       y = expression(chi[o]),
       col = "Year",
       shape = "Year",
       linetype = "Year")+
  xlim(0, 375)+
  theme_minimal(12)+
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        legend.position = 'none',
        strip.text.y = element_text(size=0)),
bts_dat %>%
  filter(lab == 'unconditioned', qnt == 0.9) %>%
  ggplot()+
  geom_ribbon(aes(x = distance*100, ymin = lower, ymax = upper, fill = year), alpha = 0.25)+
  geom_smooth(data = true_data %>% filter(lab == 'unconditioned', qnt == 0.9) , aes(x = distance*100, y = mn,col = year, linetype = year), se=F)+
  facet_wrap(~temp)+
  labs(x = "Distance (km)",
       y = expression(chi[o]),
       col = "Year",
       shape = "Year",
       linetype = "Year")+
  xlim(0, 375)+
  theme_minimal(12)+
  theme(axis.text.x = element_text(angle = 30),
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        strip.text = element_blank(),
        legend.position = 'none',
        strip.text.y = element_text(size=0)), nrow = 3, heights = c(1, 0.9, 1))


ggsave(plot = plt, paste0("output/figs/chi_data_scale_",marg_mod,".pdf"), height = 6, width =8)