# get importance sample estimates of spatial RL of high levels on bootstrapped data sets
gc()
rm(list = ls())
library(tidyverse)

get_bts_imp_samp = function(marg_mod, bts_range, var_qnt, thresh_qnt){
  for(bts in bts_range){
    
    source('src/models/marginal/new_gpd_models.R')
    
    if(file.exists(paste0("output/simulations/simulations_on_obs_grid/bts/thresh_qnt_",thresh_qnt,"_marg_mod_",marg_mod,"_bootstrap_",bts,"_run_1"))){

      
      spatial_covars = read_csv("data/processed/obs_data_mintp_winter.csv") %>%
        dplyr::select(Station, 
                      scales_96_logged, threshold_tn_l_o_96_w_coast, 
                      scales_95_logged, threshold_tn_l_o_95_w_coast,
                      scales_97_logged, threshold_tn_l_o_97_w_coast,
                      Long.projected, Lat.projected, dist_sea_logged) %>%
        unique() 
      
      obs_grid = c()
      for(i in seq(length(spatial_covars$Station))){
        print(i)
        
        obs_grid = rbind(obs_grid, 
                         tibble(Station = spatial_covars$Station[i], 
                                year = (seq(1950, 2022))))
      }
      
      temporal_covars = read_csv("data/processed/obs_data_mintp_winter.csv") %>%
        dplyr::select(year, 
                      loess_temp_anom,
                      residuals
        ) %>% unique() %>%
        mutate(residuals = quantile(residuals, var_qnt)) %>%
        unique()
      
      if(thresh_qnt == 0.95){
        thresh_ex =read_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", marg_mod, "_bts_",bts,"_monthly_clim_covar_",var_qnt,".csv"),
                            col_names = c("bts", "Station", "year", "thresh_exceedance_95_coast")) %>%
          dplyr::select(-bts) %>% unique()
        
        obs_grid = obs_grid %>%
          left_join(temporal_covars) %>%
          left_join(thresh_ex)%>%
          left_join(spatial_covars)  %>%
          drop_na()
        
        obs_grid$scales_logged = obs_grid$scales_95_logged
        
        obs_grid$thresh_exceedance = obs_grid$thresh_exceedance_95_coast
        obs_grid$threshold = obs_grid$threshold_tn_l_o_95_w_coast
        
      }else if(thresh_qnt == 0.96){
        
        thresh_ex = read_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", marg_mod, "_bts_",bts,"_monthly_clim_covar_",var_qnt,".csv"),
                             col_names = c("bts", "Station", "year", "thresh_exceedance_96_coast")) %>%
          dplyr::select(-bts) %>% unique()
        
        obs_grid = obs_grid %>%
          left_join(temporal_covars) %>%
          left_join(thresh_ex)%>%
          left_join(spatial_covars)  %>%
          drop_na()
        
        obs_grid$scales_logged = obs_grid$scales_96_logged
        obs_grid$thresh_exceedance = obs_grid$thresh_exceedance_96_coast
        obs_grid$threshold = obs_grid$threshold_tn_l_o_96_w_coast
        
        
      }else if(thresh_qnt == 0.97){
        thresh_ex = read_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", marg_mod, "_bts_",bts,"_monthly_clim_covar_",var_qnt,".csv"),
                             col_names = c("bts", "Station", "year", "thresh_exceedance_97_coast")) %>%
          dplyr::select(-bts) %>% unique()
        
        obs_grid = obs_grid %>%
          left_join(temporal_covars) %>%
          left_join(thresh_ex)%>%
          left_join(spatial_covars)  %>%
          drop_na()
        
        obs_grid$scales_logged = obs_grid$scales_97_logged
        obs_grid$thresh_exceedance = obs_grid$thresh_exceedance_97_coast
        obs_grid$threshold = obs_grid$threshold_tn_l_o_97_w_coast
        
      }else{
        asfdgsafadafdgf()
      }
      
      
      fit = read.csv(paste0('output/bts_param_est/thresh_qnt_',thresh_qnt,'_corrected_',marg_mod,".csv"), header = FALSE, sep = ' ') %>%
       filter(V1 == bts) %>% .[-1] %>% as.numeric()
      
      
      
      # fit = readRDS(paste0('output/bts_param_est/thresh_qnt_',thresh_qnt,"_model_", marg_mod))
      pred = fit %>% predict_gpd(data = obs_grid, model = marg_mod)
      obs_grid$scale = pred$scale
      obs_grid$shape = pred$shape
      

      my_simulations_extremes = c()
      for(i in seq(1, 100)){
        my_simulations_extremes = c(my_simulations_extremes, readRDS(paste0("output/simulations/simulations_on_obs_grid/bts/thresh_qnt_",thresh_qnt,"_marg_mod_",marg_mod,"_bootstrap_",bts,"_run_",i)))
      }
      
      # ---- standardise simulations to have cost = 1
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
      

      m = length(my_simulations_standardised)
      L = 50
      sampled_rs = evd::rgpd(n=L, loc = 1, scale = 1, shape = 1)
      
      res_1950 = c()
      res_2020 = c()

      for(temp_i_want in seq(5, 20, by = 0.5)){
        
        print(temp_i_want)
        
        
        obs_grid$frechet = -1/log(1-obs_grid$thresh_exceedance*(1-evd::pgpd(temp_i_want - obs_grid$threshold,  scale = obs_grid$scale, shape =  obs_grid$shape[1] )) )
        

        T_2020 = obs_grid %>% filter(year == 2020) %>% pull(frechet)
        T_2022 = obs_grid %>% filter(year == 2022) %>% pull(frechet)
        T_1950 = obs_grid %>% filter(year == 1950) %>% pull(frechet)
        
        
        T_1950[T_1950 == -Inf] = Inf
        T_2020[T_2020 == -Inf] = Inf
        T_2022[T_2022 == -Inf] = Inf
        
        b_2020 = min(T_2020/max_at_each_site)
        b_2022 = min(T_2022/max_at_each_site)
        b_1950 = min(T_1950/max_at_each_site)
        
        
        m = length(my_simulations_standardised)
        
        above_1950 = 0
        above_2020 = 0
        above_2022 = 0
        
        # for each simulation
        for(j in seq(length(my_simulations_standardised))){
          # for each imp samp
          for(k in seq(L)){
            
            tmp = my_simulations_standardised[[j]]*sampled_rs[k]
            above_2022 = above_2022 + sum(sum(tmp*b_2022 > T_2022) > 0)
            above_2020 = above_2020 + sum(sum(tmp*b_2020 > T_2020) > 0)
            above_1950 = above_1950 + sum(sum(tmp*b_1950 > T_1950) > 0)
          }
        }
        
        tibble(temp = temp_i_want,
               bts,
               p_1950 = above_1950 / (length(my_simulations_standardised)*L*b_1950),
               p_2020 = above_2020 / (length(my_simulations_standardised)*L*b_2020),
               p_2022 = above_2022 / (length(my_simulations_standardised)*L*b_2022)) %>%
          write_csv(paste0("output/prob_extreme_temp_imp_samp_thresh_quant_",thresh_qnt,"_mar_mod_",marg_mod,"_var_qnt_",var_qnt,"_bootstrapped.csv"), append = T)
        # }
      }    
    }
  }
}

# 2 hours
job::job({get_bts_imp_samp("E2", seq(1, 50), 0.5, 0.96)})
job::job({get_bts_imp_samp("E2", seq(51, 100), 0.5, 0.96)})
job::job({get_bts_imp_samp("E2", seq(101, 150), 0.5, 0.96)})
job::job({get_bts_imp_samp("E2", seq(151, 200), 0.5, 0.96)})

job::job({get_bts_imp_samp("E2", seq(1, 50), 0.5, 0.95)})
job::job({get_bts_imp_samp("E2", seq(51, 100), 0.5, 0.95)})
job::job({get_bts_imp_samp("E2", seq(101, 150), 0.5, 0.95)})
job::job({get_bts_imp_samp("E2", seq(151, 200), 0.5, 0.95)})

job::job({get_bts_imp_samp("E2", seq(1, 50), 0.5, 0.97)})
job::job({get_bts_imp_samp("E2", seq(51, 100), 0.5, 0.97)})
job::job({get_bts_imp_samp("E2", seq(101, 150), 0.5, 0.97)})
job::job({get_bts_imp_samp("E2", seq(151, 200), 0.5, 0.97)})

job::job({get_bts_imp_samp("E2", seq(1, 50), 0.1, 0.96)})
job::job({get_bts_imp_samp("E2", seq(51, 100), 0.1, 0.96)})
job::job({get_bts_imp_samp("E2", seq(101, 150), 0.1, 0.96)})
job::job({get_bts_imp_samp("E2", seq(151, 200), 0.1, 0.96)})

job::job({get_bts_imp_samp("E2", seq(1, 50), 0.9, 0.96)})
job::job({get_bts_imp_samp("E2", seq(51, 100), 0.9, 0.96)})
job::job({get_bts_imp_samp("E2", seq(101, 150), 0.9, 0.96)})
job::job({get_bts_imp_samp("E2", seq(151, 200), 0.9, 0.96)})

job::job({get_bts_imp_samp("E1", seq(1, 50), 0.5, 0.96)})
job::job({get_bts_imp_samp("E1", seq(51, 100), 0.5, 0.96)})
job::job({get_bts_imp_samp("E1", seq(101, 150), 0.5, 0.96)})
job::job({get_bts_imp_samp("E1", seq(151, 200), 0.5, 0.96)})