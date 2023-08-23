# find extreme temperatures on frechet scale
gc()
rm(list = ls())
library(tidyverse)


get_temp_frechet_scale_obs = function(marg_mod, var_qnt, thresh_qnt){
  source('src/models/marginal/gpd_models.R')
  
  
  
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
    thresh_ex = read_csv(paste0("data/processed/thresh_exceedance_lambda_mintp_",var_qnt,".csv")) %>%
      dplyr::select(Station, year, thresh_exceedance_95_coast) %>% unique()
    
    obs_grid = obs_grid %>%
      left_join(temporal_covars) %>%
      left_join(thresh_ex)%>%
      left_join(spatial_covars)  %>%
      drop_na()
    
    obs_grid$scales_logged = obs_grid$scales_95_logged
    
    obs_grid$thresh_exceedance = obs_grid$thresh_exceedance_95_coast
    obs_grid$threshold = obs_grid$threshold_tn_l_o_95_w_coast
        
  }else if(thresh_qnt == 0.96){
    
    thresh_ex = read_csv(paste0("data/processed/thresh_exceedance_lambda_mintp_",var_qnt,".csv")) %>%
      dplyr::select(Station, year, thresh_exceedance_96_coast) %>% unique()
    
    obs_grid = obs_grid %>%
      left_join(temporal_covars) %>%
      left_join(thresh_ex)%>%
      left_join(spatial_covars)  %>%
      drop_na()
    
    obs_grid$scales_logged = obs_grid$scales_96_logged
    obs_grid$thresh_exceedance = obs_grid$thresh_exceedance_96_coast
    obs_grid$threshold = obs_grid$threshold_tn_l_o_96_w_coast
    
    
  }else if(thresh_qnt == 0.97){
    thresh_ex = read_csv(paste0("data/processed/thresh_exceedance_lambda_mintp_",var_qnt,".csv")) %>%
      dplyr::select(Station, year, thresh_exceedance_97_coast) %>% unique()
    
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
  
  
  fit = readRDS(paste0('output/bts_param_est/thresh_qnt_',thresh_qnt,"_model_", marg_mod))
  pred = fit %>% predict_gpd(data = obs_grid, model = marg_mod)
  obs_grid$scale = pred$scale
  obs_grid$shape = pred$shape
  

  extreme_temps = seq(5, 20, by = 0.5)
  
  res = c()
  for(tmp in extreme_temps){
    res = cbind(res,(-1/log(1-obs_grid$thresh_exceedance*(1-evd::pgpd(q = (tmp - obs_grid$threshold),  scale = obs_grid$scale, shape =  obs_grid$shape[1])))))
  }
  
  res = as_tibble(res)
  names(res) = as.character(extreme_temps)

  res %>%
    mutate(year = obs_grid$year,
           Station = obs_grid$Station) %>%
    pivot_longer(-c(year,Station)) %>%
    rename(temp = name,
           frechet_value = value) %>%
    mutate(temp = as.numeric(temp)) %>%
    write_csv(paste0("output/obs_sites_extreme_temps_frechet_scale_thresh_qnt_",thresh_qnt,"_model_",marg_mod,"_qnt_",var_qnt,".csv"))
}

job::job({get_temp_frechet_scale_obs(marg_mod = 'E2', var_qnt = 0.5, thresh_qnt = 0.96)})
job::job({get_temp_frechet_scale_obs(marg_mod = 'E2', var_qnt = 0.5, thresh_qnt = 0.95)})
job::job({get_temp_frechet_scale_obs(marg_mod = 'E2', var_qnt = 0.5, thresh_qnt = 0.97)})

job::job({get_temp_frechet_scale_obs(marg_mod = 'E2', var_qnt = 0.1, thresh_qnt = 0.96)})
job::job({get_temp_frechet_scale_obs(marg_mod = 'E2', var_qnt = 0.9, thresh_qnt = 0.96)})
job::job({get_temp_frechet_scale_obs(marg_mod = 'E1', var_qnt = 0.5, thresh_qnt = 0.96)})




get_temp_frechet_scale_obs_bts = function(marg_mods, var_qnt, bts_range, thresh_qnt){
  
  for(marg_mod in marg_mods){
    for(bts_num in bts_range){
      

      print(bts_num)
      source('src/models/marginal/new_gpd_models.R')

      
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
                                year = c(1950, 2020, 2022)))

      }
      
      temporal_covars = read_csv("data/processed/obs_data_mintp_winter.csv") %>%
        dplyr::select(year, 
                      loess_temp_anom,
                      residuals) %>% unique() %>%
        mutate(residuals = quantile(residuals, var_qnt)) %>%
        unique()
      
      
      
      if(thresh_qnt == 0.95){
        
        thresh_ex = read_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", marg_mod, "_bts_",bts_num,"_monthly_clim_covar_",var_qnt,".csv"),
                             col_names = c('bts', 'Station', 'year' ,'thresh_exceedance_95_coast')) %>%
          dplyr::select(Station, year, thresh_exceedance_95_coast) %>% unique()
        
        obs_grid = obs_grid %>%
          left_join(temporal_covars) %>%
          left_join(thresh_ex)%>%
          left_join(spatial_covars)  %>%
          drop_na()
        
        obs_grid$scales_logged = obs_grid$scales_95_logged
        
        obs_grid$thresh_exceedance = obs_grid$thresh_exceedance_95_coast
        obs_grid$threshold = obs_grid$threshold_tn_l_o_95_w_coast
        
        
        

      }else if(thresh_qnt == 0.96){
        
        thresh_ex = read_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", marg_mod, "_bts_",bts_num,"_monthly_clim_covar_",var_qnt,".csv"),
                             col_names = c('bts', 'Station', 'year' ,'thresh_exceedance_96_coast')) %>%
          dplyr::select(Station, year, thresh_exceedance_96_coast) %>% unique()
        
        obs_grid = obs_grid %>%
          left_join(temporal_covars) %>%
          left_join(thresh_ex)%>%
          left_join(spatial_covars)  %>%
          drop_na()
        
        obs_grid$scales_logged = obs_grid$scales_96_logged
        obs_grid$thresh_exceedance = obs_grid$thresh_exceedance_96_coast
        obs_grid$threshold = obs_grid$threshold_tn_l_o_96_w_coast
        
        
      }else if(thresh_qnt == 0.97){
        
        thresh_ex = read_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", marg_mod, "_bts_",bts_num,"_monthly_clim_covar_",var_qnt,".csv"),
                             col_names = c('bts', 'Station', 'year' ,'thresh_exceedance_97_coast')) %>%
          dplyr::select(Station, year, thresh_exceedance_97_coast) %>% unique()
        
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
      
      
     
      assign(paste0('corrected_', marg_mod),read.csv(paste0('output/bts_param_est/thresh_qnt_',thresh_qnt,'_corrected_',marg_mod,".csv"), header = FALSE, sep = ' '))
      fit = get(paste0('corrected_', marg_mod)) 
      fit = fit %>% filter(V1 == bts_num) %>% as.numeric() %>% .[-1]
      
      pred = fit %>% predict_gpd(data = obs_grid, model = marg_mod)
      obs_grid$scale = pred$scale
      obs_grid$shape = pred$shape
      
      extreme_temps = seq(5, 18, by = 1)
      
      res = c()
      for(tmp in extreme_temps){
        res = cbind(res,(-1/log(1-obs_grid$thresh_exceedance*(1-evd::pgpd(q = (tmp - obs_grid$threshold),  scale = obs_grid$scale, shape =  obs_grid$shape[1])))))
      }
      
      res = as_tibble(res)
      names(res) = as.character(extreme_temps)
      
      res %>%
        mutate(bts_num = bts_num,
               year = obs_grid$year,
               Station = obs_grid$Station) %>%
        pivot_longer(-c(year,Station, bts_num)) %>%
        rename(temp = name,
               frechet_value = value) %>%
        mutate(temp = as.numeric(temp)) %>%
        
        write_csv(paste0("output/bootstrapped_obs_sites_extreme_temps_frechet_scale_thresh_qnt_",thresh_qnt,"marg_mod",marg_mod,"_qnt_",var_qnt,".csv"), append = T)
    }
  }
}

job::job({get_temp_frechet_scale_obs_bts(c('E2'), 0.5, seq(200), 0.96)})
job::job({get_temp_frechet_scale_obs_bts(c('E2'), 0.1, seq(200), 0.96)})
job::job({get_temp_frechet_scale_obs_bts(c('E2'), 0.9, seq(200), 0.96)})
# 
job::job({get_temp_frechet_scale_obs_bts(c('E2'), 0.5, seq(200), 0.95)})
job::job({get_temp_frechet_scale_obs_bts(c('E2'), 0.5, seq(200), 0.97)})
job::job({get_temp_frechet_scale_obs_bts(c('E1'), 0.5, seq(200), 0.96)})