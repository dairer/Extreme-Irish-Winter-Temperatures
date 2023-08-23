
gc()
rm(list = ls())
library(tidyverse)
library(evgam)

calc_lambda_true_clim_grid = function(clim_osc_covar){
  
  quantiles_to_estimate_bulk = c(readRDS("output/quantiles_to_estimate_bulk"), readRDS("output/quantiles_to_estimate_bulk_clim"))

  quant_reg_model_pars = rbind(read_csv(paste0("data/processed/qunatile_model_fit_pars_mintp.csv"),
                                  col_names = c('tau', 'beta_0', 'beta_1', 'beta_2', 'beta_3', "beta_4")),
                               read_csv(paste0("data/processed/qunatile_model_fit_pars_mintp_clim.csv"),
                                        col_names = c('tau', 'beta_0', 'beta_1', 'beta_2', 'beta_3', "beta_4"))) %>%
    arrange(tau)
  
  
  # # ------------ get splines on clim scale
  clim_grid = read_csv("data/processed/winter_scales_clim_grid_mintp.csv") %>% dplyr::select(-c(scales_95, scales_96, scales_97))
  clim_dat_full = read_csv(paste0("data/processed/full_clim_data_mintp.csv"))


  dist_sea = read_csv("~/data/processed/sites_clim_sea_dist.csv") %>%
    dplyr::select(Long, Lat, dist_sea) %>%
    mutate(dist_sea_logged = log(1+dist_sea)) %>%
    mutate(dist_sea_logged = dist_sea_logged - readRDS("output/dist_sea_scaling")) %>%
    dplyr::select(-dist_sea)
  
  dist_sea$Lat = dist_sea$Lat %>% signif(4)
  dist_sea$Long = dist_sea$Long %>% signif(4)
  
  # ---- get covariates for prediction
  temporal_covariates = readRDS(paste0("data/processed/obs_data_for_bulk_model_mintp.csv")) %>%
    dplyr::select(year, month, loess_temp_anom, residuals) %>%
    unique() %>%
    arrange(year) %>%
    mutate(residuals = quantile(residuals, clim_osc_covar)) %>%
    dplyr::select(-month) %>%
    unique()
  
  # estimate empiracle quantules for climate data
  clim_quantiles = clim_dat_full %>%
    group_by(id) %>%
    group_map(~{
      res = c()
      for(q in quantiles_to_estimate_bulk){
        res = rbind(res, tibble(id = .x$id[1],
                                Long = .x$Long[1],
                                Lat = .x$Lat[1],
                                Long.projected = .x$Long.projected[1],
                                Lat.projected = .x$Lat.projected[1],
                                quantile = q,
                                value = as.numeric(quantile(.x$temp, q))))
      }
      res
    }, .keep = T)%>%
    plyr::rbind.fill() %>%
    as.tibble()
  
  
  clim_quantiles_subset = clim_quantiles %>%
    group_by(id) %>%
    group_map(~{
      tibble(id = .x$id[1], quantile = list(.x$quantile), value = list(.x$value))
    }, .keep = T) %>%
    plyr::rbind.fill() %>%
    as_tibble()
  
  
  clim_grid = clim_grid %>%
    left_join(clim_quantiles_subset) %>%
    left_join(dist_sea)
  clim_date_w_quantile_mod = c()
  
  

  # here
  for(i in seq(nrow(clim_grid))){
    print(clim_grid[i,])
    
    this_data_for_pred = tibble(year = c(1950, 2020, 2022)) %>%
      left_join(temporal_covariates) %>%
      mutate(Long = clim_grid[i,]$Long,
             dist_sea_logged = clim_grid[i,]$dist_sea_logged,
             Lat = clim_grid[i,]$Lat,
             Long.projected = clim_grid[i,]$Long.projected,
             Lat.projected = clim_grid[i,]$Lat.projected,
             id = clim_grid[i,]$id,
             quantile = clim_grid[i,]$quantile,
             value = clim_grid[i,]$value)
    
    # predict quantile for each year and site
    quant_reg_pars = quant_reg_model_pars %>%
      arrange(tau)
    
    clim_vals = clim_grid[i,]$value[[1]]
    res = c()
    for(q in seq_along(quantiles_to_estimate_bulk)){
      qpars = quant_reg_pars[q,]
      
      res = rbind(res,
                  tibble(quantile =  qpars$tau,
                         year = this_data_for_pred$year,
                         loess_temp_anom = this_data_for_pred$loess_temp_anom,
                         residuals = this_data_for_pred$residuals,
                         quant_value = qpars$beta_0 + 
                           qpars$beta_1*clim_vals[q] +
                           (qpars$beta_2)*(this_data_for_pred$loess_temp_anom) +
                           (qpars$beta_3)*this_data_for_pred$dist_sea_logged + 
                           (qpars$beta_4)*this_data_for_pred$residuals))
    }
    
    res = res %>%
      group_by(year) %>%
      group_map(~{
        tibble(year = .x$year[1],
               tau_to_temp = list(splinefun(.x$quantile,.x$quant_value,  method = 'monoH.FC')),
               temp_to_tau = list(splinefun(.x$quant_value,.x$quantile,  method = 'monoH.FC')))
      }, .keep = T) %>%
      plyr::rbind.fill() %>%
      as_tibble()%>%
      mutate(id = clim_grid[i,]$id[1])
    
    clim_date_w_quantile_mod = rbind(clim_date_w_quantile_mod,
                                     this_data_for_pred %>% left_join(res))
    
  }
  
  
  clim_date_w_quantile_mod = clim_date_w_quantile_mod %>%
    left_join(read_csv('data/processed/winter_scales_clim_grid_mintp.csv') %>%
                dplyr::select(Long, Lat,
                              threshold_95, scales_95,
                              threshold_96, scales_96,
                              threshold_97, scales_97))
  
  
  
  qunatile_model_fit_95 = readRDS("output/threshold_model/threshold_model_tn_l_95_coast")
  clim_date_w_quantile_mod$threshold_95_o = predict(qunatile_model_fit_95, newdata = clim_date_w_quantile_mod)$location
  rm("qunatile_model_fit_95")
  
  qunatile_model_fit_96 = readRDS("output/threshold_model/threshold_model_tn_l_96_coast")
  clim_date_w_quantile_mod$threshold_96_o = predict(qunatile_model_fit_96, newdata = clim_date_w_quantile_mod)$location
  rm("qunatile_model_fit_96")
  
  qunatile_model_fit_97 = readRDS("output/threshold_model/threshold_model_tn_l_97_coast")
  clim_date_w_quantile_mod$threshold_97_o = predict(qunatile_model_fit_97, newdata = clim_date_w_quantile_mod)$location
  rm("qunatile_model_fit_97")
  
  
  # Calculate lambda
  lambda_thresh_ex = clim_date_w_quantile_mod %>%
    group_by(id) %>%
    group_map(~{
      
      thresh_exceedance_95 = clim_date_w_quantile_mod %>%
        filter(id == .x$id[1]) %>%
        pull(temp_to_tau) %>%
        sapply(function(x) sapply(.x$threshold_95_o[1], x))
      
      
      thresh_exceedance_96 = clim_date_w_quantile_mod %>%
        filter(id == .x$id[1]) %>%
        pull(temp_to_tau) %>%
        sapply(function(x) sapply(.x$threshold_96_o[1], x))
      
      
      thresh_exceedance_97 = clim_date_w_quantile_mod %>%
        filter(id == .x$id[1]) %>%
        pull(temp_to_tau) %>%
        sapply(function(x) sapply(.x$threshold_97_o[1], x))
      
      
      tibble(id = .x$id[1],
             year = .x$year,
             thresh_exceedance_95 = 1- thresh_exceedance_95,
             thresh_exceedance_96 = 1- thresh_exceedance_96,
             thresh_exceedance_97 = 1- thresh_exceedance_97)
    }, .keep = T) %>%
    plyr::rbind.fill() %>%
    as_tibble()
  
  lambda_thresh_ex %>%
    write_csv(paste0("data/processed/CLIM_GRID_thresh_exceedance_lambda_mintp_",clim_osc_covar,".csv"))
}


# 7 mins
job::job({calc_lambda_true_clim_grid(0.1)})
job::job({calc_lambda_true_clim_grid(0.5)})
job::job({calc_lambda_true_clim_grid(0.9)})



calc_lambda_bootstraps_clim_grid = function(bts_range, marg_mod, clim_osc_covar, thresh_qnt){
  
  clim_grid = read_csv("data/processed/winter_scales_clim_grid_mintp.csv") %>% dplyr::select(-c(scales_95, scales_96, scales_97))
  clim_dat_full = read_csv(paste0("data/processed/full_clim_data_mintp.csv"))
  quantiles_to_estimate_bulk = c(readRDS("output/quantiles_to_estimate_bulk"), readRDS("output/quantiles_to_estimate_bulk_clim"))
  
  # ---- get covariates for prediction
  temporal_covariates = readRDS(paste0("data/processed/obs_data_for_bulk_model_mintp.csv")) %>%
    dplyr::select(year,month, residuals, loess_temp_anom) %>%
    unique() %>%
    arrange(year)
  
  temporal_covariates = temporal_covariates %>%
    mutate(residuals = quantile(residuals, clim_osc_covar)) %>%
    dplyr::select(-month) %>%
    filter(year %in% c(1950, 2020, 2022)) %>%
    unique()
  
  clim_quantiles = clim_dat_full %>%
    group_by(id) %>%
    group_map(~{
      res = c()
      for(q in quantiles_to_estimate_bulk){
        res = rbind(res, tibble(id = .x$id[1],
                                Long = .x$Long[1],
                                Lat = .x$Lat[1],
                                Long.projected = .x$Long.projected[1],
                                Lat.projected = .x$Lat.projected[1],
                                quantile = q,
                                value = as.numeric(quantile(.x$temp, q))))
      }
      res
    }, .keep = T) %>%
    plyr::rbind.fill() %>%
    as.tibble()
  
  
  clim_quantiles_subset = clim_quantiles %>%
    group_by(id) %>%
    group_map(~{
      tibble(id = .x$id[1], quantile = list(.x$quantile), value = list(.x$value))
    }, .keep = T) %>%
    plyr::rbind.fill() %>%
    as_tibble()

  bootstrapped_models = rbind(read_csv(paste0("output/bts_quant_reg_thresh_qnt_",thresh_qnt,"_model_",marg_mod,"_monthly.csv"),
                                       col_names = c("bts", "tau", "beta_0", "beta_1", "beta_2", "beta_3", "beta_4")),
                              read_csv(paste0("output/bts_quant_reg_thresh_qnt_",thresh_qnt,"_model_",marg_mod,"_monthly_clim.csv"),
                                       col_names = c("bts", "tau", "beta_0", "beta_1", "beta_2", "beta_3", "beta_4")))
  
  
  for(bts_num in bts_range){

    if(!file.exists(paste0("output/bts_thresh_ex_lambda/clim_grid_bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", marg_mod, "_bts_",bts_num,"_monthly_clim_covar_",clim_osc_covar,".csv"))){
        
      clim_date_w_quantile_mod = c()
      
      clim_grid = clim_grid %>%
        left_join(clim_quantiles_subset) %>%
        left_join(readRDS(paste0("output/bootstrapped_thresh_clim/thresh_qnt_",thresh_qnt, "_model_", marg_mod,"_bts_", bts_num))) 
      
      
      # here
      for(i in seq(nrow(clim_grid))){
        print(clim_grid[i,])
        
        this_data_for_pred = tibble(year = c(1950, 2020, 2022)) %>%
          left_join(temporal_covariates) %>%
          mutate(Long = clim_grid[i,]$Long,
                 dist_sea_logged = clim_grid[i,]$dist_sea_logged,
                 Lat = clim_grid[i,]$Lat,
                 Long.projected = clim_grid[i,]$Long.projected,
                 Lat.projected = clim_grid[i,]$Lat.projected,
                 id = clim_grid[i,]$id,
                 quantile = clim_grid[i,]$quantile,
                 value = clim_grid[i,]$value)
        
        # predict quantile for each year and site
        quant_reg_pars = bootstrapped_models %>%
          filter(bts == bts_num) %>%
          dplyr::select(-bts) %>%
          arrange(tau)
        
        clim_vals = clim_grid[i,]$value[[1]]
        
        dist_sea_vals = clim_grid[i,]$dist_sea_logged %>% as.numeric() %>% unique()
        
        
        res = c()
        for(q in seq_along(quantiles_to_estimate_bulk)){
          qpars = quant_reg_pars[q,]
          
          
          res = rbind(res,
                      tibble(quantile =  qpars$tau,
                             year = this_data_for_pred$year,
                             loess_temp_anom = this_data_for_pred$loess_temp_anom,
                             residuals = this_data_for_pred$residuals,
                             quant_value = qpars$beta_0 + 
                               qpars$beta_1*clim_vals[q] +
                               (qpars$beta_2)*(this_data_for_pred$loess_temp_anom) +
                               (qpars$beta_3)*dist_sea_vals + 
                               (qpars$beta_4)*this_data_for_pred$residuals))
        }
        
        res = res %>%
          group_by(year) %>%
          group_map(~{
            tibble(year = .x$year[1],
                   tau_to_temp = list(splinefun(.x$quantile,.x$quant_value,  method = 'monoH.FC')),
                   temp_to_tau = list(splinefun(.x$quant_value,.x$quantile,  method = 'monoH.FC')))
          }, .keep = T) %>%
          plyr::rbind.fill() %>%
          as_tibble()
        
        clim_date_w_quantile_mod = rbind(clim_date_w_quantile_mod,
                                         this_data_for_pred %>% left_join(res))
        
      }
      
      
      clim_date_w_quantile_mod = clim_date_w_quantile_mod %>%
        left_join(read_csv('data/processed/winter_scales_clim_grid_mintp.csv') %>%
                    dplyr::select(Long, Lat,
                                  threshold_95, scales_95,
                                  threshold_96, scales_96,
                                  threshold_97, scales_97))
      
      
      print(paste0("Calculating lambda for bts num ", bts_num))
      this_est_of_qnt_mods = bootstrapped_models %>%
        filter(bts == bts_num)
      
      
      if(thresh_qnt == 0.95){
        clim_grid %>%
          group_by(id) %>%
          group_map(~{
            
            thresh_exceedance_95_coast = clim_date_w_quantile_mod%>%
              filter(id == .x$id[1]) %>%
              pull(temp_to_tau) %>%
              sapply(function(x) sapply(.x$threshold_tn_l_o_95_w_coast[1], x))
            
            tibble(bts_num = bts_num,
                   id = .x$id[1],
                   year = temporal_covariates$year,
                   thresh_exceedance_95_coast = 1-thresh_exceedance_95_coast)  %>%
              write_csv(paste0("output/bts_thresh_ex_lambda/clim_grid_bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", marg_mod, "_bts_",bts_num,"_monthly_clim_covar_",clim_osc_covar,".csv"), append = T)
            
          }, .keep = T)   
        
      }else if(thresh_qnt == 0.96){
        
        clim_grid %>%
          group_by(id) %>%
          group_map(~{
            
            thresh_exceedance_96_coast = clim_date_w_quantile_mod%>%
              filter(id == .x$id[1]) %>%
              pull(temp_to_tau) %>%
              sapply(function(x) sapply(.x$threshold_tn_l_o_96_w_coast[1], x))
            
            tibble(bts_num = bts_num,
                   id = .x$id[1],
                   year = temporal_covariates$year,
                   thresh_exceedance_96_coast = 1-thresh_exceedance_96_coast)  %>%
              write_csv(paste0("output/bts_thresh_ex_lambda/clim_grid_bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", marg_mod, "_bts_",bts_num,"_monthly_clim_covar_",clim_osc_covar,".csv"), append = T)
            
          }, .keep = T)    
        
      }else if(thresh_qnt == 0.97){
        
        clim_grid %>%
          group_by(id) %>%
          group_map(~{
            
            thresh_exceedance_97_coast = clim_date_w_quantile_mod%>%
              filter(id == .x$id[1]) %>%
              pull(temp_to_tau) %>%
              sapply(function(x) sapply(.x$threshold_tn_l_o_97_w_coast[1], x))
            
            tibble(bts_num = bts_num,
                   id = .x$id[1],
                   year = temporal_covariates$year,
                   thresh_exceedance_97_coast = 1-thresh_exceedance_97_coast)  %>%
              write_csv(paste0("output/bts_thresh_ex_lambda/clim_grid_bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", marg_mod, "_bts_",bts_num,"_monthly_clim_covar_",clim_osc_covar,".csv"), append = T)
            
          }, .keep = T)    
        
      }else{
        asfdgsafadafdgf()
      }
    }
  }
}

# 1.5 hours each
# job::job({calc_lambda_bootstraps_clim_grid(seq(1, 50), "m22", 0.5, 0.96)})
# job::job({calc_lambda_bootstraps_clim_grid(seq(51, 100), "m22", 0.5, 0.96)})
# job::job({calc_lambda_bootstraps_clim_grid(seq(101, 150), "m22", 0.5, 0.96)})
# job::job({calc_lambda_bootstraps_clim_grid(seq(151, 200), "m22", 0.5, 0.96)})


# job::job({calc_lambda_bootstraps_clim_grid(seq(1, 50), "m22", 0.1, 0.96)})
# job::job({calc_lambda_bootstraps_clim_grid(seq(51, 100), "m22", 0.1, 0.96)})
# job::job({calc_lambda_bootstraps_clim_grid(seq(101, 150), "m22", 0.1, 0.96)})
# job::job({calc_lambda_bootstraps_clim_grid(seq(151, 200), "m22", 0.1, 0.96)})

# job::job({calc_lambda_bootstraps_clim_grid(seq(1, 50),    "m22", 0.9, 0.96)})
# job::job({calc_lambda_bootstraps_clim_grid(seq(51, 100),    "m22", 0.9, 0.96)})
# job::job({calc_lambda_bootstraps_clim_grid(seq(101, 150),    "m22", 0.9, 0.96)})
# job::job({calc_lambda_bootstraps_clim_grid(seq(151, 200),    "m22", 0.9, 0.96)})


# job::job({calc_lambda_bootstraps_clim_grid(seq(200), "m22", 0.5, 0.95)})
# job::job({calc_lambda_bootstraps_clim_grid(seq(200), "m22", 0.5, 0.97)})
# job::job({calc_lambda_bootstraps_clim_grid(seq(200), "m21", 0.5, 0.96)})
