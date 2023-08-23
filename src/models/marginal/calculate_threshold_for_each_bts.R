# this script calcualtes the threshold for each bootstrapped sample, on observational and climate grid

rm(list=ls())
library(tidyverse)

clim_data = read_csv('data/processed/winter_scales_clim_grid_mintp.csv')

read_csv("data/processed/obs_data_mintp_winter.csv") %>%
  dplyr::select(Station, threshold_96,threshold_97, threshold_95, dist_sea_logged) %>% 
  unique() %>%
  saveRDS("output/thresh_mod_covars")



threshold_bts = function(bts_range, marg_mod, thresh_qnt){
  x_subset1 <- grep(paste0("_",marg_mod,"_"), list.files('output/bootstrapped_data_sets/'), value = TRUE) 
  x_subset2 <- grep(paste0("thresh_qnt_",thresh_qnt), x_subset1, value = TRUE) 
  files_to_read = paste0("output/bootstrapped_data_sets/", x_subset2)
  
  files_to_keep = c()
  for(b in as.character(bts_range)){
    files_to_keep = c(files_to_keep, which(endsWith(files_to_read, paste0("_",b))))
  }
  
  files_to_read = files_to_read[files_to_keep]

  
  for(file_name in files_to_read){

    bts_num = str_remove(file_name, paste0("output/bootstrapped_data_sets/thresh_qnt_",thresh_qnt, "_model_",marg_mod,"_bootstrap_"))

    dat <- readRDS(file_name)%>%
      mutate(year = lubridate::year(date),
             month = lubridate::month(date)) %>%
      left_join(readRDS("output/thresh_mod_covars"))
    
    
    if(thresh_qnt == 0.95){

      qunatile_model_95_w_coast <-  temp ~ threshold_95 + dist_sea_logged
      qunatile_model_fit_95_w_coast <- evgam::evgam(qunatile_model_95_w_coast, dat, family = "ald", ald.args = list(tau = 0.95))
      dat$threshold_tn_l_o_95_w_coast = qunatile_model_fit_95_w_coast$location$fitted
      
      dat %>%
        dplyr::select(Station, threshold_tn_l_o_95_w_coast) %>%
        unique() %>%
        saveRDS(paste0("output/bootstrap_thresh/thresh_qnt_",thresh_qnt, "_model_",marg_mod,"_bootstrap_",bts_num))
    
      }else if(thresh_qnt == 0.96){
      
      qunatile_model_96_w_coast <-  temp ~ threshold_96 + dist_sea_logged
      qunatile_model_fit_96_w_coast <- evgam::evgam(qunatile_model_96_w_coast, dat, family = "ald", ald.args = list(tau = 0.96))
      dat$threshold_tn_l_o_96_w_coast = qunatile_model_fit_96_w_coast$location$fitted
      
      dat %>%
        dplyr::select(Station, threshold_tn_l_o_96_w_coast) %>%
        unique() %>%
        saveRDS(paste0("output/bootstrap_thresh/thresh_qnt_",thresh_qnt, "_model_",marg_mod,"_bootstrap_",bts_num))
      
    }else if(thresh_qnt == 0.97){

      qunatile_model_97_w_coast <-  temp ~ threshold_97 + dist_sea_logged
      qunatile_model_fit_97_w_coast <- evgam::evgam(qunatile_model_97_w_coast, dat, family = "ald", ald.args = list(tau = 0.97))
      dat$threshold_tn_l_o_97_w_coast = qunatile_model_fit_97_w_coast$location$fitted
      
      dat %>%
        dplyr::select(Station, threshold_tn_l_o_97_w_coast) %>%
        unique() %>%
        saveRDS(paste0("output/bootstrap_thresh/thresh_qnt_",thresh_qnt, "_model_",marg_mod,"_bootstrap_",bts_num))
    
      }else{
      asfdgsafadafdgf() # cause an error, used for de bugging
    }
  }
}

job::job({threshold_bts(seq(1, 200), "E2", 0.95)})
job::job({threshold_bts(seq(1, 200), "E2", 0.96)})
job::job({threshold_bts(seq(1, 200), "E2", 0.97)})
job::job({threshold_bts(seq(1, 200), "E2", 0.96)})



threshold_bts_clim = function(bts_range, marg_mod, thresh_qnt){
  clim_grid = read_csv("data/processed/winter_scales_clim_grid_mintp.csv") %>% dplyr::select(-c(scales_95, scales_96, scales_97))
  
  dist_sea = read_csv("~/Inference for extreme spatial temperature events in a changing climate with application to Ireland/data/processed/sites_clim_sea_dist.csv") %>%
    dplyr::select(Long, Lat, dist_sea) %>%
    mutate(dist_sea_logged = log(1+dist_sea)) %>%
    mutate(dist_sea_logged = dist_sea_logged - readRDS("output/dist_sea_scaling")) %>%
    dplyr::select(-dist_sea)
  
  dist_sea$Lat = dist_sea$Lat %>% signif(4)
  dist_sea$Long = dist_sea$Long %>% signif(4)
  
  
  clim_grid = clim_grid %>%
    left_join(dist_sea)

  x_subset1 <- grep(paste0("_",marg_mod,"_"), list.files('output/bootstrapped_data_sets/'), value = TRUE) 
  x_subset2 <- grep(paste0("thresh_qnt_",thresh_qnt), x_subset1, value = TRUE) 
  files_to_read = paste0("output/bootstrapped_data_sets/", x_subset2)
  
  files_to_keep = c()
  for(b in as.character(bts_range)){
    files_to_keep = c(files_to_keep, which(endsWith(files_to_read, paste0("_",b))))
  }
  
  files_to_read = files_to_read[files_to_keep]

  for(file_name in files_to_read){
    
    bts_num = str_remove(file_name, paste0("output/bootstrapped_data_sets/thresh_qnt_",thresh_qnt, "_model_",marg_mod,"_bootstrap_"))
    
    dat <- readRDS(file_name)%>%
      mutate(year = lubridate::year(date),
             month = lubridate::month(date)) %>%
      left_join(readRDS("output/thresh_mod_covars"))
    
    if(thresh_qnt == 0.95){
      qunatile_model_95_w_coast <-  temp ~ threshold_95 + dist_sea_logged
      qunatile_model_fit_95_w_coast <- evgam::evgam(qunatile_model_95_w_coast, dat, family = "ald", ald.args = list(tau = 0.95))
      clim_grid$threshold_tn_l_o_95_w_coast = predict(qunatile_model_fit_95_w_coast, newdata = clim_grid)$location 
      clim_grid %>%
        dplyr::select(id, dist_sea_logged, threshold_tn_l_o_95_w_coast) %>%
        saveRDS(paste0("output/bootstrapped_thresh_clim/thresh_qnt_",thresh_qnt, "_model_", marg_mod,"_bts_", bts_num))
      
    }else if(thresh_qnt == 0.96){
      
      qunatile_model_96_w_coast <-  temp ~ threshold_96 + dist_sea_logged
      qunatile_model_fit_96_w_coast <- evgam::evgam(qunatile_model_96_w_coast, dat, family = "ald", ald.args = list(tau = 0.96))
      clim_grid$threshold_tn_l_o_96_w_coast = predict(qunatile_model_fit_96_w_coast, newdata = clim_grid)$location 
      clim_grid %>%
        dplyr::select(id, dist_sea_logged, threshold_tn_l_o_96_w_coast) %>%
        saveRDS(paste0("output/bootstrapped_thresh_clim/thresh_qnt_",thresh_qnt, "_model_", marg_mod,"_bts_", bts_num))
      
      
    }else if(thresh_qnt == 0.97){
      
      qunatile_model_97_w_coast <-  temp ~ threshold_97 + dist_sea_logged
      qunatile_model_fit_97_w_coast <- evgam::evgam(qunatile_model_97_w_coast, dat, family = "ald", ald.args = list(tau = 0.97))
      clim_grid$threshold_tn_l_o_97_w_coast = predict(qunatile_model_fit_97_w_coast, newdata = clim_grid)$location 
      clim_grid %>%
        dplyr::select(id, dist_sea_logged, threshold_tn_l_o_97_w_coast) %>%
        saveRDS(paste0("output/bootstrapped_thresh_clim/thresh_qnt_",thresh_qnt, "_model_", marg_mod,"_bts_", bts_num))
      
    }else{
      asfdgsafadafdgf()# cause an error and hault
    }

   }
}


job::job({threshold_bts_clim(seq(1, 200), "E2", 0.95)})
job::job({threshold_bts_clim(seq(1, 200), "E2", 0.96)})
job::job({threshold_bts_clim(seq(1, 200), "E2", 0.97)})
job::job({threshold_bts_clim(seq(1, 200), "E2", 0.96)})
