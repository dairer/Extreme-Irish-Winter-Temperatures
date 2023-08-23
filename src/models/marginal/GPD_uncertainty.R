# This script calculates uncertainties of scale and shape estimates based on bootrapped samples
rm(list = ls())
gc()
library(tidyverse)
models_to_run = c("E1", "E2")


read_csv("data/processed/obs_data_mintp_winter.csv") %>%
  dplyr::select(Station,
                Long, 
                Lat,
                year,
                dist_sea_logged,
                month,
                threshold_tn_l_o_95_w_coast,
                threshold_tn_l_o_96_w_coast,
                threshold_tn_l_o_97_w_coast,
                scales_95_logged,
                scales_96_logged,
                scales_97_logged,
                yearly_residuals,
                loess_temp_anom,
                residuals) %>%
  unique %>%
  saveRDS('output/bts_covars')



for(m in models_to_run){
  file.remove(paste0('output/bts_param_est/thresh_qnt_',0.95,'_uncorrected_',m,".csv"))
  file.remove(paste0('output/bts_param_est/thresh_qnt_',0.96,'_uncorrected_',m,".csv"))
  file.remove(paste0('output/bts_param_est/thresh_qnt_',0.97,'_uncorrected_',m,".csv"))
}


fit_uncorrected = function(models_to_run, bts_range, thresh_qnt){
  source("src/models/marginal/gpd_models.R")
  

  bts_covars = readRDS("output/bts_covars") %>%
    dplyr::select(-c(Long, Lat))

  
  if(thresh_qnt == 0.95){
    bts_covars$threshold =  bts_covars$threshold_tn_l_o_95_w_coast
    bts_covars$scales_logged =  bts_covars$scales_95_logged
  }else if(thresh_qnt == 0.96){
    bts_covars$threshold =  bts_covars$threshold_tn_l_o_96_w_coast
    bts_covars$scales_logged =  bts_covars$scales_96_logged
  }else if(thresh_qnt == 0.97){
    bts_covars$threshold =  bts_covars$threshold_tn_l_o_97_w_coast
    bts_covars$scales_logged =  bts_covars$scales_97_logged
  }else{
    lkzfhglkjsghls() # just fail
  }
  
  for(i in bts_range){
    print(paste0("Running bootstrap marginal model fit, ", i))
    for(m in models_to_run){
      this_bts = readRDS(paste0('output/bootstrapped_data_sets/thresh_qnt_',thresh_qnt,'_model_',m,"_bootstrap_", i)) %>%
        mutate(year = lubridate::year(date),
               month = lubridate::month(date)) %>%
        left_join(bts_covars)
      
      this_bts = this_bts %>%
        mutate(excess = temp - threshold) %>%
        filter(excess > 0)
      
      
     c(i,fit_gpd(this_bts, model = m)) %>% matrix() %>% t %>%
        write.table(paste0('output/bts_param_est/thresh_qnt_',thresh_qnt,'_uncorrected_',m,".csv") , append = T, col.names = F, row.names = F)
      
    }
  }
}

job::job({fit_uncorrected(models_to_run, seq(1, 200), 0.96)})
job::job({fit_uncorrected('E2', seq(1, 200), 0.95)})
job::job({fit_uncorrected('E2', seq(1, 200), 0.97)})



for(m in models_to_run){
  file.remove(paste0('output/bts_param_est/thresh_qnt_',0.95,'_corrected_',m,".csv"))
  file.remove(paste0('output/bts_param_est/thresh_qnt_',0.96,'_corrected_',m,".csv"))
  file.remove(paste0('output/bts_param_est/thresh_qnt_',0.97,'_corrected_',m,".csv"))
}



fit_corrected = function(models_to_run, bts_range, thresh_qnt){
  source("src/models/marginal/gpd_models.R")
  bts_covars = readRDS("output/bts_covars") %>%
    dplyr::select(-c(Long, Lat))
  
  if(thresh_qnt == 0.95){
    bts_covars$threshold =  bts_covars$threshold_tn_l_o_95_w_coast
    bts_covars$scales_logged =  bts_covars$scales_95_logged
  }else if(thresh_qnt == 0.96){
    bts_covars$threshold =  bts_covars$threshold_tn_l_o_96_w_coast
    bts_covars$scales_logged =  bts_covars$scales_96_logged
  }else if(thresh_qnt == 0.97){
    bts_covars$threshold =  bts_covars$threshold_tn_l_o_97_w_coast
    bts_covars$scales_logged =  bts_covars$scales_97_logged
  }else{
    lkzfhglkjsghls() # just fail
  }

  for(m in models_to_run){
    # --- fit on the actual data
    assign(paste0('pars_', m),readRDS(paste0('output/bts_param_est/thresh_qnt_',thresh_qnt,"_model_", m)))

    assign(paste0('uncorrected_', m),read.csv(paste0('output/bts_param_est/thresh_qnt_',thresh_qnt,'_uncorrected_',m,".csv"), header = FALSE, sep = ' ')[,-1])
    
    # shift bootstaped xi's by difference between mean and point estimate
    assign(paste0('correction_factor_', m),
           get(paste0('pars_', m))[length(get(paste0('pars_', m)))] - mean(get(paste0('uncorrected_', m))[,length(get(paste0('uncorrected_', m)))]))
  }
  
  
  
  for(i in bts_range){

    print(paste0("Running bootstrap, shape correction, ", i))
    for(m in models_to_run){

      this_bts = readRDS(paste0('output/bootstrapped_data_sets/thresh_qnt_',thresh_qnt,'_model_',m,"_bootstrap_", i)) %>%
        mutate(year = lubridate::year(date),
               month = lubridate::month(date)) %>%
        left_join(bts_covars)

      this_bts = this_bts %>%
        mutate(excess = temp - threshold) %>%
        filter(excess > 0)

      c(i, fit_gpd(this_bts,
              model = m,
              fix_shape = TRUE,
              initial_pars = signif(get(paste0('pars_', m))[-length(get(paste0('pars_', m)))], digits = 3),
              shape_est = as.numeric(get(paste0('uncorrected_', m))[i,][length(get(paste0('uncorrected_', m)))] + get(paste0('correction_factor_', m)))),
        as.numeric(get(paste0('uncorrected_', m))[i,][length(get(paste0('uncorrected_', m)))] + get(paste0('correction_factor_', m))))  %>%
        matrix() %>% t %>%
        write.table(paste0('output/bts_param_est/thresh_qnt_',thresh_qnt,'_corrected_',m,".csv") , append = T, col.names = F, row.names = F)
    }
  }
}

job::job({fit_corrected(models_to_run, seq(1, 200), 0.96)})
job::job({fit_corrected('E2', seq(1, 200), 0.95)})
job::job({fit_corrected('E2', seq(1, 200), 0.97)})