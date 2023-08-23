gc()
rm(list = ls())
library(tidyverse)
library(spatialsample)
library(evgam)

fit_quant_reg = T
calc_lam = T

fit_quant_regression = function(bts_range, marg_mod, thresh_qnt){
  
  quantiles_to_estimate_bulk = readRDS("output/quantiles_to_estimate_bulk_clim")
  obs_data = readRDS(paste0("data/processed/obs_data_for_bulk_model_mintp_clim.csv"))
  
  # ---- get covariates for prediction
  temporal_covariates = obs_data %>%
    dplyr::select(year, month, loess_temp_anom, residuals) %>%
    unique() %>%
    arrange(year)
  
  spatial_covariates = obs_data %>%
    dplyr::select(Station, dist_sea_logged) %>%
    unique()
  
  x_subset1 <- grep(paste0("_",marg_mod,"_"), list.files('output/bootstrapped_data_sets/'), value = TRUE) 
  x_subset2 <- grep(paste0("thresh_qnt_",thresh_qnt), x_subset1, value = TRUE) 
  files_to_read = paste0("output/bootstrapped_data_sets/", x_subset2)
  
  files_to_keep = c()
  for(b in as.character(bts_range)){
    files_to_keep = c(files_to_keep, which(endsWith(files_to_read, paste0("_",b))))
  }
  
  files_to_read = files_to_read[files_to_keep]
  
  res = c()
  dat = c()
  
  for(file_name in files_to_read){
    bts_num = str_remove(file_name, paste0("output/bootstrapped_data_sets/thresh_qnt_",thresh_qnt, "_model_",marg_mod,"_bootstrap_"))
    
    dat <- readRDS(file_name)%>%
      mutate(year = lubridate::year(date),
             month = lubridate::month(date)) %>%
      left_join(temporal_covariates) %>%
      left_join(spatial_covariates)
    
    
    for(q in seq_along(quantiles_to_estimate_bulk)){
      
      zeta = obs_data$quantile[[1]][q] # quantile to estimate
      
      print(paste0("quantile ", zeta))
      
      # spatial covariate [q^tau_c(s)]
      obs_data$thisval = obs_data$value %>% lapply(`[[`, q) %>% unlist
      
      dat = dat %>%
        dplyr::select(date,Station, temp, year, month, loess_temp_anom, dist_sea_logged, residuals) %>%
        left_join(obs_data %>% dplyr::select(Station, value = thisval) %>% unique) %>% drop_na()
      
      # obs_data_for_quant_reg$value = obs_data_for_quant_reg$value %>% lapply(`[[`, q) %>% unlist
      qunatile_model_fit <- evgam(temp ~ value + loess_temp_anom + dist_sea_logged + residuals, data = dat,
                                  family = "ald", ald.args = list(tau = zeta))
      
      # save parameter estimates w CI
      tibble(bts = bts_num,
             tau = zeta,
             beta_0 = qunatile_model_fit$location$coefficients[1],
             beta_1 = qunatile_model_fit$location$coefficients[2],
             beta_2 = qunatile_model_fit$location$coefficients[3],
             beta_3 = qunatile_model_fit$location$coefficients[4],
             beta_4 = qunatile_model_fit$location$coefficients[5]) %>%
        write_csv(paste0("output/bts_quant_reg_thresh_qnt_",thresh_qnt,"_model_",marg_mod,"_monthly_clim.csv"), append = T)
    }
  }
}


# approx 1 hour to run each
# job::job({fit_quant_regression(seq(1,25),  'E2', 0.96)})
# job::job({fit_quant_regression(seq(26,50),  'E2', 0.96)})
# job::job({fit_quant_regression(seq(51,75),  'E2', 0.96)})
# job::job({fit_quant_regression(seq(76,100),  'E2', 0.96)})
# job::job({fit_quant_regression(seq(101,125),  'E2', 0.96)})
# job::job({fit_quant_regression(seq(126,150),  'E2', 0.96)})
# job::job({fit_quant_regression(seq(151,175),  'E2', 0.96)})
# job::job({fit_quant_regression(seq(176,200),  'E2', 0.96)})
# 
# job::job({fit_quant_regression(seq(1,25),  'E2', 0.95)})
# job::job({fit_quant_regression(seq(26,50),  'E2', 0.95)})
# job::job({fit_quant_regression(seq(51,75),  'E2', 0.95)})
# job::job({fit_quant_regression(seq(76,100),  'E2', 0.95)})
# job::job({fit_quant_regression(seq(101,125),  'E2', 0.95)})
# job::job({fit_quant_regression(seq(126,150),  'E2', 0.95)})
# job::job({fit_quant_regression(seq(151,175),  'E2', 0.95)})
# job::job({fit_quant_regression(seq(176,200),  'E2', 0.95)})
# 
# job::job({fit_quant_regression(seq(1,25),  'E2', 0.97)})
# job::job({fit_quant_regression(seq(26,50),  'E2', 0.97)})
# job::job({fit_quant_regression(seq(51,75),  'E2', 0.97)})
# job::job({fit_quant_regression(seq(76,100),  'E2', 0.97)})
# job::job({fit_quant_regression(seq(101,125),  'E2', 0.97)})
# job::job({fit_quant_regression(seq(126,150),  'E2', 0.97)})
# job::job({fit_quant_regression(seq(151,175),  'E2', 0.97)})
# job::job({fit_quant_regression(seq(176,200),  'E2', 0.97)})

# job::job({fit_quant_regression(seq(1,50),  'E1', 0.96)})
# job::job({fit_quant_regression(seq(51,100),  'E1', 0.96)})
# job::job({fit_quant_regression(seq(101,150),  'E1', 0.96)})
# job::job({fit_quant_regression(seq(151,200),  'E1', 0.96)})


