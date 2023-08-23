# # Author: Dáire Healy
#
# Fits and saves quantile regression model and lambda estimates

# ---- 28 Feb. 2023

gc()
rm(list = ls())
library(tidyverse)
library(evgam)
setwd("~/chapter_6/")

quantiles_to_estimate_bulk = readRDS("output/quantiles_to_estimate_bulk_clim")
obs_data = readRDS(paste0("data/processed/obs_data_for_bulk_model_mintp_clim.csv"))




# ---- get covariates for prediction
temporal_covariates = obs_data %>%
  dplyr::select(year, month, loess_temp_anom, residuals) %>%
  unique() %>%
  arrange(year)


fit_clim_quants = T # bool, re-estimate clim quantiles?

if(fit_clim_quants){
  # file.remove(paste0("data/processed/qunatile_model_fit_pars_mintp.csv"))
  for(q in seq_along(quantiles_to_estimate_bulk)){
    
    #obs_data_for_quant_reg = obs_data
    obs_data_for_quant_reg = rlang::duplicate(obs_data, shallow = FALSE)
    zeta = obs_data_for_quant_reg$quantile[[1]][q] # quantile to estimate
    print(paste0("Fitting  tau = ", zeta))
    
    obs_data_for_quant_reg$value = obs_data_for_quant_reg$value %>% lapply(`[[`, q) %>% unlist
    qunatile_model_fit <- evgam(temp ~ value + loess_temp_anom + dist_sea_logged + residuals, obs_data_for_quant_reg,
                                family = "ald", ald.args = list(tau = zeta))
    
    # save parameter estimates
    tibble(tau = zeta,
           beta_0 = qunatile_model_fit$location$coefficients[1],
           beta_1 = qunatile_model_fit$location$coefficients[2],
           beta_2 = qunatile_model_fit$location$coefficients[3],
           beta_3 = qunatile_model_fit$location$coefficients[4],
           beta_4 = qunatile_model_fit$location$coefficients[5]) %>%
      write_csv(paste0("data/processed/qunatile_model_fit_pars_mintp_clim.csv"), append = T)
  }
}





# --- read in fitted quantile regression coefficients
quant_reg_model_pars = read_csv(paste0("data/processed/qunatile_model_fit_pars_mintp_clim.csv"),
                                col_names = c('tau', 'beta_0', 'beta_1', 'beta_2', 'beta_3', "beta_4"))





# # --- creates a tibble with each station and its quantule model
obs_smoothed_quantiles = obs_data %>%
  group_by(Station) %>%
  group_map(~{

    # .x = obs_data %>% filter(Station == 'shannon_airport')
    
    # --- get the climate qunatile estimates closest to current station
    clim_vals <<- obs_data %>%
      filter(Station == .x$Station[1]) %>%
      dplyr::select(quantile, value) %>%
      unique() %>%
      pull(value) %>%
      unlist()
    
    
    # --- spatial covariate
    dist_sea_vals = obs_data %>%
      filter(Station == .x$Station[1]) %>%
      dplyr::select(Station, dist_sea_logged) %>%
      unique() %>%
      pull(dist_sea_logged) %>%
      as.numeric()
    
    # predict quantile for each year and site
    quant_reg_pars = quant_reg_model_pars %>%
      arrange(tau)
    
    res = c()
    
    # interpolating through min and max temps
    min_and_max = read_csv("output/sites_w_min_max_obs") %>%
      filter(Station == .x$Station[1])
    
    
    res = rbind(res, tibble(quantile =  0,
                            year = temporal_covariates$year,
                            month = temporal_covariates$month,
                            loess_temp_anom = temporal_covariates$loess_temp_anom,
                            residuals = temporal_covariates$residuals,
                            quant_value = min_and_max$min_val),
                tibble(quantile =  1,
                       year = temporal_covariates$year,
                       month = temporal_covariates$month,
                       loess_temp_anom = temporal_covariates$loess_temp_anom,
                       residuals = temporal_covariates$residuals,
                       quant_value = min_and_max$max_val))
    
    
    for(q in seq_along(quantiles_to_estimate_bulk)){
      qpars = quant_reg_pars[q,]
      
      res = rbind(res,
                  tibble(quantile =  qpars$tau,
                         year = temporal_covariates$year,
                         month = temporal_covariates$month,
                         loess_temp_anom = temporal_covariates$loess_temp_anom,
                         residuals = temporal_covariates$residuals,
                         quant_value = qpars$beta_0 + 
                           qpars$beta_1*clim_vals[q] +
                           (qpars$beta_2)*(temporal_covariates$loess_temp_anom) +
                           (qpars$beta_3)*dist_sea_vals + 
                           (qpars$beta_4)*temporal_covariates$residuals))
    }
    
    print(paste0("Interpolating quantile estimates for ", .x$Station[1]))
    # interpolate quantiles over tau for each year
    
    # 
    # 
    #     # for each year at that site, interpolate over quantiles
    #     test = res %>%
    #       filter(year == 2010)
    # 
    #     tau_to_temp = splinefun(test$quantile, test$quant_value,  method = 'monoH.FC')
    #     temp_to_tau = splinefun(test$quant_value, test$quantile,  method = 'monoH.FC')
    # 
    #     tibble(x = seq(0,1, by = 0.01),
    #            y = tau_to_temp(seq(0,1, by = 0.01))) %>%
    #       ggplot()+
    #       geom_line(aes(x,y))+
    #       geom_point(data = test, aes(quantile, quant_value), col = 'red')
    # 
    # 
    #     res %>%
    #       ggplot() +
    #       geom_line(aes(quantile, quant_value, group = year, col = year))+
    #       facet_wrap(~month)+
    #       viridis::scale_color_viridis()+
    #       geom_hline(yintercept = .x$threshold_tn_l_o_96_w_coast)
    
    # res %>% saveRDS("output/shan_air")
    
    
    
    res %>%
      group_by(year, month) %>%
      group_map(~{
        tibble(year = .x$year[1],
               month = .x$month[1],
               tau_to_temp = list(splinefun(.x$quantile,.x$quant_value,  method = 'monoH.FC')),
               temp_to_tau = list(splinefun(.x$quant_value,.x$quantile,  method = 'monoH.FC')))
      }, .keep = T) %>%
      plyr::rbind.fill() %>%
      as_tibble() %>%
      mutate(Station = .x$Station[1])
    
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()



# save quantile models
obs_smoothed_quantiles %>% saveRDS(paste0("output/quant_models_mintp.csv"))
obs_smoothed_quantiles = readRDS(paste0("output/quant_models_mintp.csv"))




# Calculate lambda
lambda_thresh_ex = obs_data %>%
  group_by(Station) %>%
  group_map(~{
    
    print(.x$Station[1])
    
    
    # .x = obs_data %>% filter(Station == 'adare manor')
    thresh_exceedance_96 = obs_smoothed_quantiles%>%
      filter(Station == .x$Station[1]) %>%
      pull(temp_to_tau) %>%
      sapply(function(x) sapply(.x$threshold_tn_l_o_96[1], x))
    
    thresh_exceedance_96_coast = obs_smoothed_quantiles%>%
      filter(Station == .x$Station[1]) %>%
      pull(temp_to_tau) %>%
      sapply(function(x) sapply(.x$threshold_tn_l_o_96_w_coast[1], x))
    
    thresh_exceedance_97 = obs_smoothed_quantiles%>%
      filter(Station == .x$Station[1]) %>%
      pull(temp_to_tau) %>%
      sapply(function(x) sapply(.x$threshold_tn_l_o_97[1], x))
    
    thresh_exceedance_97_coast = obs_smoothed_quantiles%>%
      filter(Station == .x$Station[1]) %>%
      pull(temp_to_tau) %>%
      sapply(function(x) sapply(.x$threshold_tn_l_o_97_w_coast[1], x))
    
    thresh_exceedance_95 = obs_smoothed_quantiles%>%
      filter(Station == .x$Station[1]) %>%
      pull(temp_to_tau) %>%
      sapply(function(x) sapply(.x$threshold_tn_l_o_95[1], x))
    
    thresh_exceedance_95_coast = obs_smoothed_quantiles%>%
      filter(Station == .x$Station[1]) %>%
      pull(temp_to_tau) %>%
      sapply(function(x) sapply(.x$threshold_tn_l_o_95_w_coast[1], x))
    
    tibble(Station = .x$Station[1],
           month = temporal_covariates$month,
           year = temporal_covariates$year,
           thresh_exceedance_95 = 1-thresh_exceedance_95,
           thresh_exceedance_95_coast = 1-thresh_exceedance_95_coast,
           thresh_exceedance_96 = 1-thresh_exceedance_96,
           thresh_exceedance_96_coast = 1-thresh_exceedance_96_coast,
           thresh_exceedance_97 = 1-thresh_exceedance_97,
           thresh_exceedance_97_coast = 1-thresh_exceedance_97_coast)
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()

lambda_thresh_ex %>%
  write_csv(paste0("data/processed/thresh_exceedance_lambda_mintp.csv"))
# read_csv(paste0("data/processed/thresh_exceedance_lambda_num_quantiles_",30,"_mintp.csv"))




# 
