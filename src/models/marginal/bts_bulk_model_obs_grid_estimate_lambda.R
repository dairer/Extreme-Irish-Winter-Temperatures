# this script fits bulk model to each bootstrapped data set and estimates threshold exceedence parameter for each bootstrap for a number of levels of SCV and thresolds

gc()
rm(list = ls())
library(tidyverse)
library(spatialsample)
library(evgam)

fit_quant_reg = T
calc_lam = T

fit_quant_regression = function(bts_range, marg_mod, thresh_qnt){

    # num_quantiles = 30
  quantiles_to_estimate_bulk = readRDS("output/quantiles_to_estimate_bulk")
  obs_data = readRDS(paste0("data/processed/obs_data_for_bulk_model_mintp.csv"))
  
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
        write_csv(paste0("output/bts_quant_reg_thresh_qnt_",thresh_qnt,"_model_",marg_mod,"_monthly.csv"), append = T)
    }
  }
}


# approx 1 hour each
# job::job({fit_quant_regression(seq(1,25),  'E2', 0.96)})
# job::job({fit_quant_regression(seq(26,50),  'E2', 0.96)})
# job::job({fit_quant_regression(seq(51,75),  'E2', 0.96)})
# job::job({fit_quant_regression(seq(76,100),  'E2', 0.96)})
# job::job({fit_quant_regression(seq(101,125),  'E2', 0.96)})
# job::job({fit_quant_regression(seq(126,150),  'E2', 0.96)})
# job::job({fit_quant_regression(seq(151,175),  'E2', 0.96)})
# job::job({fit_quant_regression(seq(176,200),  'E2', 0.96)})

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




# plot results
marg_mod = 'E2'
thresh_qnt = 0.96
true_pars = read_csv(paste0("data/processed/qunatile_model_fit_pars_mintp.csv"),
                     col_names = c("tau", "beta_0", "beta_1", "beta_2", "beta_3", "beta_4"))
res = read_csv(paste0("output/bts_quant_reg_thresh_qnt_",thresh_qnt,"_model_",marg_mod,"_monthly.csv"),
               col_names = c('bts', "tau", "b0", "b1", "b2","b3", "b4"))
res$bts %>% unique() %>% length()


plt= gridExtra::grid.arrange(res %>%
                          group_by(tau) %>%
                          summarise(upper = quantile(b1, 0.98),
                                    lower = quantile(b1, 0.02)) %>%
                          ggplot()+
                            geom_ribbon(aes(tau, ymin = lower, ymax = upper), alpha = 0.3, fill = 'forestgreen')+
                            geom_line(aes(tau, lower), alpha = 0.5)+
                            geom_line(aes(tau, upper), alpha = 0.5)+
                            geom_line(data = true_pars, aes(tau, beta_1))+
                          coord_cartesian(xlim = c(0.1, 0.9),  ylim = c(0.3, 0.9))+
                            # scale_x_continuous(breaks = c(0.1, 0.25, 0.5, 0.75, 0.9))+
                          theme_minimal(12)+
                          labs(y = expression(beta[1]),
                               x = expression(tau))+
                          theme(axis.title.y = element_text(angle = 0, vjust = 0.5)),
                        res %>%
                          group_by(tau) %>%
                          summarise(upper = quantile(b2, 0.98),
                                    lower = quantile(b2, 0.02)) %>%
                          ggplot()+
                          geom_ribbon(aes(tau, ymin = lower, ymax = upper), alpha = 0.3, fill = 'forestgreen')+
                          geom_line(aes(tau, lower), alpha = 0.5)+
                          geom_line(aes(tau, upper), alpha = 0.5)+
                          geom_line(data = true_pars , aes(tau, beta_2))+
                          coord_cartesian(xlim = c(0.1, 0.9),   ylim = c(-2.5, 0.2))+
                          # scale_x_continuous(breaks = c(0.1, 0.25, 0.5, 0.75, 0.9))+
                          theme_minimal(12)+
                          labs(y = expression(beta[2]),
                               x = expression(tau))+
                          theme(axis.title.y = element_text(angle = 0, vjust = 0.5)),
                        res %>%
                          group_by(tau) %>%
                          summarise(upper = quantile(b3, 0.98),
                                    lower = quantile(b3, 0.02)) %>%
                          ggplot()+
                          geom_ribbon(aes(tau, ymin = lower, ymax = upper), alpha = 0.3, fill = 'forestgreen')+
                          geom_line(aes(tau, lower), alpha = 0.5)+
                          geom_line(aes(tau, upper), alpha = 0.5)+
                          geom_line(data = true_pars, aes(tau, beta_3))+
                          coord_cartesian(xlim = c(0.1, 0.9),   ylim = c(0, 0.9))+
                          # scale_x_continuous(breaks = c(0.1, 0.25, 0.5, 0.75, 0.9))+
                          theme_minimal(12)+
                          labs(y = expression(beta[3]),
                               x = expression(tau))+
                          theme(axis.title.y = element_text(angle = 0, vjust = 0.5)),
                        res %>%
                          group_by(tau) %>%
                          summarise(upper = quantile(b4, 0.98),
                                    lower = quantile(b4, 0.02)) %>%
                          ggplot()+
                          geom_ribbon(aes(tau, ymin = lower, ymax = upper), alpha = 0.3, fill = 'forestgreen')+
                          geom_line(aes(tau, lower), alpha = 0.5)+
                          geom_line(aes(tau, upper), alpha = 0.5)+
                          geom_line(data = true_pars, aes(tau, beta_4))+
                          coord_cartesian(xlim = c(0.1, 0.9),   ylim = c(1, 2.5))+
                          # scale_x_continuous(breaks = c(0.1, 0.25, 0.5, 0.75, 0.9))+
                          theme_minimal(12)+
                          labs(y = expression(beta[4]),
                               x = expression(tau))+
                          theme(axis.title.y = element_text(angle = 0, vjust = 0.5)), nrow = 1)
                        

ggsave(plt, filename = "output/figs/bulk_params.pdf", height = 2.5, width = 9)




calc_lambda_bts = function(bts_range, marg_mod, thresh_qnt){
  quantiles_to_estimate_bulk = readRDS("output/quantiles_to_estimate_bulk")
  obs_data = readRDS(paste0("data/processed/obs_data_for_bulk_model_mintp.csv"))
  
  # # # ---- get covariates for prediction
  temporal_covariates = obs_data %>%
    dplyr::select(year,month, residuals, loess_temp_anom) %>%
    unique() %>%
    arrange(year)
  
  bootstrapped_models = read_csv(paste0("output/bts_quant_reg_thresh_qnt_",thresh_qnt,"_model_",marg_mod,"_monthly.csv"),
                                 col_names = c("bts", "tau", "beta_0", "beta_1", "beta_2", "beta_3", "beta_4"))
  
  for(bts_num in bts_range){

    
    print(paste0("Calculating lambda for bts num ", bts_num))
    this_est_of_qnt_mods = bootstrapped_models %>%
      filter(bts == bts_num)
    
    # # --- creates a tibble with each station and its quantule model
    obs_smoothed_quantiles = obs_data %>%
      group_by(Station) %>%
      group_map(~{
        .x = obs_data %>% filter(Station == "shannon_airport")
        
        
        
        dist_sea_logged_var = .x$dist_sea_logged[1]
        
        clim_vals = obs_data %>%
          filter(Station == .x$Station[1]) %>%
          .[1,] %>% 
          pull(value) %>%
          unlist()
        
        # predict quantile for each year and site
        quant_reg_pars = this_est_of_qnt_mods %>%
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
          res =  rbind(res,
                       tibble(quantile =  qpars$tau,
                              year = temporal_covariates$year,
                              month = temporal_covariates$month,
                              loess_temp_anom = temporal_covariates$loess_temp_anom,
                              residuals = temporal_covariates$residuals,
                              quant_value = qpars$beta_0 + 
                                qpars$beta_1*clim_vals[q] +
                                (qpars$beta_2)*(temporal_covariates$loess_temp_anom) +
                                (qpars$beta_3)*dist_sea_logged_var + 
                                (qpars$beta_4)*temporal_covariates$residuals))
        }
        
        # interpolate quantiles over tau for each year
        res %>%
          group_by(year, month) %>% # and station???
          group_map(~{
            tibble(year = .x$year[1],
                   quant_spline = list(splinefun(.x$quant_value,.x$quantile,  method = 'monoH.FC')))
            # tau_to_temp = list(splinefun(.x$quantile,.x$quant_value,  method = 'monoH.FC')),
            # temp_to_tau = list(splinefun(.x$quant_value,.x$quantile,  method = 'monoH.FC')))
          }, .keep = T) %>%
          plyr::rbind.fill() %>%
          as_tibble() %>%
          mutate(Station = .x$Station[1])
        
      }, .keep = T) %>%
      plyr::rbind.fill() %>%
      as_tibble()
    
    # look at threshold exceedance for this model
    obs_data %>%
      group_map(~{
        
        # print(.x$Station)
        # thresh_exceedance_95_coast = obs_smoothed_quantiles%>%
        #   filter(Station == .x$Station[1]) %>%
        #   pull(quant_spline) %>%
        #   sapply(function(x) sapply(.x$threshold_tn_l_o_95_w_coast[1], x))
        # 
        # thresh_exceedance_95 = obs_smoothed_quantiles%>%
        #   filter(Station == .x$Station[1]) %>%
        #   pull(quant_spline) %>%
        #   sapply(function(x) sapply(.x$threshold_tn_l_o_95[1], x))
        
        thresh_exceedance_96_coast = obs_smoothed_quantiles%>%
          filter(Station == .x$Station[1]) %>%
          pull(quant_spline) %>%
          sapply(function(x) sapply(.x$threshold_tn_l_o_96_w_coast[1], x))
        # thresh_exceedance_96 = obs_smoothed_quantiles%>%
        #   filter(Station == .x$Station[1]) %>%
        #   pull(quant_spline) %>%
        #   sapply(function(x) sapply(.x$threshold_tn_l_o_96[1], x))
        
        # thresh_exceedance_97_coast = obs_smoothed_quantiles%>%
        #   filter(Station == .x$Station[1]) %>%
        #   pull(quant_spline) %>%
        #   sapply(function(x) sapply(.x$threshold_tn_l_o_97_w_coast[1], x))
        # thresh_exceedance_97 = obs_smoothed_quantiles%>%
        #   filter(Station == .x$Station[1]) %>%
        #   pull(quant_spline) %>%
        #   sapply(function(x) sapply(.x$threshold_tn_l_o_97[1], x))
        
        tibble(bts_num = bts_num,
               Station = .x$Station[1],
               year = temporal_covariates$year,
               # thresh_exceedance_95 = 1-thresh_exceedance_95,
               # thresh_exceedance_96 = 1-thresh_exceedance_96,
               # thresh_exceedance_97 = 1-thresh_exceedance_97,
               # thresh_exceedance_95_coast = 1-thresh_exceedance_95_coast,
               thresh_exceedance_96_coast = 1-thresh_exceedance_96_coast#,
               # thresh_exceedance_97_coast = 1-thresh_exceedance_97_coast
               )  %>%
          write_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", marg_mod, "_bts_",bts_num,"_orig.csv"), append = T)
      }, .keep = T)
  }
}



job::job({calc_lambda_bts(seq(1,25),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(26,50),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(51,75),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(76,100),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(101,125),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(126,150),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(151,175),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(176,200),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(61,70),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(71,80),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(81,90),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(91,100),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(101,110),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(111,120),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(121,130),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(131,140),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(141,150),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(151,160),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(161,170),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(171,180),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(181,190),  'E2', 0.96)})
job::job({calc_lambda_bts(seq(191,200),  'E2', 0.96)})




# plot results

files = list.files("output/bts_thresh_ex_lambda/")
files = files[grepl("monthly.csv", files) & grepl(marg_mod, files) & grepl(paste0("thresh_qnt_",thresh_qnt), files)]
files


all_bts_dat = c()
for(f in files){
  print(f)
  all_bts_dat = rbind(all_bts_dat,read_csv(paste0("output/bts_thresh_ex_lambda/", f),
                                           col_names = c("bts", "Station", "month",  "year", "thresh_exceedance_96_coast")))
}


true_res = read_csv("data/processed/thresh_exceedance_lambda_mintp.csv") %>%
  group_by(year, month) %>%
  summarise(thresh_exceedance_96_coast = mean(thresh_exceedance_96_coast)) 




res = read_csv(paste0("output/bts_quant_reg_thresh_qnt_",thresh_qnt,"_model_",marg_mod,".csv"),
               col_names = c('bts', "tau", "b0", "b1", "b2","b3"))

all_bts_dat %>%
  group_by(bts, year, month) %>%
  summarise(thresh_exceedance_96_coast = mean(thresh_exceedance_96_coast))  %>%
  ungroup() %>%
  group_by(year, month) %>%
  summarise(upper = quantile(thresh_exceedance_96_coast, 0.98),
            lower = quantile(thresh_exceedance_96_coast, 0.02)) %>%
  ungroup() %>%
  ggplot()+
  geom_ribbon(aes(year, ymin = lower, ymax = upper), alpha = 0.3) +
  facet_wrap(~month)+
  geom_line(data = true_res, aes(year, thresh_exceedance_96_coast))+
  theme_minimal(12)+
  labs(y = expression(lambda[t]),
       x = "Year")+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))


plts=gridExtra::grid.arrange(
  res %>%
    group_by(tau) %>%
    summarise(upper = quantile(b1, 0.975),
              lower = quantile(b1, 0.025)) %>%
    ggplot()+
    geom_ribbon(aes(tau, ymin = lower, ymax = upper), alpha = 0.3)+
    geom_line(data = true_pars, aes(tau, beta_1))+
    coord_cartesian(xlim = c(0.1, 0.9))+
    theme_minimal(12)+
    labs(y = expression(beta[2]),
         x = expression(tau))+
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5)),
  res %>%
    group_by(tau) %>%
    summarise(upper = quantile(b2, 0.975),
              lower = quantile(b2, 0.025)) %>%
    ggplot()+
    geom_ribbon(aes(tau, ymin = lower, ymax = upper), alpha = 0.3)+
    geom_line(data = true_pars, aes(tau, beta_2))+
    coord_cartesian(xlim = c(0.1, 0.9))+
    theme_minimal(12)+
    labs(y = expression(beta[1]),
         x = expression(tau))+
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5)),
  res %>%
    group_by(tau) %>%
    summarise(upper = quantile(b3, 0.98),
              lower = quantile(b3, 0.02)) %>%
    ggplot()+
    geom_ribbon(aes(tau, ymin = lower, ymax = upper), alpha = 0.3)+
    geom_line(data = true_pars, aes(tau, beta_3))+
    theme_minimal(12)+
    labs(y = expression(beta[3]),
         x = expression(tau))+
    coord_cartesian(ylim=c(0.175, 0.65),
                    xlim = c(0.1, 0.9))+
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5)),
  all_bts_dat %>%
    group_by(bts, year) %>%
    summarise(thresh_exceedance_96_coast = mean(thresh_exceedance_96_coast))  %>%
    group_by(year) %>%
    summarise(upper = quantile(thresh_exceedance_96_coast, 0.975),
              lower = quantile(thresh_exceedance_96_coast, 0.025)) %>%
    ungroup() %>%
    ggplot()+
    geom_ribbon(aes(year, ymin = lower, ymax = upper), alpha = 0.3) +
    geom_line(data = true_res, aes(year, thresh_exceedance))+
    theme_minimal(12)+
    labs(y = expression(lambda[t]),
         x = "Year")+
    theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+
    xlim(c(1950, 2022)), nrow = 1)

ggsave(plts,filename = "output/figs/bulk_quantile_regression.pdf", height = 3, width = 7.5)
