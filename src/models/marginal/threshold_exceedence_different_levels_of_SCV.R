# this script estimtes thresold exceedance probability for different levels of short term climatic variability (SCV) and for differet thresholds for both
# point estimtes (on observed data) and on bostrapped estimates 


gc()
rm(list = ls())
library(tidyverse)
library(evgam)


calc_lambda_true = function(clim_osc_covar){
  
  # --- read in fitted quantile regression coefficients
  quant_reg_model_pars = read_csv(paste0("data/processed/qunatile_model_fit_pars_mintp.csv"),
                                  col_names = c('tau', 'beta_0', 'beta_1', 'beta_2', 'beta_3', "beta_4"))
  
  quantiles_to_estimate_bulk = readRDS("output/quantiles_to_estimate_bulk")
  obs_data = readRDS(paste0("data/processed/obs_data_for_bulk_model_mintp.csv"))
  
  # ---- get covariates for prediction
  temporal_covariates = obs_data %>%
    dplyr::select(year,month, residuals, loess_temp_anom) %>%
    unique() %>%
    arrange(year)
  
  temporal_covariates = temporal_covariates %>%
    mutate(residuals = quantile(residuals, clim_osc_covar)) %>%
    dplyr::select(-month) %>%
    unique()
  
  
  # # --- creates a tibble with each station and its quantule model
  obs_smoothed_quantiles = obs_data %>%
    group_by(Station)%>%
    group_map(~{
      
      # .x = obs_data %>% filter(Station == 'adare manor')
      
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
                           loess_temp_anom = temporal_covariates$loess_temp_anom,
                           residuals = temporal_covariates$residuals,
                           quant_value = qpars$beta_0 + 
                             qpars$beta_1*clim_vals[q] +
                             (qpars$beta_2)*(temporal_covariates$loess_temp_anom) +
                             (qpars$beta_3)*dist_sea_vals + 
                             (qpars$beta_4)*temporal_covariates$residuals))
      }
      
      print(paste0("Interpolating quantile estimates for ", .x$Station[1]))
      
      res %>%
        group_by(year) %>%
        group_map(~{
          tibble(year = .x$year[1],
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
  obs_smoothed_quantiles %>% saveRDS(paste0("output/quant_models_mintp_y_clim_covar_",clim_osc_covar,".csv"))
  
  obs_smoothed_quantiles = readRDS(paste0("output/quant_models_mintp_y_clim_covar_",clim_osc_covar,".csv"))
  
  
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
    write_csv(paste0("data/processed/thresh_exceedance_lambda_mintp_",clim_osc_covar,".csv"))
}

# --- 3 mins
job::job({calc_lambda_true(0.1)})
job::job({calc_lambda_true(0.5)})
job::job({calc_lambda_true(0.9)})







calc_lambda_bts = function(bts_range, marg_mod, thresh_qnt, clim_osc_covar){
  quantiles_to_estimate_bulk = readRDS("output/quantiles_to_estimate_bulk")
  obs_data = readRDS(paste0("data/processed/obs_data_for_bulk_model_mintp.csv"))
  
  # # # ---- get covariates for prediction
  temporal_covariates = obs_data %>%
    dplyr::select(year,month, residuals, loess_temp_anom) %>%
    unique() %>%
    arrange(year)
  
  temporal_covariates = temporal_covariates %>%
    mutate(residuals = quantile(residuals, clim_osc_covar)) %>%
    dplyr::select(-month) %>%
    unique()
  
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
                                loess_temp_anom = temporal_covariates$loess_temp_anom,
                                residuals = temporal_covariates$residuals,
                                quant_value = min_and_max$min_val),
                    tibble(quantile =  1,
                           year = temporal_covariates$year,
                           loess_temp_anom = temporal_covariates$loess_temp_anom,
                           residuals = temporal_covariates$residuals,
                           quant_value = min_and_max$max_val))
        
        
        for(q in seq_along(quantiles_to_estimate_bulk)){
          qpars = quant_reg_pars[q,]
          res =  rbind(res,
                       tibble(quantile =  qpars$tau,
                              year = temporal_covariates$year,
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
          group_by(year) %>% # and station???
          group_map(~{
            tibble(year = .x$year[1],
                   quant_spline = list(splinefun(.x$quant_value,.x$quantile,  method = 'monoH.FC')))
          }, .keep = T) %>%
          plyr::rbind.fill() %>%
          as_tibble() %>%
          mutate(Station = .x$Station[1])
        
      }, .keep = T) %>%
      plyr::rbind.fill() %>%
      as_tibble()

    if(thresh_qnt == 0.95){
      obs_data %>%
        group_by(Station) %>%
        dplyr::select(-c(threshold_tn_l_o_95,
                         threshold_tn_l_o_95_w_coast)) %>%
        left_join(readRDS(paste0("output/bootstrap_thresh/thresh_qnt_",thresh_qnt, "_model_",marg_mod,"_bootstrap_",bts_num))) %>%
        group_map(~{
          
          thresh_exceedance_95_coast = obs_smoothed_quantiles%>%
            filter(Station == .x$Station[1]) %>%
            pull(quant_spline) %>%
            sapply(function(x) sapply(.x$threshold_tn_l_o_95_w_coast[1], x))
          
          
          tibble(bts_num = bts_num,
                 Station = .x$Station[1],
                 year = temporal_covariates$year,
                 thresh_exceedance_95_coast = 1-thresh_exceedance_95_coast)  %>%
            write_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", marg_mod, "_bts_",bts_num,"_monthly_clim_covar_",clim_osc_covar,".csv"), append = T)
        
        }, .keep = T)    
      }else if(thresh_qnt == 0.96){
      
      obs_data %>%
        group_by(Station) %>%
        dplyr::select(-c(threshold_tn_l_o_96,
                         threshold_tn_l_o_96_w_coast)) %>%
        left_join(readRDS(paste0("output/bootstrap_thresh/thresh_qnt_",thresh_qnt, "_model_",marg_mod,"_bootstrap_",bts_num))) %>%
        group_map(~{
          
          thresh_exceedance_96_coast = obs_smoothed_quantiles%>%
            filter(Station == .x$Station[1]) %>%
            pull(quant_spline) %>%
            sapply(function(x) sapply(.x$threshold_tn_l_o_96_w_coast[1], x))
          
          
          tibble(bts_num = bts_num,
                 Station = .x$Station[1],
                 year = temporal_covariates$year,
                 thresh_exceedance_96_coast = 1-thresh_exceedance_96_coast)  %>%
            write_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", marg_mod, "_bts_",bts_num,"_monthly_clim_covar_",clim_osc_covar,".csv"), append = T)
        }, .keep = T)
      
    }else if(thresh_qnt == 0.97){
      
      obs_data %>%
        group_by(Station) %>%
        dplyr::select(-c(threshold_tn_l_o_97,
                         threshold_tn_l_o_97_w_coast)) %>%
        left_join(readRDS(paste0("output/bootstrap_thresh/thresh_qnt_",thresh_qnt, "_model_",marg_mod,"_bootstrap_",bts_num))) %>%
        group_map(~{
          
          thresh_exceedance_97_coast = obs_smoothed_quantiles%>%
            filter(Station == .x$Station[1]) %>%
            pull(quant_spline) %>%
            sapply(function(x) sapply(.x$threshold_tn_l_o_97_w_coast[1], x))
          
          
          tibble(bts_num = bts_num,
                 Station = .x$Station[1],
                 year = temporal_covariates$year,
                 thresh_exceedance_97_coast = 1-thresh_exceedance_97_coast)  %>%
            write_csv(paste0("output/bts_thresh_ex_lambda/bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", marg_mod, "_bts_",bts_num,"_monthly_clim_covar_",clim_osc_covar,".csv"), append = T)
        }, .keep = T)
      
    }else{
      asfdgsafadafdgf()
    }
  }
}



job::job({calc_lambda_bts(seq(1,50),  'E2', 0.96, 0.5)})
job::job({calc_lambda_bts(seq(51,100),  'E2', 0.96, 0.5)})
job::job({calc_lambda_bts(seq(101,150),  'E2', 0.96, 0.5)})
job::job({calc_lambda_bts(seq(151,200),  'E2', 0.96, 0.5)})
job::job({calc_lambda_bts(seq(1,100),  'E2', 0.96, 0.9)})
job::job({calc_lambda_bts(seq(101,200),  'E2', 0.96, 0.9)})
job::job({calc_lambda_bts(seq(39,100),  'E2', 0.96, 0.1)})
job::job({calc_lambda_bts(seq(101,200),  'E2', 0.96, 0.1)})
job::job({calc_lambda_bts(seq(1,100),  'E1', 0.96, 0.5)})
job::job({calc_lambda_bts(seq(101,200),  'E1', 0.96, 0.5)})




job::job({calc_lambda_bts(seq(1,50),  'E2', 0.95, 0.5)})
job::job({calc_lambda_bts(seq(51,100),  'E2', 0.95, 0.5)})
job::job({calc_lambda_bts(seq(101,150),  'E2', 0.95, 0.5)})
job::job({calc_lambda_bts(seq(151,200),  'E2', 0.95, 0.5)})
job::job({calc_lambda_bts(seq(1,100),  'E2', 0.97, 0.5)})
job::job({calc_lambda_bts(seq(101,200),  'E2', 0.97, 0.5)})




# --- plot resuts



marg_mod = 'E2'
thresh_qnt = 0.96
clim_osc_covar = 0.1
files = list.files("output/bts_thresh_ex_lambda/")
files = files[(!grepl("clim_grid_bootstrapped_threshold", files) )&grepl(paste0("_monthly_clim_covar_",clim_osc_covar), files) & grepl(marg_mod, files) & grepl(paste0("thresh_qnt_",thresh_qnt), files)]
files

all_bts_dat_1 = c()
for(f in files){
  print(f)
  all_bts_dat_1 = rbind(all_bts_dat_1,vroom::vroom(paste0("output/bts_thresh_ex_lambda/", f), col_names = c("bts", "Station",  "year", "thresh_exceedance_96_coast")))
}



clim_osc_covar = 0.5
files = list.files("output/bts_thresh_ex_lambda/")
files = files[(!grepl("clim_grid_bootstrapped_threshold", files) )&grepl(paste0("_monthly_clim_covar_",clim_osc_covar), files) & grepl(marg_mod, files) & grepl(paste0("thresh_qnt_",thresh_qnt), files)]
files

all_bts_dat_5 = c()
for(f in files){
  print(f)
  all_bts_dat_5 = rbind(all_bts_dat_5,vroom::vroom(paste0("output/bts_thresh_ex_lambda/", f), col_names = c("bts", "Station",  "year", "thresh_exceedance_96_coast")))
}



clim_osc_covar = 0.9
files = list.files("output/bts_thresh_ex_lambda/")
files = files[(!grepl("clim_grid_bootstrapped_threshold", files) )&grepl(paste0("_monthly_clim_covar_",clim_osc_covar), files) & grepl(marg_mod, files) & grepl(paste0("thresh_qnt_",thresh_qnt), files)]
files

all_bts_dat_9 = c()
for(f in files){
  print(f)
  all_bts_dat_9 = rbind(all_bts_dat_9,vroom::vroom(paste0("output/bts_thresh_ex_lambda/", f), col_names = c("bts", "Station",  "year", "thresh_exceedance_96_coast")))
}


all_bts_dat = rbind(all_bts_dat_1 %>% mutate(qnt = 0.1),
                    all_bts_dat_5 %>% mutate(qnt = 0.5),
                    all_bts_dat_9 %>% mutate(qnt = 0.9))
all_bts_dat %>% saveRDS("output/lambda_bts_data")


true_res = rbind(read_csv(paste0("data/processed/thresh_exceedance_lambda_mintp_",0.1,".csv")) %>%
                   group_by(year) %>%
                   summarise(thresh_exceedance_96_coast = mean(thresh_exceedance_96_coast)) %>% 
                   mutate(qnt = 0.1),
                 read_csv(paste0("data/processed/thresh_exceedance_lambda_mintp_",0.5,".csv")) %>%
                   group_by(year) %>%
                   summarise(thresh_exceedance_96_coast = mean(thresh_exceedance_96_coast))  %>% 
                   mutate(qnt = 0.5),
                 read_csv(paste0("data/processed/thresh_exceedance_lambda_mintp_",0.9,".csv")) %>%
                   group_by(year) %>%
                   summarise(thresh_exceedance_96_coast = mean(thresh_exceedance_96_coast))  %>% 
                   mutate(qnt = 0.9))

true_res %>%
  filter(year %in% c(2022, 1950))


#(0.0183 - 0.00688)/0.00688
#(0.0399 - 0.0153)/0.0153
#(0.107 - 0.0629)/0.0629
#(0.00688 - 0.0183)/0.0183
#(0.0153 - 0.0399)/0.0399
#(0.0629 - 0.107)/0.107




# --- now set equal to 0.1, 0.5 and 0.9?
plt = readRDS("output/lambda_bts_data") %>%
  group_by(bts, year, qnt) %>%
  summarise(thresh_exceedance_96_coast = mean(thresh_exceedance_96_coast))  %>%
  ungroup() %>%
  group_by(year, qnt) %>%
  summarise(upper = quantile(thresh_exceedance_96_coast, 0.98),
            lower = quantile(thresh_exceedance_96_coast, 0.02)) %>%
  ungroup() %>%
  ggplot()+
  geom_ribbon(aes(year, ymin = lower, ymax = upper), alpha = 0.3, fill = 'forestgreen')+
  geom_line(aes(year, lower), alpha = 0.5)+
  geom_line(aes(year, upper), alpha = 0.5)+
  geom_line(data = true_res, aes(year, thresh_exceedance_96_coast))+
  facet_wrap(~qnt, scale = 'free')+
  scale_x_continuous(breaks = c(1950, 1970, 1990, 2010, 2022))+
  theme_minimal(12)+
  labs(y = expression(lambda[t]),
       x = "Year")+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
        strip.text = element_blank())
ggsave("output/figs/thresh_ex.pdf", width =9, height = 2.5)
