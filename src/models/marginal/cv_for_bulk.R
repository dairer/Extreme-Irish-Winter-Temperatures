# Cross validation for bulk (quantile regression) models
gc()
rm(list = ls())

run_bulk_cv = function(mod, seq_id, cv_type){
  library(tidyverse)
  library(spatialsample)
  library(evgam)
  
  num_quantiles = 15
  quantiles_to_estimate_bulk = seq(0.001,0.999,length.out = num_quantiles)
  
  obs_data = readRDS(paste0("data/processed/obs_data_for_bulk_model_num_quantiles_",num_quantiles,"_mintp.csv"))
  # obs_data$dist_sea_logged = log(1+obs_data$dist_sea_logged)

  yrs = obs_data$year %>% unique()
  obs_data = obs_data %>%
    group_by(year, Station) %>%
    mutate(emp_quantile = ecdf(temp)(temp)) %>%
    ungroup()

  fit_qunat_mod = function(train_set, test_set, mod){
    # ------- define quqntile reg models
    if(mod == 'ma') qunatile_model <- temp ~ 1 # m0 constatn
    if(mod == 'mb') qunatile_model <- temp ~ value # m1 spatial
    if(mod == 'mc') qunatile_model <- temp ~ value + loess_temp_anom # m2 spatio-temporal
    if(mod == 'md') qunatile_model <- temp ~ value + loess_temp_anom + dist_sea_logged
    if(mod == 'me') qunatile_model <- temp ~ loess_temp_anom 
    if(mod == 'mf') qunatile_model <- temp ~ loess_temp_anom  + dist_sea_logged
    if(mod == 'mg') qunatile_model <- temp ~ value + loess_temp_anom + dist_sea_logged + residuals
                                
    
    # ------ fit quantile model to observational data for each tau
    obs_q_mod = c()
    for(q in seq_along(quantiles_to_estimate_bulk)){
      obs_data_for_quant_reg = train_set
      zeta = obs_data_for_quant_reg$quantile[[1]][q] # quantile to estimate
      print(paste0("Fitting ", mod, " to tau = ", zeta))
      
      # spatial covariate [q^tau_c(s)]
      obs_data_for_quant_reg$value = obs_data_for_quant_reg$value %>% lapply(`[[`, q) %>% unlist
      qunatile_model_fit <- evgam(qunatile_model, obs_data_for_quant_reg,
                                  family = "ald", ald.args = list(tau = zeta))
      # save quantile models 
      obs_q_mod = rbind(obs_q_mod, tibble(q, quant_mod = list(qunatile_model_fit)))
    }
    
    if(mod == 'ma'){
      # --- constant in space and time
      res = c()
      
      for(q in seq_along(quantiles_to_estimate_bulk)){
        res = rbind(res,tibble(quantile =  quantiles_to_estimate_bulk[q],
                               value = unique(obs_q_mod[q,]$quant_mod %>% lapply(predict) %>% .[[1]] %>% .$location)))
      }
      
      quant_spline = splinefun(res$quantile,res$value,  method = 'monoH.FC')
      pred = quant_spline(test_set$emp_quantile)
      sqrt(mean((test_set$temp-pred)^2))
      
    }else if(mod == "mb"){
      
      # --- creates a tibble with each station and its quantule model
      obs_smoothed_quantiles = train_set %>%
        group_by(Station) %>%
        group_map(~{
          
          # get the climate qunatile estimates closest to current station
          clim_vals = train_set %>%
            filter(Station == .x$Station[1]) %>%
            dplyr::select(quantile, value) %>%
            unique() %>%
            pull(value) %>%
            unlist()
          
          # for each tau, for each location, get quantile estimate
          res = c()
          for(q in seq_along(quantiles_to_estimate_bulk)){
            res = rbind(res,tibble(quantile =  quantiles_to_estimate_bulk[q],
                                   value = unique(obs_q_mod[q,]$quant_mod %>% lapply(predict,tibble(value = clim_vals[q]) ) %>% .[[1]] %>% .$location)))
          }
          tibble(Station = .x$Station[1], quant_spline = list(splinefun(res$quantile, res$value,  method = 'monoH.FC')))
          
        }, .keep = T) %>%
        plyr::rbind.fill() %>%
        as_tibble()
      
      # ------- calculate RMSE
      squared_error = test_set %>%
        group_by(Station) %>%
        group_map(~{
          
          this_quant_spline = obs_smoothed_quantiles %>%
            filter(Station == .x$Station[1]) %>%
            pull(quant_spline) %>%
            .[[1]]
          
          pred = this_quant_spline(.x$emp_quantile)
          (.x$temp-pred)^2
          
        }, .keep = T) %>%
        unlist
      sqrt(mean(squared_error))
      
    }else{
      
      # --- creates a tibble with each station and its quantule model 
      obs_smoothed_quantiles = train_set %>%
        group_by(Station) %>%
        group_map(~{
          
          # --- get the climate qunatile estimates closest to current station
          clim_vals = train_set %>%
            filter(Station == .x$Station[1]) %>%
            dplyr::select(quantile, value) %>%
            unique() %>%
            pull(value) %>%
            unlist()
          
          res = c()
          
        if(mod == "mc"){
            # --- get temporal covariate
            loess_temp_anom_vals = obs_data %>%
              dplyr::select(year, loess_temp_anom) %>%
              unique() %>%
              arrange(year) %>%
              pull(loess_temp_anom) %>%
              as.numeric()
            
            # for each tau, for each location, get quantile estimate
            for(q in seq_along(quantiles_to_estimate_bulk)){
              res = rbind(res,tibble(year = yrs,
                                     quantile =  quantiles_to_estimate_bulk[q],
                                     value = obs_q_mod[q,]$quant_mod %>% lapply(predict, tibble(value = clim_vals[q], loess_temp_anom = loess_temp_anom_vals)) %>% .[[1]] %>% .$location))
            }
          }else if(mod == "md"){
            loess_temp_anom_vals = obs_data %>%
              dplyr::select(year, loess_temp_anom) %>%
              unique() %>%
              arrange(year) %>%
              pull(loess_temp_anom) %>%
              as.numeric()
          
            # --- spatial covariate
            dist_sea_logged_vals = obs_data %>%
              filter(Station == .x$Station[1]) %>%
              dplyr::select(Station, dist_sea_logged) %>%
              unique() %>%
              pull(dist_sea_logged) %>%
              as.numeric()
  
            # for each tau, for each location, get quantile estimate
            for(q in seq_along(quantiles_to_estimate_bulk)){
              res = rbind(res,tibble(year = yrs,
                                     quantile =  quantiles_to_estimate_bulk[q],
                                     value = obs_q_mod[q,]$quant_mod %>% lapply(predict, tibble(value = clim_vals[q], loess_temp_anom = loess_temp_anom_vals, dist_sea_logged = dist_sea_logged_vals)) %>% .[[1]] %>% .$location))
            }
          }else if(mod == "me"){
            loess_temp_anom_vals = obs_data %>%
              dplyr::select(year, loess_temp_anom) %>%
              unique() %>%
              arrange(year) %>%
              pull(loess_temp_anom) %>%
              as.numeric()
            
            # for each tau, for each location, get quantile estimate
            for(q in seq_along(quantiles_to_estimate_bulk)){
              res = rbind(res,tibble(year = yrs,
                                     quantile =  quantiles_to_estimate_bulk[q],
                                     value = obs_q_mod[q,]$quant_mod %>% lapply(predict, tibble(loess_temp_anom = loess_temp_anom_vals)) %>% .[[1]] %>% .$location))
            }
          }else if(mod == "mf"){
            loess_temp_anom_vals = obs_data %>%
              dplyr::select(year, loess_temp_anom) %>%
              unique() %>%
              arrange(year) %>%
              pull(loess_temp_anom) %>%
              as.numeric()
            
            # --- spatial covariate
            dist_sea_logged_vals = obs_data %>%
              filter(Station == .x$Station[1]) %>%
              dplyr::select(Station, dist_sea_logged) %>%
              unique() %>%
              pull(dist_sea_logged) %>%
              as.numeric()
            
            # for each tau, for each location, get quantile estimate
            for(q in seq_along(quantiles_to_estimate_bulk)){
              res = rbind(res,tibble(year = yrs,
                                     quantile =  quantiles_to_estimate_bulk[q],
                                     value = obs_q_mod[q,]$quant_mod %>% lapply(predict, tibble(loess_temp_anom = loess_temp_anom_vals, dist_sea_logged = dist_sea_logged_vals)) %>% .[[1]] %>% .$location))
            }
          }else if(mod == "mg"){
            
            temporal_covariates = obs_data %>%
              dplyr::select(year, month, loess_temp_anom, residuals) %>%
              unique() %>%
              arrange(year)
            
            
            
            
            # --- spatial covariate
            dist_sea_logged_vals = obs_data %>%
              filter(Station == .x$Station[1]) %>%
              dplyr::select(Station, dist_sea_logged) %>%
              unique() %>%
              pull(dist_sea_logged) %>%
              as.numeric()
            
            # for each tau, for each location, get quantile estimate
            for(q in seq_along(quantiles_to_estimate_bulk)){
              res = rbind(res,tibble(year = temporal_covariates$year,
                                     quantile =  quantiles_to_estimate_bulk[q],
                                     value = obs_q_mod[q,]$quant_mod %>% lapply(predict, newdata = tibble(value = clim_vals[q], loess_temp_anom = temporal_covariates$loess_temp_anom, 
                                                                                                          residuals = temporal_covariates$residuals,
                                                                                                          dist_sea_logged = dist_sea_logged_vals)) %>% .[[1]] %>% .$location))
            }
          }
          
          res %>% 
            group_by(year) %>%
            group_map(~{
              tibble(year = .x$year[1], quant_spline = list(splinefun(.x$quantile,.x$value,  method = 'monoH.FC')))
            }, .keep = T) %>%
            plyr::rbind.fill() %>%
            as_tibble() %>%
            mutate(Station = .x$Station[1])
          
        }, .keep = T) %>%
        plyr::rbind.fill() %>%
        as_tibble()
      
      # ------- calculate RMSE
      squared_error = test_set %>%
        group_by(Station, year) %>%
        group_map(~{
          
          this_quant_spline = obs_smoothed_quantiles %>%
            filter(Station == .x$Station[1]) %>%
            filter(year == .x$year[1]) %>%
            pull(quant_spline) %>%
            .[[1]]
          
          pred = this_quant_spline(.x$emp_quantile)
          (.x$temp-pred)^2
          
        }, .keep = T) %>%
        unlist()
      sqrt(mean(squared_error))
    }
  }

  
  if(cv_type == "ST"){
    
    set.seed(51964)
    
    num_spatial_folds  = 15
    clutered = spatial_clustering_cv(obs_data %>%
                                       dplyr::select(Station, Long, Lat) %>%
                                       unique(), coords = c(Long, Lat), v = num_spatial_folds)
    
    spatial_folds = c()
    for(i in seq(num_spatial_folds)){
      spatial_folds = rbind(spatial_folds, 
                            assessment(clutered$splits[i][[1]]) %>% 
                              dplyr::select(Station, Long, Lat) %>% 
                              unique %>%
                              mutate(spatial_fold = i))
    }
    
    
    obs_data = obs_data %>%
      left_join(spatial_folds)
    
    obs_data = obs_data %>% mutate(week = lubridate::week(date))
    
    week_chunks = list(c(48, 51, 1, 4, 7),
                       c(49, 52, 2, 5, 8),
                       c(50, 53, 3, 6, 9))
    

    
    obs_data$temporal_fold = 123456789
    obs_data[obs_data$week %in% week_chunks[[1]],]$temporal_fold=1
    obs_data[obs_data$week %in% week_chunks[[2]],]$temporal_fold=2
    obs_data[obs_data$week %in% week_chunks[[3]],]$temporal_fold=3

    for(spat_fold in seq_id){
      for(temp_fold in (obs_data$temporal_fold %>% unique())){
        test_set = obs_data %>% filter(spatial_fold == spat_fold, temporal_fold == temp_fold)
        train_set = obs_data %>% anti_join(test_set)

        tibble(mod, fit_qunat_mod(train_set, test_set, mod)) %>%
          write_csv("output/bulk_cv_res_ST.csv", append = T)
      }
    }
    
  }else{
    num_random_folds = 45
    set.seed(1234567)
    obs_data$random_fold = sample(seq(num_random_folds), size = nrow(obs_data), replace = T)
    
    for(i in seq_id){
      test_set = obs_data %>% filter(random_fold == i)
      train_set = obs_data %>% filter(random_fold != i)
      
      tibble(mod, fit_qunat_mod(train_set, test_set, mod)) %>%
        write_csv("output/bulk_cv_res_random.csv", append = T)
    }
  }
}



# ----------- define spatial clusters

job::job({run_bulk_cv("ma", seq(15), cv_type = "ST")})
job::job({run_bulk_cv("mb", seq(15), cv_type = "ST")})
job::job({run_bulk_cv("mc", seq(15), cv_type = "ST")})
job::job({run_bulk_cv("md", seq(15), cv_type = "ST")})
job::job({run_bulk_cv("me", seq(15), cv_type = "ST")})
job::job({run_bulk_cv("mf", seq(15), cv_type = "ST")})
job::job({run_bulk_cv("mg", seq(15), cv_type = "ST")})

job::job({run_bulk_cv("ma", seq(45), cv_type = "n-fold")})
job::job({run_bulk_cv("mb", seq(45), cv_type = "n-fold")})
job::job({run_bulk_cv("mc", seq(45), cv_type = "n-fold")})
job::job({run_bulk_cv("md", seq(45), cv_type = "n-fold")})
job::job({run_bulk_cv("me", seq(45), cv_type = "n-fold")})
job::job({run_bulk_cv("mf", seq(45), cv_type = "n-fold")})
job::job({run_bulk_cv("mg", seq(45), cv_type = "n-fold")})



# --- calcualtes results and create table

results = rbind(read_csv("~/output/bulk_cv_res_ST.csv",
                         col_names = c("model", "RMSE")) %>%
                  mutate(method = "SV"),
                read_csv("~/output/bulk_cv_res_random.csv",
                         col_names = c("model", "RMSE"))%>%
                  mutate(method = "n-fold"))

models = c("$\\beta^{(\\tau)}_0$",
           "$\\beta^{(\\tau)}_0 + \\beta^{(\\tau)}_1 q^{(\\tau)}_{c}(\\boldsymbol{s})$",
           "$\\beta^{(\\tau)}_0 + \\beta^{(\\tau)}_1 q^{(\\tau)}_{c}(\\boldsymbol{s}) + \\beta^{(\\tau)}_2 M^{\\text{I}}(t)$",
           "$\\beta^{(\\tau)}_0 + \\beta^{(\\tau)}_1 q^{(\\tau)}_{c}(\\boldsymbol{s}) + \\beta^{(\\tau)}_2 M^{\\text{I}}(t) + \\beta^{(\\tau)}_3 C(\\boldsymbol{s})$",
           "$\\beta^{(\\tau)}_0 + \\beta^{(\\tau)}_1 M^{\\text{I}}(t)$",
           "$\\beta^{(\\tau)}_0 + \\beta^{(\\tau)}_1 M^{\\text{I}}(t) + \\beta^{(\\tau)}_2 C(\\boldsymbol{s})$",
           "$\\beta^{(\\tau)}_0 + \\beta^{(\\tau)}_1 q^{(\\tau)}_{c}(\\boldsymbol{s}) + \\beta^{(\\tau)}_2 M^{\\text{I}}(t) + \\beta^{(\\tau)}_3 C(\\boldsymbol{s})  + \\beta^{(\\tau)}_4 M^{\\text{I}}_{r,m}(t)$")

data.frame(name = c('Ma', 'Mb', 'Mc', 'Md', 'Me', 'Mf', 'Mg'), 
           models,
             RMSE_ST = results %>%
               filter(method == "SV") %>%
               group_by(model) %>% 
               summarise(RMSE = mean(RMSE)) %>%
               arrange(model)%>%
             pull(RMSE) %>% signif(4),
           CRPS_ST = results %>%
             filter(method == "n-fold") %>%
             group_by(model) %>% 
             summarise(RMSE = mean(RMSE))%>%
             arrange(model)%>%
             pull(RMSE) %>% signif(4)) %>%
  xtable::xtable(digits = 3) %>%
  print(include.rownames = F, 
        sanitize.colnames.function = identity,
        sanitize.text.function = identity)
