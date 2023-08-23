# this script performs cross validation for gpd models
gc()
rm(list = ls())
source('src/models/marginal/gpd_models.R')

library(tidyverse)
library(spatialsample)
library(scoringRules)
library(job)

# colour palatte 
my_pal = c( 
  '#062c30', # extra dark 
  '#003f5c',
  '#2f4b7c',
  '#665191',
  '#a05195',
  '#d45087',
  '#f95d6a',
  '#ff7c43',
  '#ffa600')

obs_data = read_csv("data/processed/obs_data_mintp_winter.csv")



# map of ireland sf
ireland_sf <-rnaturalearth::ne_countries(type = "map_units",
                                         scale='large', # resolution of map
                                         returnclass = 'sf', # spatial object
                                         continent = "europe") %>%
  .[.$name %in% c('Ireland', 'N. Ireland'),]

# -------- CREATE CV Folds
set.seed(51964)

num_spatial_folds  = 15
clutered = spatial_clustering_cv(obs_data %>%
                                   dplyr::select(Station, Long, Lat) %>%
                                   unique(), coords = c(Long, Lat), 
                                 v = num_spatial_folds)

spatial_folds = c()
for(i in seq(num_spatial_folds)){
  spatial_folds = rbind(spatial_folds, 
                        assessment(clutered$splits[i][[1]]) %>% 
                          dplyr::select(Station, Long, Lat) %>% 
                          unique %>%
                          mutate(spatial_fold = i))
}

week_chunks = list(c(48, 51, 1, 4, 7),
                   c(49, 52, 2, 5, 8),
                   c(50, 53, 3, 6, 9))

obs_data$temporal_fold = 123456789
obs_data[obs_data$week %in% week_chunks[[1]],]$temporal_fold=1
obs_data[obs_data$week %in% week_chunks[[2]],]$temporal_fold=2
obs_data[obs_data$week %in% week_chunks[[3]],]$temporal_fold=3


run_cv = function(thresh_qnt,cv_method, models, obs_data, spatial_folds, num_spatial_folds = "NA", week_chunks="NA"){
  setwd("~/chapter_6/")
  
  source('src/models/marginal/gpd_model.R')
  
  my_mods = models
  get_metrics = function(orig_dat, excess, quant, threshold,  pred_scale, pred_shape){
    likelihood_vals = evd::dgpd(x = excess, loc = 0, scale = pred_scale, shape = pred_shape, log = T)
    ll = sum(likelihood_vals[!is.infinite(likelihood_vals)])
    ll_standardised = sum(likelihood_vals[!is.infinite(likelihood_vals)])/length(likelihood_vals[!is.infinite(likelihood_vals)])
    
    x = orig_dat
    y = evd::qgpd(p = quant, loc = threshold, scale = pred_scale, shape = pred_shape[1])
    
    my_rmse = sqrt(mean((x-y)^2))
    scoring = crps(y = excess, family = "gpd", location = 0, scale = pred_scale, shape = pred_shape[1], mass = 0) %>% mean
    return(paste0(ll, ",", ll_standardised, ",",my_rmse, ",", scoring))
  }
  
  
  if(cv_method == "spatial-temporal"){
    
    
    if(thresh_qnt == 0.95){
      obs_data$threshold =  obs_data$threshold_tn_l_o_95_w_coast
      obs_data$scales_logged =  obs_data$scales_95_logged
      
    }else if(thresh_qnt == 0.96){
      obs_data$threshold =  obs_data$threshold_tn_l_o_96_w_coast
      obs_data$scales_logged =  obs_data$scales_96_logged
      
    }else if(thresh_qnt == 0.97){
      obs_data$threshold =  obs_data$threshold_tn_l_o_97_w_coast
      obs_data$scales_logged =  obs_data$scales_97_logged
      
    }
    
    
    
    extreme_data = obs_data %>%
      mutate(excess = temp - threshold) %>%
      filter(excess > 0) %>%
      left_join(spatial_folds) %>%
      group_by(year, Station) %>%
      group_map(~{
        .x %>%
          mutate(quant = rank(excess)/(length(excess)+1))
      },.keep=T) %>%
      plyr::rbind.fill() %>%
      as_tibble() %>%
      ungroup()

    for(s in seq(num_spatial_folds)){
      for(t in week_chunks){
        print(s)
        print(t)

        test = extreme_data %>% 
          filter(spatial_fold == s) %>%
          filter(week %in% t)
        
        train = extreme_data %>%
          anti_join(test)
        
        for(m in my_mods){
          est_pars = fit_gpd(data = train, model = m) 
          pred = est_pars %>% predict_gpd(model = m, data = test)
          write.table(paste0(m, ",", get_metrics(test$temp, test$excess, test$quant, test$threshold,  pred$scale, pred$shape[1])),
                      file = paste0('output/cv_res/paper_spatio_temporal_cv_thresh_qnt_',thresh_qnt,'.csv'),
                      sep = ",", append = TRUE, quote = FALSE,
                      col.names = FALSE, row.names = FALSE)
        }
      }
    }
  }
  
  if(cv_method == "nfold"){

    if(thresh_qnt == 0.95){
      obs_data$threshold =  obs_data$threshold_tn_l_o_95_w_coast
      obs_data$scales_logged =  obs_data$scales_95_logged
      
    }else if(thresh_qnt == 0.96){
      obs_data$threshold =  obs_data$threshold_tn_l_o_96_w_coast
      obs_data$scales_logged =  obs_data$scales_96_logged
      
    }else if(thresh_qnt == 0.97){
      obs_data$threshold =  obs_data$threshold_tn_l_o_97_w_coast
      obs_data$scales_logged =  obs_data$scales_97_logged
    }

    extreme_data = obs_data %>%
      mutate(excess = temp - threshold) %>%
      filter(excess > 0) %>%
      group_by(year, Station) %>%
      group_map(~{
        
        .x %>%
          mutate(quant = rank(excess)/(length(excess)+1))
        
      },.keep=T) %>%
      plyr::rbind.fill() %>%
      as_tibble()%>%
      ungroup()
    
    #cv_method = 'nfold'
    # ---------- define random folds
    num_random_folds = 45
    set.seed(1234567)
    extreme_data$random_fold = sample(seq(num_random_folds), size = nrow(extreme_data), replace = T)
    
    for(i in seq(num_random_folds)){
      
      test = extreme_data %>% filter(random_fold == i)
      train = extreme_data %>% filter(random_fold != i)
      
      for(m in my_mods){
        est_pars = fit_gpd(data = train, model = m) 
        pred = est_pars %>% predict_gpd(model = m, data = test)
        write.table(paste0(m, ",", get_metrics(test$temp, test$excess, test$quant, test$threshold,  pred$scale, pred$shape[1])),
                    file = paste0('output/cv_res/n-fold_thresh_qnt_',thresh_qnt,'.csv'),
                    sep = ",", append = TRUE, quote = FALSE,
                    col.names = FALSE, row.names = FALSE)
      }
    }
  }
}    



job::job({run_cv(0.96, 'spatial-temporal', models = paste0('m',seq(1,3)), obs_data, spatial_folds, num_spatial_folds, week_chunks)}, import = c("obs_data", 'run_cv', 'spatial_folds', 'num_spatial_folds', 'week_chunks'))
job::job({run_cv(0.96, 'spatial-temporal', models = paste0('a',seq(1,13)), obs_data, spatial_folds, num_spatial_folds, week_chunks)}, import = c("obs_data", 'run_cv', 'spatial_folds', 'num_spatial_folds', 'week_chunks'))
job::job({run_cv(0.96, 'spatial-temporal', models = paste0('b',seq(1,13)), obs_data, spatial_folds, num_spatial_folds, week_chunks)}, import = c("obs_data", 'run_cv', 'spatial_folds', 'num_spatial_folds', 'week_chunks'))
job::job({run_cv(0.96, 'spatial-temporal', models = paste0('c',seq(1,13)), obs_data, spatial_folds, num_spatial_folds, week_chunks)}, import = c("obs_data", 'run_cv', 'spatial_folds', 'num_spatial_folds', 'week_chunks'))
job::job({run_cv(0.96, 'spatial-temporal', models = paste0('d',seq(1,13)), obs_data, spatial_folds, num_spatial_folds, week_chunks)}, import = c("obs_data", 'run_cv', 'spatial_folds', 'num_spatial_folds', 'week_chunks'))
job::job({run_cv(0.96, 'spatial-temporal', models = paste0('e',seq(1,13)), obs_data, spatial_folds, num_spatial_folds, week_chunks)}, import = c("obs_data", 'run_cv', 'spatial_folds', 'num_spatial_folds', 'week_chunks'))
job::job({run_cv(0.96, 'spatial-temporal', models = paste0('f',seq(1,13)), obs_data, spatial_folds, num_spatial_folds, week_chunks)}, import = c("obs_data", 'run_cv', 'spatial_folds', 'num_spatial_folds', 'week_chunks'))
job::job({run_cv(0.96, 'spatial-temporal', models = paste0('g',seq(1,13)), obs_data, spatial_folds, num_spatial_folds, week_chunks)}, import = c("obs_data", 'run_cv', 'spatial_folds', 'num_spatial_folds', 'week_chunks'))


job::job({run_cv(0.96, 'nfold', models = paste0('m',seq(1,3)), obs_data)}, import = c("obs_data", 'run_cv'))
job::job({run_cv(0.96, 'nfold', models = paste0('a',seq(1,13)), obs_data)}, import = c("obs_data", 'run_cv'))
job::job({run_cv(0.96, 'nfold', models = paste0('b',seq(1,13)), obs_data)}, import = c("obs_data", 'run_cv'))
job::job({run_cv(0.96, 'nfold', models = paste0('c',seq(1,13)), obs_data)}, import = c("obs_data", 'run_cv'))
job::job({run_cv(0.96, 'nfold', models = paste0('d',seq(1,13)), obs_data)}, import = c("obs_data", 'run_cv'))
job::job({run_cv(0.96, 'nfold', models = paste0('e',seq(1,13)), obs_data)}, import = c("obs_data", 'run_cv'))
job::job({run_cv(0.96, 'nfold', models = paste0('f',seq(1,13)), obs_data)}, import = c("obs_data", 'run_cv'))
job::job({run_cv(0.96, 'nfold', models = paste0('g',seq(1,13)), obs_data)}, import = c("obs_data", 'run_cv'))


# produce table of results

results = rbind(read_csv("~/chapter_6/output/cv_res/paper_spatio_temporal_cv_thresh_qnt_0.96.csv",
                         col_names = c("model", "LL", "LLs", "RMSE", "CRPS")) %>%
                  mutate(method = "SV"),
                read_csv("~/chapter_6/output/cv_res/paper_n-fold_thresh_qnt_0.96.csv",
                         col_names = c("model", "LL", "LLs", "RMSE", "CRPS"))%>%
                  mutate(method = "n-fold"))


results$model = factor(results$model, levels = c(paste0("m", c(1,2,3)), 
                                                 paste0("a", seq(1,13)),
                                                 paste0("b", seq(1,13)),
                                                 paste0("c", seq(1,13)),
                                                 paste0("d", seq(1,13)),
                                                 paste0("e", seq(1,13)),
                                                 paste0("f", seq(1,13)),
                                                 paste0("g", seq(1,13))))
                                                

models = c("$\\beta_0$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s})$",
           "$\\beta_1C(\\boldsymbol{s})$",
           "$\\beta_1M^{I}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2M^{I}_{r,m}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2M^{I}_{r}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2M^{G}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2M^{N}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2\\text{AO}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2\\text{AO}_{r}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2\\text{AO}_{m}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2\\text{AO}_{r,m}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2\\text{NAO}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2\\text{NAO}_{r}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2\\text{NAO}_{m}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2\\text{NAO}_{r,m}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}_{r,m}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}_{r}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{G}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{N}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3\\text{AO}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3\\text{AO}_{r}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3\\text{AO}_{m}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3\\text{AO}_{r,m}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3\\text{NAO}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3\\text{NAO}_{r}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3\\text{NAO}_{m}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3\\text{NAO}_{r,m}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s})$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4M^{I}_{r,m}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4M^{I}_{r}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4M^{G}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4M^{N}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{AO}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{AO}_{r}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{AO}_{m}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{AO}_{r,m}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{NAO}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{NAO}_{r}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{NAO}_{m}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{NAO}_{r,m}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s})$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4M^{I}_{r,m}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4M^{I}_{r}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4M^{G}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4M^{N}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{AO}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{AO}_{r}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{AO}_{m}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{AO}_{r,m}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{NAO}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{NAO}_{r}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{NAO}_{m}(t)$",
           "$\\beta_1\\sigma_c(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{NAO}_{r,m}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) \\beta_3M^{I}_{r,m}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) \\beta_3M^{I}_{r}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) \\beta_3M^{G}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) \\beta_3M^{N}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) \\beta_3\\text{AO}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) \\beta_3\\text{AO}_{r}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) \\beta_3\\text{AO}_{m}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) \\beta_3\\text{AO}_{r,m}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) \\beta_3\\text{NAO}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) \\beta_3\\text{NAO}_{r}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) \\beta_3\\text{NAO}_{m}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) \\beta_3\\text{NAO}_{r,m}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) $",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4M^{I}_{r,m}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4M^{I}_{r}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4M^{G}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4M^{N}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{AO}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{AO}_{r}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{AO}_{m}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{AO}_{r,m}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{NAO}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{NAO}_{r}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{NAO}_{m}(t)$",
           "$\\beta_1C(\\boldsymbol{s}) + \\beta_2M^{I}(t) + \\beta_3M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{NAO}_{r,m}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2M^{I}(t)C(\\boldsymbol{s}) $",
           "$\\beta_1M^{I}(t) + \\beta_2M^{I}(t)C(\\boldsymbol{s}) + \\beta_4M^{I}_{r,m}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2M^{I}(t)C(\\boldsymbol{s}) + \\beta_4M^{I}_{r}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2M^{I}(t)C(\\boldsymbol{s}) + \\beta_4M^{G}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2M^{I}(t)C(\\boldsymbol{s}) + \\beta_4M^{N}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{AO}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{AO}_{r}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{AO}_{m}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{AO}_{r,m}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{NAO}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{NAO}_{r}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{NAO}_{m}(t)$",
           "$\\beta_1M^{I}(t) + \\beta_2M^{I}(t)C(\\boldsymbol{s}) + \\beta_4\\text{NAO}_{r,m}(t)$")


name = results %>%
  filter(method == "SV") %>%
  group_by(model) %>% 
  summarise(RMSE = mean(RMSE), CRPS = mean(CRPS), LL = mean(LL)) %>%
  pull(model)


data.frame(frst = c(c("", "", ""),
                    c("A", rep("", 12)),
                    c("B", rep("", 12)),
                    c("C", rep("", 12)),
                    c("D", rep("", 12)),
                    c("E", rep("", 12)),
                    c("F", rep("", 12)),
                    c("G", rep("", 12))),
           sec = as.character( c(c(1,2,3), (rep(seq(1,13), 7)))),
           name, models,
           RMSE_ST = results %>%
             filter(method == "SV") %>%
             group_by(model) %>% 
             summarise(RMSE = mean(RMSE), CRPS = mean(CRPS), LL = mean(LL)) %>%
             pull(RMSE) %>% signif(4),
           CRPS_ST = results %>%
             filter(method == "SV") %>%
             group_by(model) %>% 
             summarise(RMSE = mean(RMSE), CRPS = mean(CRPS), LL = mean(LL)) %>%
             pull(CRPS) %>% signif(4),
           RMSE_n = results %>%
             filter(method == "n-fold") %>%
             group_by(model) %>% 
             summarise(RMSE = mean(RMSE), CRPS = mean(CRPS), LL = mean(LL)) %>%
             pull(RMSE) %>% signif(4),
           CRPS_n = results %>%
             filter(method == "n-fold") %>%
             group_by(model) %>% 
             summarise(RMSE = mean(RMSE), CRPS = mean(CRPS), LL = mean(LL)) %>%
             pull(CRPS) %>% signif(4)) %>%
  dplyr::select(-name) %>%
  xtable::xtable(digits = 3) %>%
  print(include.rownames = F, 
        sanitize.colnames.function = identity,
        sanitize.text.function = identity)
