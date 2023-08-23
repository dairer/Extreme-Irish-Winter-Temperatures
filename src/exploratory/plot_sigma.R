# this script calculates and plots scale parameter of the GPD in 2022 and its change since 1950 for three levels of climatic variability. 

rm(list = ls())
gc()
library(tidyverse)
setwd("~/chapter_6/")
source('src/models/marginal/new_gpd_models.R')

ireland_sf <-rnaturalearth::ne_countries(type = "map_units",
                                         scale='large', # resolution of map
                                         returnclass = 'sf', # spatial object, change this to "sp" if needed, this will break the plot though
                                         continent = "europe") %>%
  .[.$name %in% c('Ireland', 'N. Ireland'),]


# colours for plotting
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


# # --- temporal covariates
temporal_covariates_yearly = read_csv("data/processed/obs_data_mintp_winter.csv") %>%
  dplyr::select(year, loess_temp_anom, yearly_residuals, temp_anom_around_irel_resids, temp_anom_north_irel_resids, nao, ao, ao_resid, nao_resid) %>%
  unique()


run_sigma = function(marg_mods, quantile_of_covars, no_flux = F, thresh_qnt){
  
  library(tidyverse)
  setwd("~/chapter_6/")
  source('src/models/marginal/new_gpd_models.R')
  for(qnt in quantile_of_covars){
    
    
    if(!no_flux){
      temporal_covariates_monthly = vroom::vroom("data/processed/obs_data_mintp_winter.csv") %>%
        dplyr::select(year, month, residuals, loess_temp_anom) %>%
        unique() %>%
        mutate(residuals = quantile(residuals, qnt))
      
    }else{
      
      temporal_covariates_monthly = vroom::vroom("data/processed/obs_data_mintp_winter.csv") %>%
        dplyr::select(year, month, residuals,loess_temp_anom) %>%
        unique()
      temporal_covariates_monthly$residuals = 0
      
    }
    
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
    
    
    qunatile_model_fit_95 = readRDS("output/threshold_model/threshold_model_tn_l_95_coast")
    clim_grid$threshold_95_o = predict(qunatile_model_fit_95, newdata = clim_grid)$location
    rm("qunatile_model_fit_95")
    
    qunatile_model_fit_96 = readRDS("output/threshold_model/threshold_model_tn_l_96_coast")
    clim_grid$threshold_96_o= predict(qunatile_model_fit_96, newdata = clim_grid)$location
    rm("qunatile_model_fit_96")
    
    qunatile_model_fit_97 = readRDS("output/threshold_model/threshold_model_tn_l_97_coast")
    clim_grid$threshold_97_o= predict(qunatile_model_fit_97, newdata = clim_grid)$location
    rm("qunatile_model_fit_97")
    
    
    
    for(m in marg_mods){
      if(thresh_qnt == 0.95){
        
        
        true_change = clim_grid %>%
          left_join(read_csv(paste0("data/processed/CLIM_GRID_thresh_exceedance_lambda_mintp_",qnt,".csv"))  %>%
                      filter(year %in% c(1950, 2020, 2022))) %>%
          dplyr::select(id,  year, thresh_exceedance_95, threshold_95_o, dist_sea_logged)
        
        true_change = rbind(true_change %>% mutate(month = 1),
                            true_change %>% mutate(month = 2),
                            true_change %>% mutate(month = 12)) %>%
          left_join(temporal_covariates_monthly)  %>%
          drop_na()
        
        true_change$rl_qnt = 1 - (1/true_change$thresh_exceedance_95)/(100*90.25)
        true_change$threshold = true_change$threshold_95_o
        
      }else if(thresh_qnt == 0.96){
        
        true_change = clim_grid %>%
          left_join(read_csv(paste0("data/processed/CLIM_GRID_thresh_exceedance_lambda_mintp_",qnt,".csv"))  %>%
                      filter(year %in% c(1950, 2020, 2022))) %>%
          dplyr::select(id,  year, thresh_exceedance_96, threshold_96_o, dist_sea_logged)
        
        true_change = rbind(true_change %>% mutate(month = 1),
                            true_change %>% mutate(month = 2),
                            true_change %>% mutate(month = 12)) %>%
          left_join(temporal_covariates_monthly)  %>%
          drop_na()
        
        true_change$rl_qnt = 1 - (1/true_change$thresh_exceedance_96)/(100*90.25)
        true_change$threshold = true_change$threshold_96_o    
        
      }else if(thresh_qnt == 0.97){
        
        true_change = clim_grid %>%
          left_join(read_csv(paste0("data/processed/CLIM_GRID_thresh_exceedance_lambda_mintp_",qnt,".csv"))  %>%
                      filter(year %in% c(1950, 2020, 2022))) %>%
          dplyr::select(id,  year, thresh_exceedance_97, threshold_97_o, dist_sea_logged)
        
        true_change = rbind(true_change %>% mutate(month = 1),
                            true_change %>% mutate(month = 2),
                            true_change %>% mutate(month = 12)) %>%
          left_join(temporal_covariates_monthly)  %>%
          drop_na()
        
        true_change$rl_qnt = 1 - (1/true_change$thresh_exceedance_97)/(100*90.25)
        true_change$threshold = true_change$threshold_97_o        
      }
      
      
      assign(paste0('pars_', m), readRDS(paste0('output/bts_param_est/', m)))
      this_fit = predict_gpd(estimates_pars = get(paste0('pars_', m)), model = m, data = true_change)
      
      
      true_change$scl_est = this_fit$scale
      true_change$shp_est = this_fit$shape
      # true_change$rl = rl_gpd(estimates_pars =get(paste0('pars_', m)), model =  m, data = true_change, rl_quantile = true_change$rl_qnt)
      
      true_change %>%
        # add back in long lat
        left_join(read_csv("data/processed/winter_scales_clim_grid_mintp.csv") %>% 
                    dplyr::select(-c(scales_95, scales_96, scales_97,
                                     threshold_95, threshold_96, threshold_97))) %>%
        dplyr::select(id, Long, Lat, scl_est, year) %>%
        mutate(qnt = qnt) %>%
        unique() %>%
        write_csv(paste0("output/sigma_change/thresh_qnt_",thresh_qnt,"_model_", m, "_qnt_", qnt), append = T)
    }
  }
}


job::job({run_sigma(c("m22"),  c(0.1), thresh_qnt = 0.96)})
job::job({run_sigma(c("m22"),  c(0.5), thresh_qnt = 0.96)})
job::job({run_sigma(c("m22"),  c(0.9), thresh_qnt = 0.96)})


true = rbind(read_csv(paste0("output/sigma_change/thresh_qnt_",0.96,"_model_", 'm22', "_qnt_", 0.1),
                      col_names = c("id", "Long", "Lat", "sigma", "year", "qnt")),
             read_csv(paste0("output/sigma_change/thresh_qnt_",0.96,"_model_", 'm22', "_qnt_", 0.5),
                      col_names = c("id", "Long", "Lat", "sigma", "year", "qnt")),
             read_csv(paste0("output/sigma_change/thresh_qnt_",0.96,"_model_", 'm22', "_qnt_", 0.9),
                      col_names = c("id", "Long", "Lat", "sigma", "year", "qnt")))



gridExtra::grid.arrange(true %>%
                          filter(year == 2022) %>%
                          dplyr::select(-year) %>% 
                          unique() %>%
                          left_join(
                            true %>%
                              filter(year == 1950) %>%
                              dplyr::select(-year) %>%
                              unique() %>%
                              rename(old_sig = sigma)) %>%
                          mutate(sig_change = sigma - old_sig) %>%
                          ggplot()+
                          geom_point(aes(Long, Lat, col = sigma))+
                          facet_wrap(~qnt) + 
                          geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
                          scale_color_gradientn( colours=my_pal)+
                          labs(col = expression(sigma[o]))+
                          scale_x_continuous(breaks = c(-10, -8, -6))+
                          theme_minimal()+
                          theme(axis.text.x  = element_blank(),
                                axis.text.y  = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.title.x = element_blank(),
                                axis.title.y = element_blank()),
                        true %>%
                          filter(year == 2022) %>%
                          dplyr::select(-year) %>% 
                          unique() %>%
                          left_join(
                            true %>%
                              filter(year == 1950) %>%
                              dplyr::select(-year) %>%
                              unique() %>%
                              rename(old_sig = sigma)) %>%
                          mutate(sig_change = sigma - old_sig) %>%
                          ggplot()+
                          geom_point(aes(Long, Lat, col = sig_change))+
                          facet_wrap(~qnt) + 
                          geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
                          scale_color_gradientn( colours=my_pal)+
                          labs(col = expression(paste(nabla, sigma[o])))+
                          scale_x_continuous(breaks = c(-10, -8, -6))+
                          theme_minimal()+
                          theme(axis.text.x  = element_blank(),
                                axis.text.y  = element_blank(),
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                axis.title.x = element_blank(),
                                axis.title.y = element_blank()), nrow = 2)


ggsave(plt, filename = 'output/figs/scale_changing.pdf', height = 4.8, width = 8)

