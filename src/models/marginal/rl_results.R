#This script calculates marginal return levels and associated uncertainties from bootstrapped data sets for diffferet thresholds and different levels of SCV


rm(list=ls())
source('src/models/marginal/gpd_models.R')
library(tidyverse)

library(rnaturalearth)

# map for plotting
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

# # --- point estimate
temporal_covariates_yearly = read_csv("data/processed/obs_data_mintp_winter.csv") %>%
  dplyr::select(year, loess_temp_anom, yearly_residuals, temp_anom_around_irel_resids, temp_anom_north_irel_resids, nao, ao, ao_resid, nao_resid) %>%
  unique()



run_rl = function(marg_mods, quantile_of_covars, no_flux = F, thresh_qnt){
  
  # thresh_qnt = 0.96
  library(tidyverse)
  source('src/models/marginal/gpd_models.R')
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

    dist_sea = read_csv("~/data/processed/sites_clim_sea_dist.csv") %>%
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
      true_change$rl = rl_gpd(estimates_pars =get(paste0('pars_', m)), model =  m, data = true_change, rl_quantile = true_change$rl_qnt)
      
      
      true_change %>%
        filter(is.nan(rl))

      true_change %>%
        # add back in long lat
        left_join(read_csv("data/processed/winter_scales_clim_grid_mintp.csv") %>% 
                    dplyr::select(-c(scales_95, scales_96, scales_97,
                                     threshold_95, threshold_96, threshold_97))) %>%
        dplyr::select(id, Long, Lat, rl, year) %>%
        mutate(qnt = qnt) %>%
        unique() %>%
        write_csv(paste0("output/rl_marg/thresh_qnt_",thresh_qnt,"_model_", m, "_qnt_", qnt), append = T)
    }
  }
}


job::job({run_rl(c("E2"),  c(0.1), thresh_qnt = 0.96)})
job::job({run_rl(c("E2"),  c(0.5), thresh_qnt = 0.96)})
job::job({run_rl(c("E2"),  c(0.9), thresh_qnt = 0.96)})
job::job({run_rl(c("E2"),  c(0.5), thresh_qnt = 0.95)})
job::job({run_rl(c("E2"),  c(0.5), thresh_qnt = 0.97)})
job::job({run_rl(c("E1"),  c(0.5), thresh_qnt = 0.96)})


run_rl_uncert = function(marg_mods, bts_range, quantile_of_covars,thresh_qnt, no_flux = F){
  library(tidyverse)
  source('src/models/marginal/gpd_models.R')
  
  
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
    
    for(m in marg_mods){
      for(bts in bts_range){
        
        
        
        if(thresh_qnt == 0.95){
          
          blk_md = clim_grid %>%
            left_join(read_csv(paste0("output/bts_thresh_ex_lambda/clim_grid_bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", m, "_bts_",bts,"_monthly_clim_covar_",qnt,".csv"),
                               col_names = c('bts', 'id', 'year', 'thresh_exceedance_95'))  %>%
                        filter(year %in% c(1950, 2020, 2022))) %>%
            dplyr::select(id,  year, thresh_exceedance_95) %>%
            left_join((readRDS(paste0("output/bootstrapped_thresh_clim/thresh_qnt_",thresh_qnt, "_model_", m,"_bts_", bts))))
          
          blk_md = rbind(blk_md %>% mutate(month = 1),
                              blk_md %>% mutate(month = 2),
                              blk_md %>% mutate(month = 12)) %>%
            left_join(temporal_covariates_monthly)  %>%
            drop_na()
          
          blk_md$rl_qnt = 1 - (1/blk_md$thresh_exceedance_95)/(100*90.25)
          blk_md$threshold = blk_md$threshold_tn_l_o_95_w_coast
          
        }else if(thresh_qnt == 0.96){
          
          blk_md = clim_grid %>%
            left_join(read_csv(paste0("output/bts_thresh_ex_lambda/clim_grid_bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", m, "_bts_",bts,"_monthly_clim_covar_",qnt,".csv"),
                               col_names = c('bts', 'id', 'year', 'thresh_exceedance_96'))  %>%
                        filter(year %in% c(1950, 2020, 2022))) %>%
            dplyr::select(id,  year, thresh_exceedance_96) %>%
            left_join((readRDS(paste0("output/bootstrapped_thresh_clim/thresh_qnt_",thresh_qnt, "_model_", m,"_bts_", bts))))
          
          blk_md = rbind(blk_md %>% mutate(month = 1),
                              blk_md %>% mutate(month = 2),
                              blk_md %>% mutate(month = 12)) %>%
            left_join(temporal_covariates_monthly)  %>%
            drop_na()
          
          blk_md$rl_qnt = 1 - (1/blk_md$thresh_exceedance_96)/(100*90.25)
          blk_md$threshold = blk_md$threshold_tn_l_o_96_w_coast
          
        }else if(thresh_qnt == 0.97){
          
          blk_md = clim_grid %>%
            left_join(read_csv(paste0("output/bts_thresh_ex_lambda/clim_grid_bootstrapped_threshold_exceedence_lambda_thresh_qnt_",thresh_qnt,"_marg_mod_", m, "_bts_",bts,"_monthly_clim_covar_",qnt,".csv"),
                               col_names = c('bts', 'id', 'year', 'thresh_exceedance_97'))  %>%
                        filter(year %in% c(1950, 2020, 2022))) %>%
            dplyr::select(id,  year, thresh_exceedance_97) %>%
            left_join((readRDS(paste0("output/bootstrapped_thresh_clim/thresh_qnt_",thresh_qnt, "_model_", m,"_bts_", bts))))
          
          blk_md = rbind(blk_md %>% mutate(month = 1),
                         blk_md %>% mutate(month = 2),
                         blk_md %>% mutate(month = 12)) %>%
            left_join(temporal_covariates_monthly)  %>%
            drop_na()
          
          blk_md$rl_qnt = 1 - (1/blk_md$thresh_exceedance_97)/(100*90.25)
          blk_md$threshold = blk_md$threshold_tn_l_o_97_w_coast
        }
        


        assign(paste0('pars_', m), (read.table(paste0("output/bts_param_est/thresh_qnt_",thresh_qnt,"_corrected_",m,".csv")) %>% filter(V1 == bts) %>% unlist %>% .[-1] %>% as.numeric()))
        this_fit = predict_gpd(estimates_pars = get(paste0('pars_', m)), model = m, data = blk_md)
        
        blk_md$scl_est = this_fit$scale
        blk_md$shp_est = this_fit$shape
        blk_md$rl = rl_gpd(estimates_pars =get(paste0('pars_', m)), model =  m, data = blk_md, rl_quantile = blk_md$rl_qnt)

        blk_md %>%
          left_join(read_csv("data/processed/winter_scales_clim_grid_mintp.csv") %>% 
                      dplyr::select(-c(scales_95, scales_96, scales_97,
                                       threshold_95, threshold_96, threshold_97))) %>%
          dplyr::select(id, Long, Lat, rl, year) %>%
          mutate(bts = bts,
                 qnt = qnt) %>%
          unique() %>%
          write_csv(paste0("output/rl_marg_bts/thresh_qnt_",thresh_qnt,"_model_", m, "_qnt_", qnt), append = T)
      }
    }
  }
}
  

# 
# 
job::job({run_rl_uncert(c("E2"), seq(200), 0.1, 0.96)})
job::job({run_rl_uncert(c("E2"), seq(200), 0.9, 0.96)})
job::job({run_rl_uncert(c("E1"), seq(200), 0.5, 0.96)})
job::job({run_rl_uncert(c("E2"), seq(200), 0.5, 0.95)})
job::job({run_rl_uncert(c("E2"), seq(200), 0.5, 0.97)})




true = rbind(read_csv(paste0("output/rl_marg/thresh_qnt_",0.96,"_model_", 'E2', "_qnt_", 0.1),
                      col_names = c("id", "Long", "Lat", "rl", "year", "qnt")),
             read_csv(paste0("output/rl_marg/thresh_qnt_",0.96,"_model_", 'E2', "_qnt_", 0.5),
                      col_names = c("id", "Long", "Lat", "rl", "year", "qnt")),
             read_csv(paste0("output/rl_marg/thresh_qnt_",0.96,"_model_", 'E2', "_qnt_", 0.9),
                      col_names = c("id", "Long", "Lat", "rl", "year", "qnt")))

bts = rbind(read_csv(paste0("output/rl_marg_bts/thresh_qnt_",0.96,"_model_", 'E2', "_qnt_", 0.1),
                     col_names = c("id", "Long", "Lat", "rl", "year", 'bts', "qnt")),
            read_csv(paste0("output/rl_marg_bts/thresh_qnt_",0.96,"_model_", 'E2', "_qnt_", 0.5),
                     col_names = c("id", "Long", "Lat", "rl", "year", 'bts', "qnt")),
            read_csv(paste0("output/rl_marg_bts/thresh_qnt_",0.96,"_model_", 'E2', "_qnt_", 0.9),
                     col_names = c("id", "Long", "Lat", "rl", "year", 'bts', "qnt")))

get_diff_rl_true = function(data){
  data %>%
    filter(year == 1950) %>%
    rename(rl_100_1950 = rl) %>%
    dplyr::select(Long, Lat, id, qnt, rl_100_1950) %>%
    unique() %>%
    left_join(data %>%
                filter(year == 2022) %>%
                rename(rl_100_2022 = rl) %>%
                unique() %>%
                dplyr::select(Long, Lat, id, qnt, rl_100_2022))%>%
    mutate(rl_100_diff = (-rl_100_2022) - (-rl_100_1950))
}


get_diff_rl_bts = function(data){
  data %>%
    filter(year == 1950) %>%
    rename(rl_100_1950 = rl) %>%
    dplyr::select(Long, Lat, id, qnt, bts, rl_100_1950) %>%
    unique() %>%
    left_join(data %>%
                filter(year == 2022) %>%
                rename(rl_100_2022 = rl) %>%
                unique() %>%
                dplyr::select(Long, Lat, id, qnt, bts, rl_100_2022))%>%
    mutate(rl_100_diff = (-rl_100_2022) - (-rl_100_1950))
}




plt = gridExtra::grid.arrange(
  get_diff_rl_true(true) %>%
    filter(qnt == 0.1) %>%
    ggplot()+
    geom_point(aes(Long, Lat, col = rl_100_diff))+
    geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
    scale_color_gradientn( colours=my_pal)+
    labs(col = expression(paste('°C')))+
    scale_x_continuous(breaks = c(-10, -8, -6))+
    theme_minimal()+
    theme(axis.text.x  = element_blank(),
          axis.text.y  = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()),
  get_diff_rl_true(true) %>%
    filter(qnt == 0.5) %>%
    ggplot()+
    geom_point(aes(Long, Lat, col = rl_100_diff))+
    geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
    scale_color_gradientn( colours=my_pal)+
    labs(col = expression(paste('°C')))+
    scale_x_continuous(breaks = c(-10, -8, -6))+
    theme_minimal()+
    theme(axis.text.x  = element_blank(),
          axis.text.y  = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()),
  get_diff_rl_true(true) %>%
    filter(qnt == 0.9) %>%
    ggplot()+
    geom_point(aes(Long, Lat, col = rl_100_diff))+
    geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
    scale_color_gradientn( colours=my_pal)+
    labs(col = expression(paste('°C')))+
    scale_x_continuous(breaks = c(-10, -8, -6))+
    theme_minimal()+
    theme(axis.text.x  = element_blank(),
          axis.text.y  = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()), nrow = 1)

ggsave(plt, filename = 'output/figs/rl_E2_96.pdf', height = 2.5, width = 8)


plt = true %>%
    filter(year == 2022) %>%
    ggplot()+
    geom_point(aes(Long, Lat, col = -rl))+
    geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
    scale_color_gradientn(colours=my_pal, limits = c(-15,-2.8))+
    facet_wrap(~qnt, nrow = 1)+
    labs(col = expression(paste('°C')))+
    scale_x_continuous(breaks = c(-10, -8, -6))+
    theme_minimal()+
    theme(strip.text = element_blank(),
      axis.text.x  = element_blank(),
      axis.text.y  = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank())


ggsave(plt, filename = 'output/figs/rl_E2_96_point_est.pdf', height = 2.5, width = 8)


plt = gridExtra::grid.arrange(read_csv(paste0("output/rl_marg/thresh_qnt_",0.96,"_model_", 'E1', "_qnt_", 0.5),
                                       col_names = c("id", "Long", "Lat", "rl", "year", "qnt")) %>%
                               filter(year == 2022) %>%
                               ggplot()+
                               geom_point(aes(Long, Lat, col = -rl))+
                               geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
                               scale_color_gradientn(colours=my_pal, breaks = -c(8, 11, 14, 17))+
                               # scale_color_gradientn(colours=my_pal, limits = c(-15,-2.8))+
                               facet_wrap(~qnt, nrow = 1)+
                               labs(col = expression(paste('°C')))+
                               scale_x_continuous(breaks = c(-10, -8, -6))+
                               theme_minimal()+
                               theme(strip.text = element_blank(),
                                     axis.text.x  = element_blank(),
                                     axis.text.y  = element_blank(),
                                     panel.grid.major = element_blank(),
                                     panel.grid.minor = element_blank(),
                                     axis.title.x = element_blank(),
                                     axis.title.y = element_blank()),
                              read_csv(paste0("output/rl_marg/thresh_qnt_",0.96,"_model_", 'E1', "_qnt_", 0.5),
                                       col_names = c("id", "Long", "Lat", "rl", "year", "qnt")) %>%
                                get_diff_rl_true %>%
                                ggplot()+
                                geom_point(aes(Long, Lat, col = rl_100_diff))+
                                geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
                                scale_color_gradientn(colours=my_pal)+
                                facet_wrap(~qnt, nrow = 1)+
                                labs(col = expression(paste('°C')))+
                                scale_x_continuous(breaks = c(-10, -8, -6))+
                                theme_minimal()+
                                theme(strip.text = element_blank(),
                                      axis.text.x  = element_blank(),
                                      axis.text.y  = element_blank(),
                                      panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      axis.title.x = element_blank(),
                                      axis.title.y = element_blank()), widths = c(1.01, 1), nrow = 1)

ggsave(plt, filename = 'output/figs/rl_E1_96_point_est.pdf', height = 2.5, width = 5.34)


# ---- DIFFERNECE IN 100 year level from 0.95 and 0.97 to 0.96
true = rbind(read_csv(paste0("output/rl_marg/thresh_qnt_",0.95,"_model_", 'E2', "_qnt_", 0.5),
                      col_names = c("id", "Long", "Lat", "rl", "year", "qnt")) %>%
               mutate(thresh = 'a'),
             read_csv(paste0("output/rl_marg/thresh_qnt_",0.97,"_model_", 'E2', "_qnt_", 0.5),
                      col_names = c("id", "Long", "Lat", "rl", "year", "qnt")) %>%
               mutate(thresh = 'b'),
             read_csv(paste0("output/rl_marg/thresh_qnt_",0.96,"_model_", 'E2', "_qnt_", 0.5),
                      col_names = c("id", "Long", "Lat", "rl", "year", "qnt")) %>%
               mutate(thresh = 'tr'))

bts = read_csv(paste0("output/rl_marg_bts/thresh_qnt_",0.96,"_model_", 'E2', "_qnt_", 0.5),
               col_names = c("id", "Long", "Lat", "rl", "year", 'bts', "qnt"))



plt = gridExtra::grid.arrange(
  bts %>%
    filter(year == 2022) %>%
    group_by(id, Lat, Long) %>%
    summarise(upper = quantile(rl, 0.975, na.rm = T),
              lower = quantile(rl, 0.025, na.rm = T)) %>%
    mutate(diff = (upper-lower)/2) %>%
    ggplot()+
    geom_point(aes(Long, Lat, col = diff))+
    geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
    labs(col = expression(paste('°C')))+
    scale_x_continuous(breaks = c(-10, -8, -6))+
    scale_color_gradientn( colours=my_pal, breaks = c(1, 1.25, 1.5))+
    theme_minimal()+
    theme(strip.text = element_blank(),
          axis.text.x  = element_blank(),
          axis.text.y  = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()),rbind(true %>%
        filter(year == 2022) %>%
        filter(thresh == 'tr') %>%
        dplyr::select(-thresh) %>%
        rename(rl_true = rl) %>%
        left_join(true %>%
                    filter(year == 2022) %>%
                    filter(thresh == 'a') %>%
                    dplyr::select(-thresh) %>%
                    rename(rl_a = rl)) %>%
        mutate(diff = rl_a - rl_true) %>%
        dplyr::select(-rl_a)%>%
        mutate(labs = "diff_95"),
      true %>%
        filter(year == 2022) %>%
        filter(thresh == 'tr') %>%
        dplyr::select(-thresh) %>%
        rename(rl_true = rl) %>%
        left_join(true %>%
                    filter(year == 2022) %>%
                    filter(thresh == 'b') %>%
                    dplyr::select(-thresh) %>%
                    rename(rl_b = rl)) %>%
        mutate(diff = rl_b - rl_true) %>%
        dplyr::select(-rl_b)%>%
        mutate(labs = "diff_97")) %>%
  ggplot()+
  geom_point(aes(Long, Lat, col = diff))+
  geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
  labs(col = expression(paste('°C')))+
  facet_wrap(~labs)+
  scale_x_continuous(breaks = c(-10, -8, -6))+
  scale_color_gradientn( colours=my_pal)+
  theme_minimal()+
  theme(strip.text = element_blank(),
        axis.text.x  = element_blank(),
        axis.text.y  = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank()), nrow = 1, widths = c(1,1.65))

ggsave(plt, filename = 'output/figs/diff_in_rl_diff_thresholds.pdf', height = 2.5, width = 8)