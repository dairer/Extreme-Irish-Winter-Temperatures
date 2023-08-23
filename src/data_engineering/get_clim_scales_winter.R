# calculate climate scale covairate at different threshold levels

gc()
rm(list = ls())
library(tidyverse)
library(vroom)
library(mgcv)
library(raster) # package for netcdf manipulation

var_modelled = "mintp"

# months to include in the analysis
months_to_study = c(12,1,2) # winter months we consider

refit_clim_shape_param = T
my_pal = c( 
  '#062c30',
  '#003f5c',
  '#2f4b7c',
  '#665191',
  '#a05195',
  '#d45087',
  '#f95d6a',
  '#ff7c43',
  '#ffa600')

# map of ireland sf
ireland_sf <-rnaturalearth::ne_countries(type = "map_units",
                                         scale='large', # resolution of map
                                         returnclass = 'sf', # spatial object
                                         continent = "europe") %>%
  .[.$name %in% c('Ireland', 'N. Ireland'),]

if(var_modelled == "mintp"){
  clim_data = vroom("~/chapter_6/data/processed/mintp_r12i1p1-ICHEC-EC-EARTH-(Ireland)CLMcom-CLM-CCLM4-8-17(EU).csv") %>%
    filter(lubridate::month(date) %in% months_to_study) %>%
    dplyr::rename(temp = mintp) %>%
    dplyr::mutate(temp = -1*temp)
}else{
  clim_data = vroom("~/JRSS_organised_code/data/raw_data/clim/r12i1p1-ICHEC-EC-EARTH-(Ireland)CLMcom-CLM-CCLM4-8-17 (EU).csv") %>%
  filter(lubridate::month(date) %in% months_to_study) %>%
    dplyr::rename(temp = maxtp)
}


clim_data$Lat = clim_data$Lat %>% signif(4)
clim_data$Long = clim_data$Long %>% signif(4)

clim_thresh = clim_data %>%
  group_by(id) %>% 
  summarise(threshold_95 = quantile(temp, 0.95),
            threshold_96 = quantile(temp, 0.96),
            threshold_97 = quantile(temp, 0.97))

clim_thresh %>% write_csv(paste0("data/processed/clim_thresh_",var_modelled,".csv"))
clim_data = clim_data %>% left_join(clim_thresh)

# project coordinates climate data 
my_coords = clim_data %>% dplyr::select(Long, Lat)
sp::coordinates(my_coords)<- c("Long", "Lat")
sp::proj4string(my_coords) <- sp::CRS("+proj=longlat +datum=WGS84")
proj_cords <- sp::spTransform(my_coords,
                              sp::CRS(paste0("+proj=utm +zone=29 ellps=WGS84")))
proj_cords = proj_cords %>% sp::coordinates()/100000
clim_data$Long.projected = proj_cords[,1]
clim_data$Lat.projected = proj_cords[,2]
clim_data$id = clim_data %>% group_indices(Long, Lat)
clim_data %>% write_csv(paste0("data/processed/full_clim_data_",var_modelled,".csv"))
clim_data = read_csv(paste0("data/processed/full_clim_data_",var_modelled,".csv"))

# ========  ========  Tail Model  ========  ========  ======== 
print("model climate tail")


ngll = function(par){
  if(par <= 0) return(2^30)
  if(shape_param < 0){
    if(any(par > -1/shape_param)) return(2^30)
    if(any((1+shape_param*this.dat/par)< 0)) return(2^30)
  }
  -sum(evd::dgpd(x = this.dat, loc=0, scale = par, shape=shape_param, log=T))
}

estiamte_scale_fixed_shape = function(x,shape_c){
  this.dat <<- x
  shape_param <<- shape_c
  optim(par = c(0), fn = ngll, method = 'Brent', lower=0, upper = 5)
}

clim_data_extreme_96 = clim_data %>%
  group_by(id) %>%
  mutate(excess = temp - threshold_96) %>%
  filter(excess > 0) %>%
  ungroup()
clim_data_extreme_97 = clim_data %>%
  group_by(id) %>%
  mutate(excess = temp - threshold_97) %>%
  filter(excess > 0) %>%
  ungroup()
clim_data_extreme_95 = clim_data %>%
  group_by(id) %>%
  mutate(excess = temp - threshold_95) %>%
  filter(excess > 0) %>%
  ungroup()

num_sites = clim_data_extreme_95$id %>% unique()

est_1 = evd::fpot(clim_data_extreme_95$excess, threshold = 0,std.err = FALSE)$par[2]
est_2 = evd::fpot(clim_data_extreme_96$excess, threshold = 0,std.err = FALSE)$par[2]
est_3 = evd::fpot(clim_data_extreme_97$excess, threshold = 0,std.err = FALSE)$par[2]

potential_shape_values_climate = seq(min(est_1, est_2, est_3) - 0.1, max(est_1, est_2, est_3) + 0.1, length.out = 50)

if(refit_clim_shape_param){
  
  scales_95 = c()
  scales_96 = c()
  scales_97 = c()
  

  loglik_sum_95 = c()
  loglik_sum_96 = c()
  loglik_sum_97 = c()

  for(potential_shape in potential_shape_values_climate){
    print(potential_shape)
    
    loglik_95 = c()
    loglik_96 = c()
    loglik_97 = c()
    
    for(i in num_sites){
      

      # 0.9 thresh
      this_clim_95 = clim_data_extreme_95 %>% filter(id == i) %>% pull(excess)
      model_fit_95 = estiamte_scale_fixed_shape(this_clim_95, potential_shape)
      scales_95 = c(scales_95, model_fit_95$par)
      loglik_95 = c(loglik_95, model_fit_95$value)
      
      # 0.96 thresh
      this_clim_96 = clim_data_extreme_96 %>% filter(id == i) %>% pull(excess)
      model_fit_96 = estiamte_scale_fixed_shape(this_clim_96, potential_shape)
      scales_96 = c(scales_96, model_fit_96$par)
      loglik_96 = c(loglik_96, model_fit_96$value)

      # 0.97 thresh
      this_clim_97 = clim_data_extreme_97 %>% filter(id == i) %>% pull(excess)
      model_fit_97 = estiamte_scale_fixed_shape(this_clim_97, potential_shape)
      scales_97 = c(scales_97, model_fit_97$par)
      loglik_97 = c(loglik_97, model_fit_97$value)
      
    }
    
    loglik_sum_95 = c(loglik_sum_95, sum(loglik_95))
    loglik_sum_96 = c(loglik_sum_96, sum(loglik_96))
    loglik_sum_97 = c(loglik_sum_97, sum(loglik_97))
    
  }

  optimal_shape_95 = potential_shape_values_climate[which.min(loglik_sum_95)]
  optimal_shape_96 = potential_shape_values_climate[which.min(loglik_sum_96)]
  optimal_shape_97 = potential_shape_values_climate[which.min(loglik_sum_97)]
  
  optimal_shape_95 %>% saveRDS(paste0("output/optimal_shape_clim_95_", var_modelled))
  optimal_shape_96 %>% saveRDS(paste0("output/optimal_shape_clim_96_", var_modelled))
  optimal_shape_97 %>% saveRDS(paste0("output/optimal_shape_clim_97_", var_modelled))
  
}else{
  optimal_shape_95 = readRDS(paste0("output/optimal_shape_clim_95_", var_modelled))
  optimal_shape_96 = readRDS(paste0("output/optimal_shape_clim_96_", var_modelled))
  optimal_shape_97 = readRDS(paste0("output/optimal_shape_clim_97_", var_modelled))
}

scales_95 = c()
scales_96 = c()
scales_97 = c()

clim_ids = clim_data_extreme_95 %>%
  dplyr::select(Long, Lat, Long.projected, Lat.projected, id) %>%
  unique()

for(i in clim_ids$id){
  this_clim_95 = clim_data_extreme_95 %>% filter(id == i) %>% pull(excess)
  model_fit_95 = estiamte_scale_fixed_shape(this_clim_95, optimal_shape_95)
  scales_95 = c(scales_95, model_fit_95$par)

  this_clim_97 = clim_data_extreme_97 %>% filter(id == i) %>% pull(excess)
  model_fit_97 = estiamte_scale_fixed_shape(this_clim_97, optimal_shape_97)
  scales_97 = c(scales_97, model_fit_97$par)
  
  this_clim_96 = clim_data_extreme_96 %>% filter(id == i) %>% pull(excess)
  model_fit_96 = estiamte_scale_fixed_shape(this_clim_96, optimal_shape_96)
  scales_96 = c(scales_96, model_fit_96$par)
}

clim_ids$scales_95 = scales_95
clim_ids$scales_96 = scales_96
clim_ids$scales_97 = scales_97


clim_c_plt = rbind(clim_ids %>% dplyr::select(Long, Lat, sigma = scales_95) %>% mutate(model = "0.95"),
      clim_ids %>% dplyr::select(Long, Lat, sigma = scales_96) %>% mutate(model = "0.96"),
      clim_ids %>% dplyr::select(Long, Lat, sigma = scales_97) %>% mutate(model = "0.97")) %>%
  ggplot()+
  geom_point(aes(Long, Lat, col = sigma), size = 2.5)+
   geom_sf(data = ireland_sf, alpha = 0, col = 'black')+
   ggplot2::scale_color_gradientn(colors = my_pal)+
  scale_x_continuous(breaks = c(-10, -8, -6))+
   theme_minimal(12)+
  theme(strip.text = element_blank())+
   labs(col = expression(sigma[c]), x = '', y ='')+
  facet_wrap(~model, nrow = 1)

ggsave(clim_c_plt, filename = 'output/figs/clim_c_different_thresholds.pdf', height = 4, width = 9)


clim_thresh = read_csv(paste0("data/processed/clim_thresh_",var_modelled,".csv"))
clim_ids %>%
  left_join(clim_thresh) %>%
  write_csv(paste0("data/processed/winter_scales_clim_grid_",var_modelled,".csv"))
