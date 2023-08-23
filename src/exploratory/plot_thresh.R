
# this script calculates and plots thresold using different methods at different quntiles levels
# specifically, with and without coastal proximity covairaita


gc()
rm(list = ls())
library(tidyverse)

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

ireland_sf <-rnaturalearth::ne_countries(type = "map_units",
                                         scale='large', # resolution of map
                                         returnclass = 'sf', # spatial object
                                         continent = "europe") %>%
  .[.$name %in% c('Ireland', 'N. Ireland'),]


clim_data = read_csv('~/data/processed/winter_scales_clim_grid_mintp.csv')

dist_sea = read_csv("~/data/processed/sites_clim_sea_dist.csv") %>%
  dplyr::select(Long, Lat, dist_sea) %>%
  mutate(dist_sea_logged = log(1+dist_sea)) %>%
  mutate(dist_sea_logged = dist_sea_logged - readRDS("output/dist_sea_scaling")) %>%
  dplyr::select(-dist_sea)

dist_sea$Lat = dist_sea$Lat %>% signif(4)
dist_sea$Long = dist_sea$Long %>% signif(4)


clim_data = clim_data %>%
  left_join(dist_sea)

clim_data$thresh_95_coast = predict(readRDS("output/threshold_model/threshold_model_tn_l_95_coast"), newdata = clim_data)$location
clim_data$thresh_96_coast = predict(readRDS("output/threshold_model/threshold_model_tn_l_96_coast"), newdata = clim_data)$location
clim_data$thresh_97_coast = predict(readRDS("output/threshold_model/threshold_model_tn_l_97_coast"), newdata = clim_data)$location

clim_data$thresh_95 = predict(readRDS("output/threshold_model/threshold_model_tn_l_95"), newdata = clim_data)$location
clim_data$thresh_96 = predict(readRDS("output/threshold_model/threshold_model_tn_l_96"), newdata = clim_data)$location
clim_data$thresh_97 = predict(readRDS("output/threshold_model/threshold_model_tn_l_97"), newdata = clim_data)$location


clim_data = clim_data %>%
  dplyr::select(Long, Lat, thresh_95, thresh_96, thresh_97, thresh_95_coast, thresh_96_coast, thresh_97_coast) %>%
  pivot_longer(-c(Long, Lat)) %>%
  rename(thresh_val = name) 

clim_data$thresh_val = clim_data$thresh_val %>% factor( levels = c('thresh_95', 'thresh_96', 'thresh_97', 'thresh_95_coast', 'thresh_96_coast', 'thresh_97_coast'))


clim_thresh = clim_data %>%
  ggplot()+
  geom_point(aes(Long, Lat, col = -value))+
  facet_wrap(~thresh_val, nrow = 2)+
  scale_x_continuous(breaks = c(-10, -8, -6))+
  geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
  ggplot2::scale_color_gradientn(colors = my_pal)+
  theme_minimal(12)+
  labs(x = "", y = "", col = 'u (°C)')+
  theme(strip.text = element_blank())


ggsave(clim_thresh, filename = 'output/figs/clim_thresholds.pdf', height = 4.8, width = 8)

