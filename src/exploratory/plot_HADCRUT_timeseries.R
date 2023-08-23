# this script plots temperature anomolies on a global and irish scale 

gc()
rm(list = ls())
library(tidyverse)
library(vroom)
library(mgcv)
library(raster) # package for netcdf manipulation


# months to include in the analysis
months_to_study = c(12,1,2) # winter months we consider
refit_clim_shape_param = F

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

# map of ireland sf
ireland_sf <-rnaturalearth::ne_countries(type = "map_units",
                                         scale='large', # resolution of map
                                         returnclass = 'sf', # spatial object
                                         continent = "europe") %>%
  .[.$name %in% c('Ireland', 'N. Ireland'),]

# ----- HADCRUT
hadcrut_data = raster::brick("~/data/raw/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc")
hadcrut_loc_dat = raster::coordinates(hadcrut_data) %>% as_tibble()
names(hadcrut_loc_dat) = c('Long', 'Lat')
hadcrut_data = hadcrut_data %>% raster::values() %>% as_tibble()
hadcrut_data$Long = hadcrut_loc_dat$Long
hadcrut_data$Lat = hadcrut_loc_dat$Lat
names(hadcrut_data) = names(hadcrut_data) %>% str_remove_all("X")

hadcrut_data = hadcrut_data %>%
  tidyr::pivot_longer(c(-Long, -Lat), names_to = "date", values_to = "temp_anom") %>% # Making each row an observation
  mutate(date = date %>% lubridate::ymd())

irealnd_hadcrut = hadcrut_data %>%
  filter(Long > -11, Long < -6,   Lat > 50, Lat < 56) %>%
  mutate(year = lubridate::year(date)) %>%
  mutate(month = lubridate::month(date)) %>%
  filter(month %in% months_to_study) 

# store all data in one tibble for plotting
all_had = rbind(irealnd_hadcrut %>%
        drop_na() %>% 
        mutate(year = lubridate::year(date)) %>%
        filter(lubridate::month(date) %in% months_to_study) %>%
        group_by(year) %>%
        summarise(temp_anom = mean(temp_anom)) %>%
        unique() %>% mutate(lab = 'b'),
      hadcrut_data %>%
        drop_na() %>%
        filter(lubridate::month(date) %in% months_to_study) %>%
        mutate(year = lubridate::year(date)) %>%
        # filter(year>1930) %>%
        group_by(year) %>%
        summarise(temp_anom = mean(temp_anom)) %>%
        ungroup()%>%
        unique() %>% mutate(lab = 'a'))

plt = gridExtra::grid.arrange(all_had %>%
                          filter(lab == 'a')%>%
                          filter(year>1930) %>%
                          ggplot()+
                          geom_point(data = all_had %>%
                                       filter(lab == 'a') %>%
                                       filter(year == 2010), 
                                     aes(year, temp_anom), col = 'magenta', size = 2)+
                          geom_line(aes(year, temp_anom), size = 0.5, alpha = 0.7)+
                          theme_minimal(10)+
                          geom_hline(yintercept = 0, linetype = 'longdash', col = 'red')+
                          scale_x_continuous(breaks = c(1930, 1950, 1970, 1990, 2010))+
                          scale_y_continuous(limits = c(-1.7, 1.4),
                                             breaks = c( -1,  0, 1),
                                             label = paste0(c(-1,  0, 1),"°C"))+
                          labs(y = "Temperature\nanomaly      ", x = "Year")+
                          theme(strip.text = element_blank(),
                                axis.title.y = element_text(angle = 0, vjust = 0.5)),
                        all_had %>%
                          filter(lab == 'b')%>%
                          filter(year>1930) %>%
                          ggplot()+
                          geom_point(data = all_had %>%
                                       filter(lab == 'b') %>%
                                       filter(year == 2010), 
                                     aes(year, temp_anom), col = 'magenta', size = 2)+
                          geom_line(aes(year, temp_anom), size = 0.5, alpha = 0.7)+
                          theme_minimal(10)+
                          geom_hline(yintercept = 0, linetype = 'longdash', col = 'red')+
                          scale_x_continuous(breaks = c(1930, 1950, 1970, 1990, 2010))+
                          scale_y_continuous(limits = c(-1.7, 1.4),
                                             breaks = c( -1,  0, 1),
                                             label = paste0(c(-1,  0, 1),"°C"))+
                          labs(y = "", x = "Year")+
                          theme(strip.text = element_blank(),
                                axis.title.y = element_text(angle = 0, vjust = 0.5)), nrow = 1, widths = c(1, 0.85))
  
ggsave(plot = plt, filename = "output/figs/hadcrut_glob_v_irel.pdf", height = 2.5, width = 7)
