# this script pre process data for caluclate return levels of climate grid,
# over a series of very high quantiles, in order to accuratly estimates very high RLs and threshold exceedence properties
gc()
rm(list = ls())
library(tidyverse)
library(vroom)
library(evgam)
library(mgcv)
library(raster) # package for netcdf manipulation
# 

quantiles_to_estimate_bulk = seq(0.001,0.999,length.out = num_quantiles)

quantiles_to_estimate_bulk = c(0.995, 0.999, 0.9995, 0.9999)
quantiles_to_estimate_bulk = sort(quantiles_to_estimate_bulk)

quantiles_to_estimate_bulk %>% saveRDS("output/quantiles_to_estimate_bulk_clim")



obs_data = read_csv('data/processed/obs_data_mintp_winter.csv')

clim_data = read_csv("data/processed/full_clim_data_mintp.csv") %>% filter(lubridate::month(date) %in% unique(obs_data$month))

clim_quantiles = clim_data %>%
  group_by(id) %>%
  group_map(~{
    res = c()
    for(q in quantiles_to_estimate_bulk){
      res = rbind(res, tibble(id = .x$id[1],
                              Long = .x$Long[1],
                              Lat = .x$Lat[1],
                              Long.projected = .x$Long.projected[1],
                              Lat.projected = .x$Lat.projected[1],
                              quantile = q,
                              value = as.numeric(quantile(.x$temp, q))))
    }
    res
  }, .keep = T)%>%
  plyr::rbind.fill() %>%
  as.tibble()


clim_quantiles_subset = clim_quantiles %>%
  filter(id %in% obs_data$id) %>%
  group_by(id) %>%
  group_map(~{
    
    tibble(id = .x$id[1], quantile = list(.x$quantile), value = list(.x$value))
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()

clim_quantiles_subset %>%
  saveRDS(paste0("data/processed/clim_data_for_bulk_model_mintp_clim.csv"))

obs_data = obs_data %>%
  left_join(clim_quantiles_subset)

# add threshold to obs data
obs_data %>%
  saveRDS(paste0("data/processed/obs_data_for_bulk_model_mintp_clim.csv"))

