
# Description: This script prepares observational temperature data for marginal 
# modelling. This involves incorporating all covariates as descriped in the 
# paper associated to this work.
# Note: bulk modelling in preformed in this script. 
gc()
rm(list = ls())
library(tidyverse)
library(vroom)
library(mgcv)
library(raster) # package for netcdf manipulation


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

clim_data = read_csv('data/processed/winter_scales_clim_grid_mintp.csv')

obs_data = read_csv('data/processed/obs_all_data_mintp.csv') %>%
                   dplyr::select(-loc_id) %>%
  mutate(year = lubridate::year(date)) %>%
  # filter(mintp >-19.2)  %>%
  filter(lubridate::month(date) %in% months_to_study) %>%
  # mutate(mintp = -1*mintp) %>%
  rename(temp = mintp) %>%
  mutate(temp = -1*temp) %>%
  left_join(read_csv("~/data/processed/obs_data_dist_to_sea.csv")) %>%
  # left_join(read_csv("~/chapter_6/data/processed/mean_wd.csv")) %>%
  drop_na() 

# remove anomolous temepratures 
anoms = obs_data %>%
  group_by(date) %>%
  group_map(~{
    x = .x$temp
    window = 4*sd(x)
    avg = mean(x)
    
    .x[(x > (avg+window) |  x < (avg-window)),]
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()%>%
  drop_na()

obs_data = obs_data %>% anti_join(anoms)


obs_data$dist_sea_logged = log(1+obs_data$dist_sea)
mean(obs_data$dist_sea_logged) %>% saveRDS("output/dist_sea_scaling")
obs_data$dist_sea_logged = obs_data$dist_sea_logged - mean(obs_data$dist_sea_logged)

obs_data %>%
  group_by(Station, Long, Lat) %>%
  summarise(n = n()/121) %>%
  ggplot()+
  geom_point(aes(Long, Lat, col = n, size = n))+
  geom_sf(data = ireland_sf, col = 'black', alpha = 0)+
  theme_minimal()+
  viridis::scale_color_viridis()


# ----- HADCRUT covaraites
# 1c. Add temporal covariates (HADCRUT5)
hadcrut_data = raster::brick("~/max_min_temp_dependence/data/raw/HadCRUT.5.0.1.0.analysis.anomalies.ensemble_mean.nc")
hadcrut_loc_dat = raster::coordinates(hadcrut_data) %>% as_tibble()
names(hadcrut_loc_dat) = c('Long', 'Lat')
hadcrut_data = hadcrut_data %>% raster::values() %>% as_tibble()
hadcrut_data$Long = hadcrut_loc_dat$Long
hadcrut_data$Lat = hadcrut_loc_dat$Lat
names(hadcrut_data) = names(hadcrut_data) %>% str_remove_all("X")

hadcrut_data = hadcrut_data %>%
  tidyr::pivot_longer(c(-Long, -Lat), names_to = "date", values_to = "temp_anom") %>% # Making each row an observation
  mutate(date = date %>% lubridate::ymd())


ireland_sf <-rnaturalearth::ne_countries(type = "map_units",
                                         scale='large', # resolution of map
                                         returnclass = 'sf', # spatial object
                                         continent = "europe") 

europe = sf::st_union(ireland_sf)

irealnd_hadcrut = hadcrut_data %>%
  filter(Long > -11, Long < -6,   Lat > 50, Lat < 56) %>%
  mutate(year = lubridate::year(date)) %>%
  mutate(month = lubridate::month(date)) %>%
  filter(month %in% months_to_study) 

loess_fit_irel = predict(loess(temp_anom ~ year, data = irealnd_hadcrut), se=T)
irealnd_hadcrut$loess_temp_anom = loess_fit_irel$fit
irealnd_hadcrut$loess_temp_anom_u = loess_fit_irel$fit+1.96*(loess_fit_irel$se.fit)
irealnd_hadcrut$loess_temp_anom_l = loess_fit_irel$fit-1.96*(loess_fit_irel$se.fit)
irealnd_hadcrut$residuals = irealnd_hadcrut$loess_temp_anom - irealnd_hadcrut$temp_anom

irealnd_hadcrut_yrly = irealnd_hadcrut %>%
  group_by(year) %>%
  summarise(yearly_residuals = mean(residuals)) 


plt = irealnd_hadcrut %>%
  filter(year>=1950) %>%
  ggplot()+
  geom_ribbon(aes(x = year, ymin = loess_temp_anom_l, ymax = loess_temp_anom_u), alpha = 0.3, fill = 'forestgreen')+
  geom_line(aes(year, loess_temp_anom))+
  theme_minimal()+
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1))+
  labs(y = "Temmperature\n anomaly(°C)   ", x = "Year")+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

ggsave(plt, filename = "output/figs/winter_temp_change.pdf", height = 3, width = 5)
obs_data = obs_data %>%
  mutate(month = lubridate::month(date)) %>%
  left_join(irealnd_hadcrut %>%
               dplyr::select(year, month,loess_temp_anom, temp_anom, residuals) %>%
               unique()) %>%
  left_join(irealnd_hadcrut_yrly)


# --- get hadcrut around ireland
around_irel = hadcrut_data %>%
  drop_na() %>%
  filter(lubridate::month(date) %in% months_to_study) %>%
  filter(Long != -7.5 | Lat!= 52.5) %>%
  filter(Lat > 45, Long > -15, Long < 0, Lat<60) %>%
  mutate(year = lubridate::year(date)) %>%
  filter(year>1930) %>%
  group_by(year) %>%
  summarise(temp_anom_around_irel = mean(temp_anom))


north_of_irel = hadcrut_data %>%
  drop_na() %>%
  filter(lubridate::month(date) %in% months_to_study) %>%
  filter(Long != -7.5 | Lat!= 52.5) %>%
  filter(Lat > 55, Long > -15, Long < 0, Lat<60) %>%
  mutate(year = lubridate::year(date)) %>%
  filter(year>1930) %>%
  group_by(year) %>%
  summarise(temp_anom_north_irel = mean(temp_anom))


temp_anom_w_loess = obs_data %>%
  left_join(around_irel) %>%
  left_join(north_of_irel) %>%
  dplyr::select(year, temp_anom_around_irel, temp_anom_north_irel, loess_temp_anom) %>%
  unique() %>%
  drop_na

mod = lm(temp_anom_around_irel~loess_temp_anom, data = temp_anom_w_loess)
temp_anom_w_loess$temp_anom_around_irel_resids = (temp_anom_w_loess$temp_anom_around_irel - predict(mod))

mod = lm(temp_anom_north_irel~loess_temp_anom, data = temp_anom_w_loess)
temp_anom_w_loess$temp_anom_north_irel_resids = (temp_anom_w_loess$temp_anom_north_irel - predict(mod))


plot(temp_anom_w_loess$temp_anom_north_irel_resids)
irealnd_hadcrut_yrly$yearly_residuals


irealnd_hadcrut_yrly %>%
  left_join(temp_anom_w_loess %>%
              dplyr::select(year, temp_anom_around_irel_resids, temp_anom_north_irel_resids)) %>%
  drop_na() %>%
  ggplot()+
  geom_point(aes(yearly_residuals, temp_anom_around_irel_resids))

temp_anom_w_loess$year %>% unique() %>% sort
hadcrut_data %>%
  drop_na() %>%
  filter(lubridate::month(date) %in% months_to_study) %>%
  filter(Long != -7.5 | Lat!= 52.5) %>%
  filter(Lat > 45, Long > -15, Long < 0, Lat<60) %>%
  mutate(year = lubridate::year(date)) %>%
  filter(year>1930) %>%
  group_by(year) %>%
  summarise(temp_anom = mean(temp_anom)) %>%
  ungroup() %>%
  ggplot()+
  geom_point(aes(year, temp_anom))+
  geom_vline(xintercept = 2010)+
  theme_minimal()

obs_data = obs_data %>%
  left_join(temp_anom_w_loess%>%
              dplyr::select(year, temp_anom_around_irel_resids, temp_anom_north_irel_resids)) %>%
  mutate(week = lubridate::week(date))


# ----- add AO and NAO covariates
ao_data = read_table('https://www.cpc.ncep.noaa.gov/products/precip/CWlink/daily_ao_index/monthly.ao.index.b50.current.ascii',
                     col_names = c('year', 'month', 'ao')) %>%
  filter(month %in% months_to_study) %>%
  rename(ao_monthly = ao)

ao_data_yrly = ao_data %>%
  group_by(year) %>%
  summarise(ao = mean(ao_monthly)) %>% ungroup()

nao_data = read_csv("data/raw/nao.csv") %>%
  dplyr::select(year, Jan, Feb, Dec) %>%
  pivot_longer(-year)

nao_data$name[nao_data$name == "Jan"] = 1
nao_data$name[nao_data$name == "Feb"] = 2
nao_data$name[nao_data$name == "Dec"] = 12

nao_data$name = as.numeric(nao_data$name)

nao_data = nao_data %>%
  rename(nao_monthly = value,
         month = name)

nao_data_yrly = nao_data %>%
  group_by(year) %>%
  summarise(nao = mean(nao_monthly)) 

obs_data = obs_data %>%
  left_join(nao_data)%>%
  left_join(nao_data_yrly)%>%
  left_join(ao_data)%>%
  left_join(ao_data_yrly) %>% 
  drop_na()

# ---- get AO and NAO residuals 
loess_temp_anom_and_ao = obs_data %>%
  dplyr::select(year, loess_temp_anom, ao) %>%
  unique()%>%
  drop_na()

mod = lm(ao~loess_temp_anom, data = loess_temp_anom_and_ao)
loess_temp_anom_and_ao$ao_resid = (loess_temp_anom_and_ao$ao - predict(mod))

loess_temp_anom_and_ao_monthly = obs_data %>%
  dplyr::select(year, loess_temp_anom, ao_monthly) %>%
  unique() %>%
  drop_na()

mod = lm(ao_monthly~loess_temp_anom, data = loess_temp_anom_and_ao_monthly)
loess_temp_anom_and_ao_monthly$ao_monthly_resid = (loess_temp_anom_and_ao_monthly$ao_monthly - predict(mod))

loess_temp_anom_and_nao = obs_data %>%
  dplyr::select(year, loess_temp_anom, nao) %>%
  unique()%>%
  drop_na()

mod = lm(nao~loess_temp_anom, data = loess_temp_anom_and_nao)
loess_temp_anom_and_nao$nao_resid = (loess_temp_anom_and_nao$nao - predict(mod))

loess_temp_anom_and_nao_monthly = obs_data %>%
  dplyr::select(year, loess_temp_anom, nao_monthly) %>%
  unique()%>%
  drop_na()

mod = lm(nao_monthly~loess_temp_anom, data = loess_temp_anom_and_nao_monthly)
loess_temp_anom_and_nao_monthly$nao_monthly_resid = (loess_temp_anom_and_nao_monthly$nao_monthly - predict(mod))

obs_data = obs_data %>%
  left_join(loess_temp_anom_and_ao) %>%
  left_join(loess_temp_anom_and_ao_monthly) %>%
  left_join(loess_temp_anom_and_nao) %>%
  left_join(loess_temp_anom_and_nao_monthly)



# get climate subset that matches up with observational stations 
#  ------- calculate the distance between two points
dist_h <- function(long1, lat1, long2, lat2) {
  sqrt((long1 - long2)^2 + (lat1 - lat2)^2)
}

obs_sites = obs_data %>%
  dplyr::select(Station, Long, Lat, Long.projected, Lat.projected) %>%
  unique()

# itterate through all stations and find climate grid point
closest_irel = c()
for(i in seq(nrow(obs_sites))){
  smallest_dist = 9999999999
  id_of_smallest_dist_irel = 9999999999
  
  for(j in seq(nrow(clim_data))){
    this_dist = dist_h(obs_sites[i,]$Long.projected,
                       obs_sites[i,]$Lat.projected,
                       clim_data[j,]$Long.projected,
                       clim_data[j,]$Lat.projected)
    
    if(this_dist < smallest_dist){
      smallest_dist = this_dist
      id_of_smallest_dist_irel = clim_data$id[j]
    }
  }
  closest_irel = c(closest_irel, id_of_smallest_dist_irel)
}

# add column "id" to obs data which is the id of the closest clim grid point
obs_sites$id = closest_irel

obs_data = obs_data %>%
  left_join(obs_sites)

obs_data = obs_data %>%
  left_join(clim_data %>% dplyr::select(id,
                                        scales_95, scales_96, scales_97,
                                        threshold_95, threshold_96, threshold_97))


# ----- Fit threshold models
# ---- threshold = 0.95
qunatile_model_95 <-  temp ~ threshold_95
qunatile_model_fit_95 <- evgam::evgam(qunatile_model_95, obs_data, family = "ald", ald.args = list(tau = 0.95))
obs_data$threshold_tn_l_o_95 = qunatile_model_fit_95$location$fitted
qunatile_model_fit_95 %>% saveRDS("output/threshold_model/threshold_model_tn_l_95")

qunatile_model_95_w_coast <-  temp ~ threshold_95 + dist_sea_logged
qunatile_model_fit_95_w_coast <- evgam::evgam(qunatile_model_95_w_coast, obs_data, family = "ald", ald.args = list(tau = 0.95))
qunatile_model_fit_95_w_coast %>% saveRDS("output/threshold_model/threshold_model_tn_l_95_coast")
obs_data$threshold_tn_l_o_95_w_coast = qunatile_model_fit_95_w_coast$location$fitted

# ---- threshold = 0.96
qunatile_model_96 <-  temp ~ threshold_96
qunatile_model_fit_96 <- evgam::evgam(qunatile_model_96, obs_data, family = "ald", ald.args = list(tau = 0.96))
obs_data$threshold_tn_l_o_96 = qunatile_model_fit_96$location$fitted
qunatile_model_fit_96 %>% saveRDS("output/threshold_model/threshold_model_tn_l_96")

qunatile_model_96_w_coast <-  temp ~ threshold_96 + dist_sea_logged
qunatile_model_fit_96_w_coast <- evgam::evgam(qunatile_model_96_w_coast, obs_data, family = "ald", ald.args = list(tau = 0.96))
qunatile_model_fit_96_w_coast %>% saveRDS("output/threshold_model/threshold_model_tn_l_96_coast")
obs_data$threshold_tn_l_o_96_w_coast = qunatile_model_fit_96_w_coast$location$fitted

# ---- threshold = 0.97
qunatile_model_97 <-  temp ~ threshold_97
qunatile_model_fit_97 <- evgam::evgam(qunatile_model_97, obs_data, family = "ald", ald.args = list(tau = 0.97))
obs_data$threshold_tn_l_o_97 = qunatile_model_fit_97$location$fitted
qunatile_model_fit_97 %>% saveRDS("output/threshold_model/threshold_model_tn_l_97")

qunatile_model_97_w_coast <-  temp ~ threshold_97 + dist_sea_logged
qunatile_model_fit_97_w_coast <- evgam::evgam(qunatile_model_97_w_coast, obs_data, family = "ald", ald.args = list(tau = 0.97))
qunatile_model_fit_97_w_coast %>% saveRDS("output/threshold_model/threshold_model_tn_l_97_coast")
obs_data$threshold_tn_l_o_97_w_coast = qunatile_model_fit_97_w_coast$location$fitted

plt = gridExtra::grid.arrange(obs_data %>%
                          group_by(Station) %>%
                          summarise(threshold_tn_l_o_96,
                                    emp = quantile(temp, 0.96)) %>%
                          unique() %>%
                          ggplot()+
                          geom_point(aes(emp, threshold_tn_l_o_96))+
                          geom_abline()+
                          theme_minimal(12)+
                            xlim(-0.1, 7)+
                            ylim(-0.1, 7)+
                            
                          theme(axis.title = element_blank()),
                          obs_data %>%
                            group_by(Station) %>%
                            summarise(threshold_tn_l_o_96_w_coast,
                                      emp = quantile(temp, 0.96)) %>%
                            unique() %>%
                            ggplot()+
                            geom_point(aes(emp, threshold_tn_l_o_96_w_coast))+
                            geom_abline()+
                            theme_minimal(12)+
                            xlim(-0.1, 7)+
                            ylim(-0.1, 7)+
                            
                            theme(axis.title = element_blank()),
                          nrow = 1,
                        bottom = "Empiricle quantile",
                        left = "Estimated quantile")

ggsave(plt, filename = "output/figs/emp_vs_pred_qnt.pdf", width = 8, height = 3)



obs_data %>%
  dplyr::select(-c(date,temp)) %>%
  unique() %>% 
  write_csv("data/processed/obs_sites_and_closest_clim_mintp.csv")


# ---- centering covariates
obs_data$loess_temp_anom = obs_data$loess_temp_anom - mean(obs_data$loess_temp_anom)
obs_data$residuals = obs_data$residuals - mean(obs_data$residuals)
obs_data$yearly_residuals = obs_data$yearly_residuals - mean(obs_data$yearly_residuals)

obs_data$scales_95_logged = log(obs_data$scales_95)
mean(obs_data$scales_95_logged)%>% saveRDS("output/scales_95_scaling")
obs_data$scales_95_logged = obs_data$scales_95_logged - mean(obs_data$scales_95_logged)

obs_data$scales_96_logged = log(obs_data$scales_96)
mean(obs_data$scales_96_logged)%>% saveRDS("output/scales_96_scaling")
obs_data$scales_96_logged = obs_data$scales_96_logged - mean(obs_data$scales_96_logged)

obs_data$scales_97_logged = log(obs_data$scales_97)
mean(obs_data$scales_97_logged)%>% saveRDS("output/scales_97_scaling")
obs_data$scales_97_logged = obs_data$scales_97_logged - mean(obs_data$scales_97_logged)

obs_data$ao = obs_data$ao - mean(obs_data$ao)
obs_data$ao_monthly = obs_data$ao_monthly - mean(obs_data$ao_monthly)
obs_data$ao_resid = obs_data$ao_resid - mean(obs_data$ao_resid)
obs_data$nao = obs_data$nao - mean(obs_data$nao)
obs_data$nao_monthly = obs_data$nao_monthly - mean(obs_data$nao_monthly)
obs_data$temp_anom_around_irel_resids = obs_data$temp_anom_around_irel_resids - mean(obs_data$temp_anom_around_irel_resids)
obs_data$nao_resid = obs_data$nao_resid - mean(obs_data$nao_resid)
obs_data$nao_monthly_resid = obs_data$nao_monthly_resid - mean(obs_data$nao_monthly_resid)
obs_data$ao_monthly_resid = obs_data$ao_monthly_resid - mean(obs_data$ao_monthly_resid)

# remove any sites with no thresh exceedences
sites_to_rem = obs_data %>%
  group_by(Station) %>%
  group_map(~{
    
    
    tibble(Station = .x$Station[1], 
           num_extreme = .x %>%
             filter(temp > threshold_tn_l_o_97_w_coast) %>%
             nrow())
  }, .keep = T) %>%
  plyr::rbind.fill() %>% 
  as_tibble() %>%
  arrange((num_extreme)) %>%
  filter(num_extreme < 5) %>%
  pull(Station)

obs_data %>%
  filter(!(Station %in% sites_to_rem)) %>%
  write_csv("data/processed/obs_data_mintp_winter.csv")