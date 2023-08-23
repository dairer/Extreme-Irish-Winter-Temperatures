# Generate plots of covaraites

gc()
rm(list = ls())
library(tidyverse)
library(vroom)
library(mgcv)
library(raster) # package for netcdf manipulation

months_to_study = c(12,1,2) # winter months we consider

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

obs_data = read_csv('data/processed/obs_all_data_mintp.csv') %>%
  dplyr::select(-loc_id) %>%
  mutate(year = lubridate::year(date)) %>%
  filter(mintp >-19.2)  %>%
  filter(lubridate::month(date) %in% months_to_study) %>%
  rename(temp = mintp) %>%
  mutate(temp = -1*temp) %>%
  left_join(read_csv("~/data/processed/obs_data_dist_to_sea.csv")) %>%
  drop_na() 

# ----- HADCRUT data
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

ireland_sf <-rnaturalearth::ne_countries(type = "map_units",
                                         scale='large', # resolution of map
                                         returnclass = 'sf', # spatial object
                                         continent = "europe") 
europe = sf::st_union(ireland_sf)


# plot time series of residuals of temeprature anomolies from HADCRUT smoothing on a monthly and (averaged) yearly
irealnd_hadcrut = hadcrut_data %>%
  filter(Long > -11, Long < -6,   Lat > 50, Lat < 56) %>%
  mutate(year = lubridate::year(date)) %>%
  mutate(month = lubridate::month(date)) %>%
  filter(month %in% months_to_study) 

loess_fit_irel = predict(loess(temp_anom ~ year, data = irealnd_hadcrut), se=T)
irealnd_hadcrut$loess_temp_anom = loess_fit_irel$fit
irealnd_hadcrut$residuals = irealnd_hadcrut$loess_temp_anom - irealnd_hadcrut$temp_anom

obs_data = obs_data %>%
  mutate(month = lubridate::month(date)) %>%
  left_join(irealnd_hadcrut %>%
              dplyr::select(year, month,loess_temp_anom, temp_anom, residuals) %>%
              unique())

irealnd_hadcrut_yrly = irealnd_hadcrut %>%
  group_by(year) %>%
  summarise(yearly_residuals = mean(residuals)) 

plt = gridExtra::grid.arrange(irealnd_hadcrut %>%
                          ggplot() + 
                          theme_minimal()+
                          geom_hline(yintercept = 0, linetype = 'longdash', col = 'red')+
                          geom_line(aes(date, residuals), size = 0.5, alpha = 0.7)+
                          theme_minimal(12)+
                          ylim(-1.8, 2.8)+
                          xlim(lubridate::ymd("1950-01-01"), lubridate::ymd("2022-12-31"))+
                          labs( y = expression(M[paste('r,m')]^I))+
                          theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
                                axis.title.x = element_blank()),
                        irealnd_hadcrut_yrly %>%
                          ggplot() + 
                          theme_minimal()+
                          geom_hline(yintercept = 0, linetype = 'longdash', col = 'red')+
                          geom_line(aes(year, yearly_residuals), size = 0.5, alpha = 0.7)+
                          theme_minimal(12)+
                          ylim(-1.8, 2.8)+
                          xlim(1950,2022)+
                          labs( y = expression(M[r]^I))+
                          theme(axis.title.y = element_text(angle = 0, vjust = 0.5), 
                                axis.title.x = element_blank()),
                        nrow = 1, bottom ='Year')

ggsave(plot = plt, filename = 'output/figs/MIR.pdf', height = 2.5, width = 8)

# plot HADCRUT grid along with time series of residuals of temeprature anomolies from HADCRUT smoothing 
# from grid points surrounding Ireland and to the north of Ireland
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
temp_anom_w_loess$temp_anom_around_irel_resids = (predict(mod) - temp_anom_w_loess$temp_anom_around_irel)

mod = lm(temp_anom_north_irel~loess_temp_anom, data = temp_anom_w_loess)
temp_anom_w_loess$temp_anom_north_irel_resids = (predict(mod) - temp_anom_w_loess$temp_anom_north_irel)

ireland_sf <-rnaturalearth::ne_countries(type = "map_units",
                                         scale='small', # resolution of map
                                         returnclass = 'sf', # spatial object
                                         continent = "europe") 

mg_plt = gridExtra::grid.arrange(hadcrut_data %>%
                                   filter(date == "2020-02-15") %>%
                                   drop_na() %>%
                                   filter(lubridate::month(date) %in% months_to_study) %>%
                                   filter(Long != -7.5 | Lat!= 52.5)%>%
                                   filter(Lat > 45, Long > -15, Long < 0, Lat<60) %>%
                                   ggplot()+
                                   geom_tile(aes(Long, Lat), fill = 'forestgreen', col = 'black', alpha = 0.3, size = 0.25)+
                                   geom_polygon(data = positions <- data.frame(
                                     x = c(-10, -10, -5, -5, -10),
                                     y = c(55, 60, 60, 55, 55)), aes(x,y), alpha = 0.8, col = 'black', fill = 'forestgreen')+
                                   viridis::scale_fill_viridis()+
                                   geom_sf(data = europe, alpha = 0, col = 'black', fill = 'black')+
                                   ylim(42, 63)+
                                   scale_x_continuous(limits = c(-18, 3), breaks=c(-15, -10, -5, 0))+
                                   geom_line(data = data.frame(x = c(-15, -15), y = c(42, 63)), aes(x, y), col = 'black')+
                                   geom_line(data = data.frame(x = c(-10, -10), y = c(42, 63)), aes(x, y), col = 'black')+
                                   geom_line(data = data.frame(x = c(-5, -5), y = c(42, 63)), aes(x, y), col = 'black')+
                                   geom_line(data = data.frame(x = c(-0, -0), y = c(42, 63)), aes(x, y), col = 'black')+
                                   geom_line(data = data.frame(x = c(-18, 3), y = c(45,45)), aes(x, y), col = 'black')+
                                   geom_line(data = data.frame(x = c(-18, 3), y = c(50,50)), aes(x, y), col = 'black')+
                                   geom_line(data = data.frame(x = c(-18, 3), y = c(55,55)), aes(x, y), col = 'black')+
                                   geom_line(data = data.frame(x = c(-18, 3), y = c(60,60)), aes(x, y), col = 'black')+
                                   theme_minimal(12)+
                                   theme(legend.position = 'none')+
                                   theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank()),
                          temp_anom_w_loess %>%
                            ggplot()+
                            geom_hline(yintercept = 0, linetype = 'longdash', col = 'red')+
                            geom_line(aes(year, temp_anom_around_irel_resids), size = 0.5)+
                            theme_minimal(12)+
                            xlim(1950, 2022)+
                            ylim(c(-0.8,1.2))+
                            labs(x = "Year", y = expression(M^G))+
                            theme(axis.title.y = element_text(angle = 0, vjust = 0.5)),
                          temp_anom_w_loess %>%
                            ggplot()+
                            geom_hline(yintercept = 0, linetype = 'longdash', col = 'red')+
                            geom_line(aes(year, temp_anom_north_irel_resids), size = 0.5)+
                            theme_minimal(12)+
                            xlim(1950, 2022)+
                            ylim(c(-0.8,1.2))+
                            labs(x = "Year", y = expression(M^N))+
                            theme(axis.title.y = element_text(angle = 0, vjust = 0.5)),
                        widths = c(1, 2, 2))

ggsave(plot = mg_plt, filename = "output/figs/mg.pdf", height = 3, width = 10)




# Plot NAO and AO
obs_data = obs_data %>%
  left_join(temp_anom_w_loess%>%
              dplyr::select(year, temp_anom_around_irel_resids, temp_anom_north_irel_resids)) %>%
  mutate(week = lubridate::week(date))

# - add AO and NAO covariates
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



loess_temp_anom_and_ao$ao = loess_temp_anom_and_ao$ao - mean(loess_temp_anom_and_ao$ao)
loess_temp_anom_and_nao$nao = loess_temp_anom_and_nao$nao - mean(loess_temp_anom_and_nao$nao)


plt = gridExtra::grid.arrange(loess_temp_anom_and_nao %>%
                          ggplot() + 
                          theme_minimal()+
                          geom_hline(yintercept = 0, linetype = 'longdash', col = 'red')+
                          geom_line(aes(year, nao), size = 0.5, alpha = 0.7)+
                          theme_minimal(12)+
                          xlim(1950, 2022)+
                          labs(x = "Year", y = "NAO")+
                          theme(axis.title.y = element_text(angle = 0, vjust = 0.5))+
                          ylim(-3,2.2),
                        loess_temp_anom_and_ao %>%
                          ggplot() + 
                          theme_minimal()+
                          geom_hline(yintercept = 0, linetype = 'longdash', col = 'red')+
                          geom_line(aes(year, ao), size = 0.5, alpha = 0.7)+
                          theme_minimal(12)+
                          xlim(1950, 2022)+
                          ylim(-3,2.2)+
                          labs(x = "Year", y = "AO")+
                          theme(axis.title.y = element_text(angle = 0, vjust = 0.5)), nrow = 1)

ggsave(plot = plt, filename = 'output/figs/nao_ao.pdf', height = 3, width =10)