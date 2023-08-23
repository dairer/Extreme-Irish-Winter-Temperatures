
# ---- this sccript plots the number of observations per month for winter minima

obs_data = read_csv('data/processed/obs_all_data_mintp.csv') %>%
  dplyr::select(-loc_id) %>%
  mutate(year = lubridate::year(date)) %>%
  mutate(mintp = -1*mintp) 

extreme_obs_95 = obs_data %>% 
  group_by(Station) %>%
  mutate(thresh = quantile(mintp, 0.95)) %>%
  ungroup %>%
  filter(mintp>thresh) %>%
  group_by(Month = as.character(lubridate::month(date)))

extreme_obs_95$Month = factor(extreme_obs_95$Month, levels = as.character(c(8,9,10,11,12,1,2,3,4,5,6,7)))
num_extreme_obs_95 = nrow(extreme_obs_95)

extreme_obs_99 = obs_data %>% 
  group_by(Station) %>%
  mutate(thresh = quantile(mintp, 0.99)) %>%
  ungroup %>%
  filter(mintp>thresh) %>%
  group_by(Month = as.character(lubridate::month(date)))
num_extreme_obs_99 = nrow(extreme_obs_99)


extreme_obs_99 = extreme_obs_99 %>%
  summarise(n = n())

extreme_obs_99 = rbind(extreme_obs_99, c(8,0), c(9,0)) # include months with no data in plot
extreme_obs_99$Month = factor(extreme_obs_99$Month, levels = as.character(c(8,9,10,11,12,1,2,3,4,5,6,7)))


plt = gridExtra::grid.arrange(extreme_obs_95 %>%
                          summarise(n = n()) %>%
                          ggplot()+
                          geom_bar(aes(Month, n/num_extreme_obs_95), stat = 'identity')+
                          theme_minimal()+
                          labs(y = "Proportion of extreme events")+
                          ylim(0, 0.33),
                        extreme_obs_99 %>%
                          ggplot()+
                          geom_bar(aes(Month, n/num_extreme_obs_99), stat = 'identity')+
                          theme_minimal()+
                          labs(y = "Proportion of extreme events")+
                          ylim(0, 0.33), nrow = 1)
ggsave(plt, filename = "output/figs/prop_obs_per_month_mintp_winter.pdf", width = 6, height = 2.75)




