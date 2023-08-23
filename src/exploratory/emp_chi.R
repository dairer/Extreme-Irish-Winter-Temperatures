# calcualte 

rm(list = ls())
library(tidyverse)
m = "E2"

# quantile to estimate chi for
my_qnt = 0.9

# data set with three columns, first two columns contains unique pairs of sites, the thrid contians the respective distance between them.
site_pairs = readRDS("~/data/processed/obs_site_pairs")


# calculate extremal dependence between two vectors on uniform scale 
chi.emp=function(U,V,q ){
  sum((U>=q)&(V>=q))/sum((U>=q))
}
  
# --- Observational data
# itterate over bootstraps, caclulating chi 
for(bts_num in seq(200)){
  print(bts_num)
  this_dat = readRDS(paste0("output/bootstrapped_data_sets/model_","m70","_bootstrap_", bts_num))

  site_pairs%>%
    mutate(bin = ntile(distances, 25)) %>%
    group_by(bin)%>%
    group_map(~{
      
      print(.x$bin[1])
      x1_weighted = c()
      x2_weighted = c()
      
      for(s in seq(nrow(.x))){
        s1 = this_dat %>% filter(Station == .x[s,]$site_1)
        s2 = this_dat %>% filter(Station == .x[s,]$site_2)
        x1_weighted = c(x1_weighted,  s1 %>% filter(date %in% s2$date) %>% arrange(date) %>% pull(unif))
        x2_weighted = c(x2_weighted,  s2 %>% filter(date %in% s1$date) %>% arrange(date) %>% pull(unif))
      }
      
      tibble(bts = bts_num, bin = .x$bin[1], dist = mean(.x$distances), 
             chi = chi.emp(x1_weighted, x2_weighted, my_qnt)) %>%
        write_csv(paste0("output/new_chi_obs_bts_25_bins", my_qnt,".csv"), append = T)
      
      tibble()
    }, .keep = T)
}






# --- Climate data


clim_data = vroom("~/data/processed/mintp_r12i1p1-ICHEC-EC-EARTH-(Ireland)CLMcom-CLM-CCLM4-8-17(EU).csv") %>%
  dplyr::rename(temp = mintp) %>%
  dplyr::mutate(temp = -1*temp)

clim_data = clim_data %>%
  group_by(id)%>%
  group_map(~{
    dat = .x$temp
    this_dates = .x$date
    unif = ecdf(dat)(dat)
    tibble(id = .x$id[1], date = this_dates, unif = unif, Long = .x$Long[1], Lat = .x$Lat[1], Long.projected = .x$Long.projected[1],Lat.projected = .x$Lat.projected[1])
  }, .keep = T)%>%
  plyr::rbind.fill()%>%
  as_tibble()


site_pairs = readRDS("~/data/processed/empirical_chi/clim_sites_and_distances")

chi.emp=function(U,V,q ){
  sum((U>=q)&(V>=q))/sum((U>=q))
}

site_pairs %>%
  mutate(bin = ntile(distances, 25)) %>%
  group_by(bin) %>%
  group_map(~{
    
    for(i in seq(100)){
      print(i)
      this_sample = .x %>% sample_n(size = nrow(.x), replace = T)
      
      tibble(bts = i, bin = .x$bin[1], dist = mean(.x$distances), 
             chi = chi.emp(clim_data[lapply(as.numeric(this_sample$site_1), function(x) which(clim_data$id %in% x)) %>% unlist(),]$unif,
                           clim_data[lapply(as.numeric(this_sample$site_2), function(x) which(clim_data$id %in% x)) %>% unlist(),]$unif, my_qnt)) %>%
        write_csv(paste0("output/chi_clim_bts_25_bins", my_qnt,".csv"), append = T)
    }    
  }, .keep = T)




# plot obs and clim chi
plt = rbind(read_csv(paste0("output/new_chi_obs_bts_25_bins", 0.8,".csv"),
               col_names = c("i", "bts", "dist", "chi")) %>% mutate(d = "obs",p = "p = 0.8"),
      read_csv(paste0("output/new_chi_obs_bts_25_bins", 0.85,".csv"),
               col_names = c("i", "bts", "dist", "chi")) %>% mutate(d = "obs",p = "p = 0.85"),
      read_csv(paste0("output/new_chi_obs_bts_25_bins", 0.9,".csv"),
               col_names = c("i", "bts", "dist", "chi")) %>% mutate(d = "obs",p = "p = 0.9"),
      read_csv(paste0("output/chi_clim_bts_25_bins", 0.8,".csv"),
               col_names = c("i", "bts", "dist", "chi")) %>% mutate(d = "clim", p = "p = 0.8"),
      read_csv(paste0("output/chi_clim_bts_25_bins", 0.85,".csv"),
               col_names = c("i", "bts", "dist", "chi")) %>% mutate(d = "clim",p = "p = 0.85"),
      read_csv(paste0("output/chi_clim_bts_25_bins", 0.9,".csv"),
               col_names = c("i", "bts", "dist", "chi")) %>% mutate(d = "clim",p = "p = 0.9")
      )%>%
  group_by(dist, p, d) %>%
  summarise(upper = quantile(chi, 0.025, na.rm = T),
            lower = quantile(chi, 0.975, na.rm = T),
            chi = mean(chi, na.rm = T)) %>%
  ggplot()+
  geom_segment(aes(x = dist*100, xend = dist*100, y = lower, yend = upper, group = d, col = d))+
  geom_point(aes(x = dist*100, y=chi, shape = d, group = d, col = d))+
  facet_wrap(~p, nrow = 1)+
  theme_minimal(12)+
  labs(y = expression(chi),
       x = "Distance (km)")+
  theme(legend.position = 'none',
        axis.title.y = element_text(angle=0, vjust=.5))+
  ylim(0,1)

ggsave(filename = 'output/figs/chi_emp.pdf', width = 7, height = 2.5)
