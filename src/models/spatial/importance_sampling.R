# get importance sample estimates of spatial RL of high levels 

gc()
rm(list = ls())
library(tidyverse)


run_imp_samp = function(marg_mods, var_qnt, thresh_qnt){
  
  
  for(marg_mod in marg_mods){
    source('src/models/marginal/new_gpd_models.R')
    
    spatial_covars = read_csv("data/processed/obs_data_mintp_winter.csv") %>%
      dplyr::select(Station, 
                    scales_96_logged, threshold_tn_l_o_96_w_coast, 
                    scales_95_logged, threshold_tn_l_o_95_w_coast,
                    scales_97_logged, threshold_tn_l_o_97_w_coast,
                    Long.projected, Lat.projected, dist_sea_logged) %>%
      unique() 
    
    obs_grid = c()
    for(i in seq(length(spatial_covars$Station))){
      print(i)
      
      obs_grid = rbind(obs_grid, 
                       tibble(Station = spatial_covars$Station[i], 
                              year = (seq(1950, 2022))))
    }
    
    temporal_covars = read_csv("data/processed/obs_data_mintp_winter.csv") %>%
      dplyr::select(year, 
                    # month,
                    loess_temp_anom,
                    residuals
                    ) %>% unique() %>%
      mutate(residuals = quantile(residuals, var_qnt)) %>%
      unique()
    
    
    
    if(thresh_qnt == 0.95){
      thresh_ex = read_csv(paste0("data/processed/thresh_exceedance_lambda_mintp_",var_qnt,".csv")) %>%
        dplyr::select(Station, year, thresh_exceedance_95_coast) %>% unique()
      
      obs_grid = obs_grid %>%
        left_join(temporal_covars) %>%
        left_join(thresh_ex)%>%
        left_join(spatial_covars)  %>%
        drop_na()
      
      obs_grid$scales_logged = obs_grid$scales_95_logged
      
      obs_grid$thresh_exceedance = obs_grid$thresh_exceedance_95_coast
      obs_grid$threshold = obs_grid$threshold_tn_l_o_95_w_coast
      
      
      
      
    }else if(thresh_qnt == 0.96){
      
      thresh_ex = read_csv(paste0("data/processed/thresh_exceedance_lambda_mintp_",var_qnt,".csv")) %>%
        dplyr::select(Station, year, thresh_exceedance_96_coast) %>% unique()
      
      obs_grid = obs_grid %>%
        left_join(temporal_covars) %>%
        left_join(thresh_ex)%>%
        left_join(spatial_covars)  %>%
        drop_na()
      
      obs_grid$scales_logged = obs_grid$scales_96_logged
      obs_grid$thresh_exceedance = obs_grid$thresh_exceedance_96_coast
      obs_grid$threshold = obs_grid$threshold_tn_l_o_96_w_coast
      
      
    }else if(thresh_qnt == 0.97){
      thresh_ex = read_csv(paste0("data/processed/thresh_exceedance_lambda_mintp_",var_qnt,".csv")) %>%
        dplyr::select(Station, year, thresh_exceedance_97_coast) %>% unique()
      
      obs_grid = obs_grid %>%
        left_join(temporal_covars) %>%
        left_join(thresh_ex)%>%
        left_join(spatial_covars)  %>%
        drop_na()
      
      obs_grid$scales_logged = obs_grid$scales_97_logged
      obs_grid$thresh_exceedance = obs_grid$thresh_exceedance_97_coast
      obs_grid$threshold = obs_grid$threshold_tn_l_o_97_w_coast
      
    }else{
      asfdgsafadafdgf()
    }
    
    
    
    fit = readRDS(paste0('output/bts_param_est/thresh_qnt_',thresh_qnt,"_model_", marg_mod))
    pred = fit %>% predict_gpd(data = obs_grid, model = marg_mod)
    obs_grid$scale = pred$scale
    obs_grid$shape = pred$shape
    
    
    # ---- Read in simulations
    my_simulations_extremes = c()
    for(i in seq(1, 100)){
      my_simulations_extremes = c(my_simulations_extremes, readRDS(paste0("output/simulations/simulations_on_obs_grid/true/thresh_qnt_", thresh_qnt,"_marg_mod_",marg_mod,"_run_",i)))
    }
    
    my_simulations_extremes = my_simulations_extremes[seq(25000)]
    
    # ---- standardise simulations to have cost = 1
    my_simulations_standardised = list()
    for (i in 1:length(my_simulations_extremes)) {
      this_cost = mean(my_simulations_extremes[[i]])
      my_simulations_standardised[[i]] = my_simulations_extremes[[i]]/this_cost
    }
    
    # --- maximum frechet margin at each of my simulation sites
    max_at_each_site = c()
    for(s in seq(length(my_simulations_standardised[[1]]))){
      max_at_each_site = c(max_at_each_site, lapply(my_simulations_standardised, "[[", s) %>% unlist %>% max)
    }
    
    m = length(my_simulations_standardised)
    
    
    r_thresh = readRDS(paste0("output/r_thres_", marg_mod, "_thresh_quant_",thresh_qnt))
    m = length(my_simulations_standardised)
    L = 300
    sampled_rs = evd::rgpd(n=L, loc = 1, scale = 1, shape = 1)
    
    res_1950 = c()
    res_2020 = c()
    
    for(temp_i_want in seq(5, 20, by = 0.5)){
      
      print(temp_i_want)
      
      
      obs_grid$frechet = -1/log(1-obs_grid$thresh_exceedance*(1-evd::pgpd(temp_i_want - obs_grid$threshold,  scale = obs_grid$scale, shape =  obs_grid$shape[1] )) )
      
      # for(mnt in c(1,2,12)){
      
      T_2020 = obs_grid %>% filter(year == 2020) %>% pull(frechet)
      T_2022 = obs_grid %>% filter(year == 2022) %>% pull(frechet)
      T_1950 = obs_grid %>% filter(year == 1950) %>% pull(frechet)
      
      
      T_1950[T_1950 == -Inf] = Inf
      T_2020[T_2020 == -Inf] = Inf
      T_2022[T_2022 == -Inf] = Inf
      
      b_2020 = min(T_2020/max_at_each_site)
      b_2022 = min(T_2022/max_at_each_site)
      b_1950 = min(T_1950/max_at_each_site)
      
      
      m = length(my_simulations_standardised)
      
      above_1950 = 0
      above_2020 = 0
      above_2022 = 0
      
      for(j in seq(length(my_simulations_standardised))){
        # for each imp samp
        for(k in seq(L)){
          
          tmp = my_simulations_standardised[[j]]*sampled_rs[k]
          above_2022 = above_2022 + sum(sum(tmp*b_2022 > T_2022) > 0)
          above_2020 = above_2020 + sum(sum(tmp*b_2020 > T_2020) > 0)
          above_1950 = above_1950 + sum(sum(tmp*b_1950 > T_1950) > 0)
        }
      }
      
      tibble(temp = temp_i_want,
             p_1950 = above_1950 / (length(my_simulations_standardised)*L*b_1950),
             p_2020 = above_2020 / (length(my_simulations_standardised)*L*b_2020),
             p_2022 = above_2022 / (length(my_simulations_standardised)*L*b_2022)) %>%
        write_csv(paste0("output/prob_extreme_temp_imp_samp_thresh_quant_",thresh_qnt,"_mar_mod_",marg_mod,"_var_qnt_",var_qnt,".csv"), append = T)
      # }
    }
  }
}


# --- 30-35 mins
# job::job({run_imp_samp(marg_mods = c('E2'), var_qnt = 0.5, thresh_qnt = 0.96)})
# job::job({run_imp_samp(marg_mods = c('E2'), var_qnt = 0.5, thresh_qnt = 0.95)})
# job::job({run_imp_samp(marg_mods = c('E2'), var_qnt = 0.5, thresh_qnt = 0.97)})
# 
# job::job({run_imp_samp(marg_mods = c('E2'), var_qnt = 0.1, thresh_qnt = 0.96)})
# job::job({run_imp_samp(marg_mods = c('E2'), var_qnt = 0.9, thresh_qnt = 0.96)})
# job::job({run_imp_samp(marg_mods = c('E1'), var_qnt = 0.5, thresh_qnt = 0.96)})



#plot results


plt = rbind(read_csv(paste0("output/prob_extreme_temp_imp_samp_thresh_quant_",0.96,"_mar_mod_",'E2',"_var_qnt_",0.1,"_bootstrapped.csv"),
         col_names = c('temp', 'bts', 'p_1950', 'p_2020', 'p2022')) %>%
        mutate(qnt = "0.1"),
      read_csv(paste0("output/prob_extreme_temp_imp_samp_thresh_quant_",0.96,"_mar_mod_",'E2',"_var_qnt_",0.5,"_bootstrapped.csv"),
         col_names = c('temp', 'bts', 'p_1950', 'p_2020', 'p2022'))%>%
        mutate(qnt = "0.5"),
      read_csv(paste0("output/prob_extreme_temp_imp_samp_thresh_quant_",0.96,"_mar_mod_",'E2',"_var_qnt_",0.9,"_bootstrapped.csv"),
         col_names = c('temp', 'bts', 'p_1950', 'p_2020', 'p2022'))%>%
        mutate(qnt = "0.9")) %>%
  mutate(temp = -temp) %>%
  group_by(temp, qnt) %>%
  summarise(p_1950_lower = quantile(p_1950, 0.025),
            p_1950_upper = quantile(p_1950, 0.975),
            p_2022_lower = quantile(p2022, 0.025),
            p_2022_upper = quantile(p2022, 0.975)) %>%
  ggplot()+
  geom_ribbon(aes(temp,  ymin = 1/(90.25*p_1950_lower), ymax = 1/(90.25*p_1950_upper), fill = '1950'), alpha = 0.25)+
  geom_ribbon(aes(temp,  ymin = 1/(90.25*p_2022_lower), ymax = 1/(90.25*p_2022_upper), fill = '2022'), alpha = 0.25)+
  geom_line(data = read_csv(paste0("output/prob_extreme_temp_imp_samp_thresh_quant_",0.96,"_mar_mod_",'E2',"_var_qnt_",0.1,".csv"),
                            col_names = c('temp', 'p_1950', 'p_2020', 'p2022')) %>% mutate(qnt = "0.1"), 
            aes(-temp, 1/(90.25*p2022), col = '2022'), linetype = 'dashed', size=1)+
  geom_line(data = read_csv(paste0("output/prob_extreme_temp_imp_samp_thresh_quant_",0.96,"_mar_mod_",'E2',"_var_qnt_",0.1,".csv"),
                            col_names = c('temp', 'p_1950', 'p_2020', 'p2022')) %>% filter(temp<17) %>% mutate(qnt = "0.1"), 
            aes(-temp, 1/(90.25*p_1950), col = '1950'), size=1)+
  geom_line(data = read_csv(paste0("output/prob_extreme_temp_imp_samp_thresh_quant_",0.96,"_mar_mod_",'E2',"_var_qnt_",0.5,".csv"),
                            col_names = c('temp', 'p_1950', 'p_2020', 'p2022')) %>% mutate(qnt = "0.5"), 
            aes(-temp, 1/(90.25*p2022), col = '2022'), linetype = 'dashed', size=1)+
  geom_line(data = read_csv(paste0("output/prob_extreme_temp_imp_samp_thresh_quant_",0.96,"_mar_mod_",'E2',"_var_qnt_",0.5,".csv"),
                            col_names = c('temp', 'p_1950', 'p_2020', 'p2022')) %>% mutate(qnt = "0.5"), 
            aes(-temp, 1/(90.25*p_1950), col = '1950'), size=1)+
  geom_line(data = read_csv(paste0("output/prob_extreme_temp_imp_samp_thresh_quant_",0.96,"_mar_mod_",'E2',"_var_qnt_",0.9,".csv"),
                            col_names = c('temp', 'p_1950', 'p_2020', 'p2022')) %>% mutate(qnt = "0.9"), 
            aes(-temp, 1/(90.25*p2022), col = '2022'), linetype = 'dashed', size=1)+
  geom_line(data = read_csv(paste0("output/prob_extreme_temp_imp_samp_thresh_quant_",0.96,"_mar_mod_",'E2',"_var_qnt_",0.9,".csv"),
                            col_names = c('temp', 'p_1950', 'p_2020', 'p2022')) %>% mutate(qnt = "0.9"), 
            aes(-temp, 1/(90.25*p_1950), col = '1950'), size=1)+
  scale_x_continuous(limits = c(-20, -5),
                     breaks = c(-5, -10,  -15, -20),
                     label = paste0(c(-5, -10,  -15, -20),"°C"))+
  scale_y_log10(breaks = c(0.1, 1, 10,  100,  1000),
                labels = c(0.1, 1, 10,  100,  1000))+
  facet_wrap(~qnt)+
  coord_flip()+
  coord_flip(ylim = c(0.01, 1000))+
  labs(x = "Temperature",
       y = "Return period (Years)",
       fill = 'Year', col = "Year")+
  theme_minimal(12)+
  theme(legend.position = 'none',
        strip.text.x = element_blank(),
        axis.title.y = element_text(angle = 0, vjust = 0.5))

ggsave(plot = plt, filename = "output/figs/spatial_E2.pdf", height = 3, width = 11)
