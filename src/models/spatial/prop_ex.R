# calcualte proportion of exceedence of sites, and uncondition by scaling probability

rm(list=ls())
library(tidyverse)

prop_ex = function(marg_mod, temp_conditioned_on, var_qnt, thresh_qnt){
  
    obs_sites = read_csv("data/processed/obs_data_mintp_winter.csv") %>%
      dplyr::select(Station, Long.projected, Lat.projected, dist_sea_logged) %>%
      unique()
    
    grid_simulated = as.data.frame(read_csv("data/processed/obs_grid_simulated_on.csv")) %>%
      left_join(obs_sites) %>% as.tibble()
    
    
    for(tmp in temp_conditioned_on){
      yr = 2022
      
      if(file.exists(paste0("output/simulations/rescaled_simulations_on_obs_grid/true/thresh_qnt_", thresh_qnt,"_marg_mod_",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp,"_qnt_", var_qnt))){
        my_simulations = readRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/true/thresh_qnt_", thresh_qnt,"_marg_mod_",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp,"_qnt_", var_qnt))
        
        prop_ex_2022 = c()
        
        
        frechet_val = grid_simulated %>%
          left_join(read_csv(paste0("output/obs_sites_extreme_temps_frechet_scale_thresh_qnt_",thresh_qnt,"_model_",marg_mod,"_qnt_",var_qnt,".csv")) %>% 
                      filter(temp == tmp, year == yr) %>%
                      unique()) %>% 
          pull(frechet_value)
        
        
        frechet_val[frechet_val == -Inf] = Inf
        
        num_exceed_tmp = my_simulations %>%
          map(~{ sum(.x > frechet_val)}) %>%
          unlist()
        
        prop_ex_2022 =  mean(num_exceed_tmp[num_exceed_tmp>0]/108)
        
        
        yr = 1950
        my_simulations = readRDS(paste0("output/simulations/rescaled_simulations_on_obs_grid/true/thresh_qnt_", thresh_qnt,"_marg_mod_",marg_mod,"_yr_",yr, "_min_temp_conditioned_on_",tmp,"_qnt_", var_qnt))
        
        frechet_val = grid_simulated %>%
          left_join(read_csv(paste0("output/obs_sites_extreme_temps_frechet_scale_thresh_qnt_",thresh_qnt,"_model_",marg_mod,"_qnt_",var_qnt,".csv")) %>% 
                      filter(temp == tmp, year == yr) %>%
                      unique()) %>% 
          pull(frechet_value)
        
        frechet_val[frechet_val == -Inf] = Inf
        
        num_exceed_tmp = my_simulations %>% 
          map(~{ sum(.x > frechet_val)}) %>%
          unlist()
        
        prop_ex_1950 = mean(num_exceed_tmp[num_exceed_tmp>0]/108)
        
        print("---------")
        print(prop_ex_1950)
        
        tibble(temp = tmp, prop_ex_1950, prop_ex_2022) %>%
          write_csv(paste0("output/simulations/simulation_summary/new_prop_exceedance_model_true_thresh_qnt_",thresh_qnt,"_marg_mod_",marg_mod,"_qnt_", var_qnt,".csv"), append = T)
        
      }
    }
}


# job::job({prop_ex('m22', seq(6, 12), 0.1, 0.96)})
# job::job({prop_ex('m22', seq(6, 12), 0.5, 0.96)})
# job::job({prop_ex('m22', seq(6, 12), 0.9, 0.96)})

# job::job({prop_ex('m22', seq(6, 12), 0.5, 0.95)})
# job::job({prop_ex('m22', seq(6, 12), 0.5, 0.97)})
# job::job({prop_ex('m21', seq(6, 12), 0.5, 0.96)})


read_csv(paste0("output/simulations/simulation_summary/new_prop_exceedance_model_true_thresh_qnt_",0.95,"_marg_mod_",'m22',"_qnt_", 0.9,".csv"),
         col_names = c('temp', 'p1950', 'p2022')) %>%
  pivot_longer(-temp) %>%
  ggplot()+
  geom_line(aes(temp, value, col = name))


marg_mod = 'm22'
thresh_qnt = 0.96
var_qnt = 0.5
bts_dat = rbind(read_csv(paste0("output/simulations/simulation_summary/new_bootstrapped_prop_exceedance_model_true_thresh_qnt_",thresh_qnt,"_marg_mod_",marg_mod,"_qnt_", 0.1,".csv"),
                         col_names = c('bts', 'temp', 'prop_ex_1950', 'prop_ex_2022')) %>% mutate(qnt = 0.1),
                read_csv(paste0("output/simulations/simulation_summary/new_bootstrapped_prop_exceedance_model_true_thresh_qnt_",thresh_qnt,"_marg_mod_",marg_mod,"_qnt_", 0.5,".csv"),
                         col_names = c('bts', 'temp', 'prop_ex_1950', 'prop_ex_2022')) %>% mutate(qnt = 0.5),
                read_csv(paste0("output/simulations/simulation_summary/new_bootstrapped_prop_exceedance_model_true_thresh_qnt_",thresh_qnt,"_marg_mod_",marg_mod,"_qnt_", 0.9,".csv"),
                         col_names = c('bts', 'temp', 'prop_ex_1950', 'prop_ex_2022')) %>% mutate(qnt = 0.9))


bts_dat = bts_dat %>%
  group_by(temp, qnt) %>%
  summarise(prop_ex_2022_u = quantile(prop_ex_2022, 0.975, na.rm = T),
            prop_ex_2022_l = quantile(prop_ex_2022, 0.025, na.rm = T),
            prop_ex_1950_u = quantile(prop_ex_1950, 0.975, na.rm = T),
            prop_ex_1950_l = quantile(prop_ex_1950, 0.025, na.rm = T))


true_dat = rbind(read_csv(paste0("output/simulations/simulation_summary/new_prop_exceedance_model_true_thresh_qnt_",thresh_qnt,"_marg_mod_",marg_mod,"_qnt_", 0.1,".csv"),
                    col_name = c('temp', 'prop_ex_1950', 'prop_ex_2022')) %>% mutate(qnt = 0.1),
                 read_csv(paste0("output/simulations/simulation_summary/new_prop_exceedance_model_true_thresh_qnt_",thresh_qnt,"_marg_mod_",marg_mod,"_qnt_", 0.5,".csv"),
                          col_name = c('temp', 'prop_ex_1950', 'prop_ex_2022')) %>% mutate(qnt = 0.5),
                 read_csv(paste0("output/simulations/simulation_summary/new_prop_exceedance_model_true_thresh_qnt_",thresh_qnt,"_marg_mod_",marg_mod,"_qnt_", 0.9,".csv"),
                          col_name = c('temp', 'prop_ex_1950', 'prop_ex_2022')) %>% mutate(qnt = 0.9))


prob_observing_T = rbind(read_csv(paste0("output/prob_extreme_temp_imp_samp_thresh_quant_",thresh_qnt,"_mar_mod_",marg_mod,"_var_qnt_",0.1,".csv"),
                                  col_names = c('temp', 'p_1950', 'p_2020', 'p_2022')) %>% mutate(qnt = 0.1),
                         read_csv(paste0("output/prob_extreme_temp_imp_samp_thresh_quant_",thresh_qnt,"_mar_mod_",marg_mod,"_var_qnt_",0.5,".csv"),
                                  col_names = c('temp', 'p_1950', 'p_2020', 'p_2022')) %>% mutate(qnt = 0.5),
                         read_csv(paste0("output/prob_extreme_temp_imp_samp_thresh_quant_",thresh_qnt,"_mar_mod_",marg_mod,"_var_qnt_",0.9,".csv"),
                                  col_names = c('temp', 'p_1950', 'p_2020', 'p_2022')) %>% mutate(qnt = 0.9))


prob_observing_T = prob_observing_T%>%
  dplyr::select(temp,p_1950,p_2022, qnt) %>%
  pivot_longer(-c(temp,qnt)) %>%
  mutate(year = str_remove(name,"p_")) 

prob_observing_T = rbind(prob_observing_T %>% filter(name == 'p_1950') %>% 
                           mutate(typ = 'actual') %>% dplyr::select(-name),
                         prob_observing_T %>% filter(name == 'p_2022') %>% 
                           mutate(typ = 'actual') %>% dplyr::select(-name)) %>%
  mutate(year = as.numeric(year)) %>%
  rename(unconditional_factor = value)


bts_dat = bts_dat %>%
  pivot_longer(-c(temp,qnt)) %>%
  mutate(year = str_remove(name,"prop_ex_")) %>%
  mutate(year = str_remove(year,"_u")) %>%
  mutate(year = str_remove(year,"_l"))



bts_dat = rbind(bts_dat %>% filter(name == 'prop_ex_2022_u') %>% mutate(typ = 'upper') %>% dplyr::select(-name),
                bts_dat %>% filter(name == 'prop_ex_2022_l') %>% mutate(typ = 'lower') %>% dplyr::select(-name),
                bts_dat %>% filter(name == 'prop_ex_1950_u') %>% mutate(typ = 'upper') %>% dplyr::select(-name),
                bts_dat %>% filter(name == 'prop_ex_1950_l') %>% mutate(typ = 'lower') %>% dplyr::select(-name)) %>%
  mutate(year = as.numeric(year))


bts_dat = bts_dat %>%
  left_join(prob_observing_T %>% dplyr::select(-typ) %>% unique, by= c('year', 'temp', 'qnt'))

true_dat = true_dat %>%
  pivot_longer(-c(temp,qnt)) %>%
  mutate(year = str_remove(name,"prop_ex_")) 


#true_dat$typ = 'actual'
true_dat$year = as.numeric(true_dat$year)
true_dat = true_dat %>%
  left_join(prob_observing_T)


bts_dat$unconditioned = bts_dat$value*bts_dat$unconditional_factor
true_dat$unconditioned = true_dat$value*true_dat$unconditional_factor


bts_dat = rbind(bts_dat %>%
                  dplyr::select(temp, value, year, typ, qnt) %>%
                  pivot_wider(names_from = c(year, typ), values_from = value) %>%
                  mutate(lab = 'conditioned'),
                bts_dat %>%
                  dplyr::select(temp, unconditioned, year, typ, qnt) %>%
                  pivot_wider(names_from = c(year, typ), values_from = unconditioned) %>%
                  mutate(lab = 'unconditioned'))



true_dat = rbind(true_dat %>%
                   dplyr::select(temp, value, year, typ, qnt) %>%
                   pivot_wider(names_from = c(year, typ), values_from = value) %>%
                   mutate(lab = 'conditioned'),
                 true_dat %>%
                   dplyr::select(temp, unconditioned, year, typ, qnt) %>%
                   pivot_wider(names_from = c(year, typ), values_from = unconditioned) %>%
                   mutate(lab = 'unconditioned'))


# plot results
library(scales)
library(latex2exp)

plt = gridExtra::grid.arrange(bts_dat %>% 
                                filter(lab == 'conditioned') %>%
                                ggplot()+
                                geom_ribbon(aes(-temp, ymin = `2022_lower`, ymax = `2022_upper`, fill = '2022'), alpha = 0.25)+
                                geom_ribbon(aes(-temp, ymin = `1950_lower`, ymax = `1950_upper`, fill = '1950'), alpha = 0.25)+
                                geom_line(data = true_dat %>% filter(lab == 'conditioned'),
                                            aes(-temp, `1950_actual`, col = '1950', linetype = '1950'))+
                                geom_line(data = true_dat %>% filter(lab == 'conditioned'),
                                            aes(-temp, `2022_actual`, col = '2022', linetype = '2022'))+
                                facet_wrap(~qnt)+
                                scale_x_continuous(limits = c(-12, -6),
                                                   breaks = c(-12, -10, -8, -6),
                                                   label = (paste0(c(-12, -10, -8, -6),"°C")))+
                                scale_x_reverse()+
                                labs(x = "Temperature",
                                     y = expression(E[o]),
                                     col = "Year",
                                     fill = "Year",
                                     linetype = "Year")+
                                theme_minimal(12)+
                                theme(axis.text.x = element_blank(),
                                      axis.title = element_blank(),
                                      axis.title.y = element_text(angle = 0, vjust = 0.5),
                                      legend.position = 'none',
                                      strip.text.x = element_text(size=0)),
                              bts_dat %>% 
                                filter(lab == 'unconditioned') %>%
                                ggplot()+
                                geom_ribbon(aes(-temp, ymin = `2022_lower`, ymax = `2022_upper`, fill = '2022'), alpha = 0.25)+
                                geom_ribbon(aes(-temp, ymin = `1950_lower`, ymax = `1950_upper`, fill = '1950'), alpha = 0.25)+
                                geom_line(data = true_dat %>% filter(lab == 'unconditioned'),
                                          aes(-temp, `1950_actual`, col = '1950', linetype = '1942'))+
                                geom_line(data = true_dat %>% filter(lab == 'unconditioned'),
                                          aes(-temp, `2022_actual`, col = '2022', linetype = '2020'))+
                                facet_wrap(~qnt)+
                                scale_x_continuous(limits = c(-12, -6),
                                                   breaks = c(-12, -10, -8, -6),
                                                   label = paste0(c(-12, -10, -8, -6),"°C"))+
                                scale_x_reverse()+
                                theme_minimal(12)+
                                theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
                                      legend.position = 'none',
                                      strip.text.x = element_text(size=0))+
                                labs(x = "Temperature",
                                     y = expression(E[o]),
                                     col = "Year",
                                     fill = "Year",
                                     linetype = "Year")+
                                scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                              labels = trans_format("log10", math_format(10^.x))),heights = c(0.9,1), nrow = 2)

ggsave(plot =  plt, paste0("output/figs/prop_exceedance_",marg_mod,".pdf"), height = 5, width =8)




