# code to test/debug uncertainty of GPD estimtes

library(tidyverse)

models_to_run = c("E1", "E2")
thresh_qnt = 0.96


for(m in models_to_run){
  assign(paste0('pars_', m),c(0, readRDS(paste0('output/bts_param_est/thresh_qnt_',thresh_qnt,"_model_", m))))
  assign(paste0('uncorrected_', m),read.csv(paste0('output/bts_param_est/thresh_qnt_',thresh_qnt,'_uncorrected_',m,".csv"), header = FALSE, sep = ' '))
  assign(paste0('corrected_', m),read.csv(paste0('output/bts_param_est/thresh_qnt_',thresh_qnt,'_corrected_',m,".csv"), header = FALSE, sep = ' '))
}


for(m in models_to_run){
  gridExtra::grid.arrange(get(paste0('uncorrected_', m)) %>%
                            dplyr::select(-V1) %>%
                            pivot_longer(everything()) %>%
                            ggplot()+
                            geom_density(aes(value))+
                            geom_vline(data = get(paste0('pars_', m)) %>% t %>%  as_tibble() %>%dplyr::select(-V1) %>%
                                         pivot_longer(everything()), aes(xintercept = value))+
                            facet_wrap(~name, scale = 'free', nrow = 1)+
                            theme_minimal(10)+
                            labs(title = m),
                          get(paste0('corrected_', m)) %>%
                            dplyr::select(-V1)%>%
                            pivot_longer(everything()) %>%
                            ggplot()+
                            geom_density(aes(value))+
                            geom_vline(data = get(paste0('pars_', m)) %>% t %>%  as_tibble() %>%dplyr::select(-V1) %>%
                                         pivot_longer(everything()), aes(xintercept = value))+
                            facet_wrap(~name, scale = 'free', nrow = 1)+
                            labs(y= "Density", x = "Parameter value")+
                            theme_minimal(10)+theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1)), heights = c(0.87,1))
}
