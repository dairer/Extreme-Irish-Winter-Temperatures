# this script plots a qqplot on standardised exponential margins, comparing models E2 and E1



rm(list = ls())
source('src/models/marginal/gpd_models.R')
library(tidyverse)
library(gridtext)

ex_dat = read_csv("data/processed/obs_data_mintp_winter.csv") 

# get exceedence data
ex_dat = ex_dat %>%
  mutate(excess = temp - threshold_tn_l_o_96_w_coast)%>%
  filter(excess>0) 

ex_dat$scales_logged = ex_dat$scales_96


# fit model E1 and transform to uniform using probability integral transform and then to exponential
fit = fit_gpd(data = ex_dat, model = "E1") %>%
  predict_gpd(estimates_pars = ., model = "E1", data = ex_dat)
ex_dat$scale_E1 =  fit$scale
ex_dat$shape_E1 =  fit$shape
ex_dat$pred_E1  = evd::pgpd(ex_dat$excess, loc = 0, scale = ex_dat$scale_E1, shape = ex_dat$shape_E1[1])
ex_dat$pred_E1 = sort(qexp(ex_dat$pred_E1))

# fit model E2 and transform to uniform using probability integral transform and then to exponential
fit = fit_gpd(data = ex_dat, model = "E2") %>%
  predict_gpd(estimates_pars = ., model = "E2", data = ex_dat)
ex_dat$scale_E2=  fit$scale
ex_dat$shape_E2 =  fit$shape
ex_dat$pred_E2  = evd::pgpd(ex_dat$excess, loc = 0, scale = ex_dat$scale_E2, shape = ex_dat$shape_E2[1])
sum(evd::dgpd(ex_dat$excess, loc = 0, scale = ex_dat$scale_E2, shape = ex_dat$shape_E2[1], log = T))
ex_dat$LL = (evd::dgpd(ex_dat$excess, loc = 0, scale = ex_dat$scale_E2, shape = ex_dat$shape_E2[1], log = T))
ex_dat$pred_E2 = sort(qexp(ex_dat$pred_E2))

# calcualte uncertainties of quantiles
num_reps = nrow(ex_dat)
exp_ci = c()
for(i in seq(100)){
  exp_ci = rbind(exp_ci, sort(rexp(n = num_reps)))
}
lower_ci = exp_ci %>% apply(MARGIN = 2, FUN = quantile, 0.02)
upper_ci = exp_ci %>% apply(MARGIN = 2, FUN = quantile, 0.98) 
ex_dat$rank = exp_ci %>% apply(MARGIN = 2, FUN = mean) 
ex_dat$lower_ci = lower_ci
ex_dat$upper_ci = upper_ci

# plot both qq plots
plt = gridExtra::grid.arrange(ex_dat %>%
                          ggplot()+
                          geom_line(aes(y = sort(pred_E1),x = sort(rank)))+
                          geom_ribbon(aes(y = pred_E1, xmin= lower_ci, xmax = upper_ci), alpha = 0.3, fill = 'forestgreen')+
                          geom_line(aes(y = sort(pred_E1),x =  sort(lower_ci)), alpha = 0.5)+
                          geom_line(aes(y = sort(pred_E1), x =  sort(upper_ci)), alpha = 0.5)+
                          geom_abline(col = 'magenta',linetype = 'longdash')+
                          theme_minimal()+
                          theme(axis.title = element_blank())+
                          xlim(c(0,13))+
                          ylim(c(0,13)),
                        ex_dat %>%
                          ggplot()+
                          geom_line(aes(y = sort(pred_E2),x = sort(rank)))+
                          geom_ribbon(aes(y = sort(pred_E2), xmin= sort(lower_ci), xmax = sort(upper_ci)), alpha = 0.3, fill = 'forestgreen')+
                          geom_line(aes(y = sort(pred_E2), x =  sort(lower_ci)), alpha = 0.5)+
                          geom_line(aes(y = sort(pred_E2), x =  sort(upper_ci)), alpha = 0.5)+
                          geom_abline(col = 'magenta',linetype = 'longdash')+
                          theme_minimal()+
                          theme(axis.title = element_blank())+
                          xlim(c(0,13))+
                          ylim(c(0,13)), nrow = 1, bottom = "Ideal quantiles", left = grid::textGrob("Estimated\nquantiles  ", rot = 0))


ggsave(plot = plt, filename = "output/figs/qqplot_mod_comp.pdf", width = 7, height = 3)