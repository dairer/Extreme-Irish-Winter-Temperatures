# This script carries out a log liklihood test to investigate if a constant shape parameter is supported

obs_dat = read_csv('data/processed/obs_data_mintp_winter.csv') 
obs_dat$threshold =  obs_dat$threshold_tn_l_o_96_w_coast

source('src/models/marginal/gpd_models.R')


# log liklihood of the GPD
gpdll <- function(par){
  if(!fix_shape) shape_est =  par[length(par)]  
  # model E2
  scale_est = exp(par[1] + par[2]*extreme_data$dist_sea_logged + par[3]*extreme_data$loess_temp_anom + par[4]*extreme_data$residuals)
  if(any(scale_est <= 0)) return(-2^100)
  
  if((shape_est < 0) & any(excess_dat > (-scale_est/shape_est))) return(-2^100)
  
  sum(evd::dgpd(x = excess_dat, loc=0, scale = scale_est, shape=shape_est, log=T))
}



fit_gpd = function(data, initial_pars = NULL, fix_shape = FALSE, shape_est = NA){
  # covaraites for GPD fir
  excess_dat <<- data$excess
  loess_temp_anom <<- data$loess_temp_anom
  dist_sea <<- data$dist_sea_logged
  resids <<- data$residuals
  
  if(fix_shape){
    fix_shape <<- TRUE
    shape_est <<- shape_est
  }else{
    fix_shape <<- FALSE
  }

  if(fix_shape == TRUE){
    initial_pars = initial_pars[-length(initial_pars)]
  }

  fit = optim(par = initial_pars, fn = gpdll, control = list(fnscale = -1))
  fit
}

# get extreme observations
extreme_data = obs_dat %>%
  mutate(excess = temp - threshold) %>%
  filter(excess > 0)

# fit gpd to each site, allowing for a fixed and free shape parameter estimate
all_fits = extreme_data %>%
  group_by(Station) %>%
  group_map(~{
    
    print(.x$Station[1])
    .x %>%
      dplyr::select(Station, Long, Lat) %>%
      unique() %>%
      mutate(fixed = fit_gpd(.x, c(1,0,0,0,0), fix_shape = T, shape_est = -0.07934607)$value,
             free = fit_gpd(.x, c(1,0,0,0,0))$value)

  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()

# calculate LL ratio for each site
all_fits$LL_ratio = -2*log(all_fits$free/all_fits$fixed)


# plot ratios
all_fits %>%
  left_join(extreme_data %>%
              group_by(Station) %>%
              summarise(n = n())) %>%
  ggplot()+
  geom_point(aes(n, LL_ratio), size = 1)+
  geom_hline(yintercept = 0.05, col = 'red')