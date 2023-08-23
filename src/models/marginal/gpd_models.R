# This script contains code to fit GPD models as well as calculate estimated paramters (scale and shape) and return levels.

scale_models = function(par, model){
  # miscl
  if(model == 'M1')  scale_est = par[1]
  if(model == 'M2')  scale_est = par[1] + par[2]*clim_scale
  if(model == 'M3')  scale_est = par[1] + par[2]*dist_sea
  
  
  # --- Group A
  if(model == 'A1')   scale_est = par[1] + par[2]*loess_temp_anom
  if(model == 'A2')   scale_est = par[1] + par[2]*loess_temp_anom + par[3]*resids
  if(model == 'A3')   scale_est = par[1] + par[2]*loess_temp_anom + par[3]*yearly_resids
  if(model == 'A4')   scale_est = par[1] + par[2]*loess_temp_anom + par[3]*temp_anom_around_irel_resids
  if(model == 'A5')   scale_est = par[1] + par[2]*loess_temp_anom + par[3]*temp_anom_north_irel_resids
  if(model == 'A6')   scale_est = par[1] + par[2]*loess_temp_anom + par[3]*ao
  if(model == 'A7')   scale_est = par[1] + par[2]*loess_temp_anom + par[3]*ao_resid
  if(model == 'A8')   scale_est = par[1] + par[2]*loess_temp_anom + par[3]*ao_monthly
  if(model == 'A9')   scale_est = par[1] + par[2]*loess_temp_anom + par[3]*ao_monthly_resid
  if(model == 'A10')  scale_est = par[1] + par[2]*loess_temp_anom + par[3]*nao
  if(model == 'A11')  scale_est = par[1] + par[2]*loess_temp_anom + par[3]*nao_resid
  if(model == 'A12')  scale_est = par[1] + par[2]*loess_temp_anom + par[3]*nao_monthly
  if(model == 'A13')  scale_est = par[1] + par[2]*loess_temp_anom + par[3]*nao_monthly_resid
  
  # --- Group B
  if(model == 'B1')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom
  if(model == 'B2')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*resids
  if(model == 'B3')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*yearly_resids
  if(model == 'B4')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*temp_anom_around_irel_resids
  if(model == 'B5')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*temp_anom_north_irel_resids
  if(model == 'B6')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*ao
  if(model == 'B7')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*ao_resid
  if(model == 'B8')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*ao_monthly
  if(model == 'B9')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*ao_monthly_resid
  if(model == 'B10')  scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*nao
  if(model == 'B11')  scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*nao_resid
  if(model == 'B12')  scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*nao_monthly
  if(model == 'B13')  scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*nao_monthly_resid
  
  # --- Group C
  if(model == 'C1')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*dist_sea
  if(model == 'C2')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*dist_sea + par[5]*resids
  if(model == 'C3')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*dist_sea + par[5]*yearly_resids
  if(model == 'C4')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*dist_sea + par[5]*temp_anom_around_irel_resids
  if(model == 'C5')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*dist_sea + par[5]*temp_anom_north_irel_resids
  if(model == 'C6')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*dist_sea + par[5]*ao
  if(model == 'C7')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*dist_sea + par[5]*ao_resid
  if(model == 'C8')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*dist_sea + par[5]*ao_monthly
  if(model == 'C9')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*dist_sea + par[5]*ao_monthly_resid
  if(model == 'C10')  scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*dist_sea + par[5]*nao
  if(model == 'C11')  scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*dist_sea + par[5]*nao_resid
  if(model == 'C12')  scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*dist_sea + par[5]*nao_monthly
  if(model == 'C13')  scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*dist_sea + par[5]*nao_monthly_resid
  
  # --- Group D
  if(model == 'D1')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea
  if(model == 'D2')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*resids
  if(model == 'D3')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*yearly_resids
  if(model == 'D4')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*temp_anom_around_irel_resids
  if(model == 'D5')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*temp_anom_north_irel_resids
  if(model == 'D6')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*ao
  if(model == 'D7')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*ao_resid
  if(model == 'D8')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*ao_monthly
  if(model == 'D9')   scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*ao_monthly_resid
  if(model == 'D10')  scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*nao
  if(model == 'D11')  scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*nao_resid
  if(model == 'D12')  scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*nao_monthly
  if(model == 'D13')  scale_est = par[1] + par[2]*clim_scale + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*nao_monthly_resid
  
  # --- Group E
  if(model == 'E1')   scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom
  if(model == 'E2')   scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*resids
  if(model == 'E3')   scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*yearly_resids
  if(model == 'E4')   scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*temp_anom_around_irel_resids
  if(model == 'E5')   scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*temp_anom_north_irel_resids
  if(model == 'E6')   scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*ao
  if(model == 'E7')   scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*ao_resid
  if(model == 'E8')   scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*ao_monthly
  if(model == 'E9')   scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*ao_monthly_resid
  if(model == 'E10')  scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*nao
  if(model == 'E11')  scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*nao_resid
  if(model == 'E12')  scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*nao_monthly
  if(model == 'E13')  scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*nao_monthly_resid
  
  
  # --- Group F
  if(model == 'F1')   scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea
  if(model == 'F2')   scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*resids
  if(model == 'F3')   scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*yearly_resids
  if(model == 'F4')   scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*temp_anom_around_irel_resids
  if(model == 'F5')   scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*temp_anom_north_irel_resids
  if(model == 'F6')   scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*ao
  if(model == 'F7')   scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*ao_resid
  if(model == 'F8')   scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*ao_monthly
  if(model == 'F9')   scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*ao_monthly_resid
  if(model == 'F10')  scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*nao
  if(model == 'F11')  scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*nao_resid
  if(model == 'F12')  scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*nao_monthly
  if(model == 'F13')  scale_est = par[1] + par[2]*dist_sea + par[3]*loess_temp_anom + par[4]*loess_temp_anom*dist_sea + par[5]*nao_monthly_resid

  # --- Group G
  if(model == 'G1')   scale_est = par[1] + par[2]*loess_temp_anom + par[3]*loess_temp_anom*dist_sea
  if(model == 'G2')   scale_est = par[1] + par[2]*loess_temp_anom + par[3]*loess_temp_anom*dist_sea + par[4]*resids
  if(model == 'G3')   scale_est = par[1] + par[2]*loess_temp_anom + par[3]*loess_temp_anom*dist_sea + par[4]*yearly_resids
  if(model == 'G4')   scale_est = par[1] + par[2]*loess_temp_anom + par[3]*loess_temp_anom*dist_sea + par[4]*temp_anom_around_irel_resids
  if(model == 'G5')   scale_est = par[1] + par[2]*loess_temp_anom + par[3]*loess_temp_anom*dist_sea + par[4]*temp_anom_north_irel_resids
  if(model == 'G6')   scale_est = par[1] + par[2]*loess_temp_anom + par[3]*loess_temp_anom*dist_sea + par[4]*ao
  if(model == 'G7')   scale_est = par[1] + par[2]*loess_temp_anom + par[3]*loess_temp_anom*dist_sea + par[4]*ao_resid
  if(model == 'G8')   scale_est = par[1] + par[2]*loess_temp_anom + par[3]*loess_temp_anom*dist_sea + par[4]*ao_monthly
  if(model == 'G9')   scale_est = par[1] + par[2]*loess_temp_anom + par[3]*loess_temp_anom*dist_sea + par[4]*ao_monthly_resid
  if(model == 'G10')  scale_est = par[1] + par[2]*loess_temp_anom + par[3]*loess_temp_anom*dist_sea + par[4]*nao
  if(model == 'G11')  scale_est = par[1] + par[2]*loess_temp_anom + par[3]*loess_temp_anom*dist_sea + par[4]*nao_resid
  if(model == 'G12')  scale_est = par[1] + par[2]*loess_temp_anom + par[3]*loess_temp_anom*dist_sea + par[4]*nao_monthly
  if(model == 'G13')  scale_est = par[1] + par[2]*loess_temp_anom + par[3]*loess_temp_anom*dist_sea + par[4]*nao_monthly_resid
  
  return(exp(scale_est))
}


# log liklihood of GPD
gpdll <- function(par){
  if(!fix_shape) shape_est =  par[length(par)]  
  scale_est = scale_models(par, model)
  if(any(scale_est <= 0)) return(-2^100)
  
  if((shape_est < 0) & any(excess_dat > (-scale_est/shape_est))) return(-2^100)
  
  sum(evd::dgpd(x = excess_dat, loc=0, scale = scale_est, shape=shape_est, log=T))
}



fit_gpd = function(data, model, initial_pars = NULL, fix_shape = FALSE, shape_est = NA, hess = F){
  excess_dat <<- data$excess
  clim_scale <<- data$scales_logged
  loess_temp_anom <<- data$loess_temp_anom
  dist_sea <<- data$dist_sea_logged
  resids <<- data$residuals
  
  ao <<- data$ao
  ao_monthly <<- data$ao_monthly
  ao_resid <<- data$ao_resid
  ao_monthly_resid <<- data$ao_monthly_resid
  
  nao <<- data$nao
  nao_monthly <<- data$nao_monthly
  nao_resid <<- data$nao_resid
  nao_monthly_resid <<- data$nao_monthly_resid
  yearly_resids <<- data$yearly_residuals

  
  temp_anom_around_irel_resids <<- data$temp_anom_around_irel_resids
  temp_anom_north_irel_resids <<- data$temp_anom_north_irel_resids
  model <<- model
  
  
  if(fix_shape){
    fix_shape <<- TRUE
    shape_est <<- shape_est
  }else{
    fix_shape <<- FALSE
  }
  
  # set starting point of optimisation
  if(is.null(initial_pars)){
    
    if(model == 'M1')  initial_pars = c(1)
    if(model == 'M2')  initial_pars = c(1, 0)
    if(model == 'M3')  initial_pars = c(1, 0)
    
    
    
    
    if(grepl("A", model)){
      if(model == 'A1') initial_pars = c(1, 0)
      else initial_pars = c(1, 0, 0)
    }
    
    
    if(grepl("B", model)){
      if(model == 'B1') initial_pars = c(1, 0, 0)
      else initial_pars = c(1, 0, 0, 0)
    }
    
    
    if(grepl("C", model)){
      if(model == 'C1') initial_pars = c(1, 0, 0, 0)
      else initial_pars = c(1, 0, 0, 0, 0)
    }
    
    
    if(grepl("D", model)){
      if(model == 'D1') initial_pars  = c(1, 0, 0, 0)
      else initial_pars = c(1, 0, 0, 0, 0)
    }
    
    
    if(grepl("E", model)){
      if(model == 'E1') initial_pars = c(1, 0, 0)
      else initial_pars = c(1, 0, 0, 0)
    }
    
    
    if(grepl("F", model)){
      if(model == 'F1') initial_pars = c(1, 0, 0, 0)
      else initial_pars = c(1, 0, 0, 0, 0)
    }
    
    
    if(grepl("G", model)){
      if(model == 'G1') initial_pars = c(1, 0, 0)
      else initial_pars = c(1, 0, 0, 0)
    }
    

    initial_pars = c(initial_pars, 0) # add param for shape
    
    if(fix_shape == TRUE){
      # remove last digit from initial pars
      initial_pars = initial_pars[-length(initial_pars)]
    }
    
  }
  
  
  # calculate uncertainty from hessian matrix
  if(hess){
    fit = optim(par = initial_pars, fn = gpdll, control = list(fnscale = -1), hessian = T)
    
    # --- estimate se from hessian
    
    se = -fit$hessian %>%
      solve() %>%
      diag %>%
      sqrt()
    se = 1.96*se
    
    fit = list(fit$par, se = 1.96*se)
    
  }else{
    fit = optim(par = initial_pars, fn = gpdll, control = list(fnscale = -1))$par
  }
  return(fit)
}


# calcualte scale and shape from estimated GPD
predict_gpd = function(estimates_pars, model, data, fix_shape = F, shape_est = NA){
  
  clim_scale <<- data$scales_logged
  loess_temp_anom <<- data$loess_temp_anom
  dist_sea <<- data$dist_sea_logged
  resids <<- data$residuals
  
  ao <<- data$ao
  ao_monthly <<- data$ao_monthly
  ao_resid <<- data$ao_resid
  ao_monthly_resid <<- data$ao_monthly_resid
  yearly_resids <<- data$yearly_residuals
  
  nao <<- data$nao
  nao_monthly <<- data$nao_monthly
  nao_resid <<- data$nao_resid
  nao_monthly_resid <<- data$nao_monthly_resid
  temp_anom_around_irel_resids <<- data$temp_anom_around_irel_resids
  temp_anom_north_irel_resids <<- data$temp_anom_north_irel_resids
  

  model <<- model
  
  if(!fix_shape) shape_est =  estimates_pars[length(estimates_pars)]
  else shape_est = shape_est
  tibble(scale = scale_models(estimates_pars, model), shape = shape_est)
}


# calculate return level for estimated GPD
rl_gpd = function(estimates_pars, model, data, rl_quantile){
  clim_scale <<- data$scales_logged
  loess_temp_anom <<- data$loess_temp_anom
  dist_sea <<- data$dist_sea_logged
  resids <<- data$residuals
  
  
  threshold = data$threshold
  ao <<- data$ao
  ao_monthly <<- data$ao_monthly
  ao_resid <<- data$ao_resid
  ao_monthly_resid <<- data$ao_monthly_resid
  yearly_resids <<- data$yearly_residuals
  
  nao <<- data$nao
  nao_monthly <<- data$nao_monthly
  nao_resid <<- data$nao_resid
  nao_monthly_resid <<- data$nao_monthly_resid
  temp_anom_around_irel_resids <<- data$temp_anom_around_irel_resids
  temp_anom_north_irel_resids <<- data$temp_anom_north_irel_resids
  
  interaction1 <<- data$interaction1
  interaction2 <<- data$interaction2
  
  model <<- model
  
  shape =  estimates_pars[length(estimates_pars)]
  scale = scale_models(estimates_pars, model)
  return(data$threshold + scale * ((1-rl_quantile)^(-shape) - 1)/shape)
}
