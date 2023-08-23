# code used to test if dependence is temporally stationary
rm(list=ls())
library(tidyverse)
library(doParallel)
library(foreach)

fit_rparp_true = function(marg_mod, num_clusters, nu, thresh_qnt){
  
  
  HRD_ll = function(dt, lcs, vr, conditioned.site = 1, zeta, gamma, eta){
    
    nlocs = nrow(lcs) # number of sites
    loc.id.pairs = expand.grid(seq(nlocs) ,seq(nlocs))
    loc.pairs = cbind(lcs[loc.id.pairs[,1],], lcs[loc.id.pairs[,2],]) # all pairs of locations
    
    dstncs = loc.pairs %>%
      apply(MARGIN = 1, FUN = function(x) sqrt((x[3] - x[1])^2 + (x[4] - x[2])^2))
    
    Lambda = (vr(dstncs)) %>% matrix(nrow = nlocs, byrow = T)
    conditioned.site.lambda = Lambda[conditioned.site,]
    Lambda = Lambda[-conditioned.site, -conditioned.site]
    
    psi = outer(conditioned.site.lambda[-1], conditioned.site.lambda[-1], "+") - Lambda # "covar" matrix
    inversePsi = psi %>% solve()  # inverse of "covarianve" matrix
    
    detPsi = determinant(psi)$modulus[1] # calculates log det by default
    
    omegas = log(t(sapply(dt, function(x) x[-conditioned.site]))/sapply(dt, "[[", conditioned.site)) %>% sweep(2, conditioned.site.lambda[-1], "+")
    
    if(nlocs == 2){
      omegas = t(omegas)
    }
    
    summ = apply(omegas, MARGIN = 1, function(x){t(x) %*% inversePsi %*% (x)}) %>% sum
    length(dt)*detPsi+ summ # LL
  }
  
  # negative likelihood function
  ngll = function(pars){
    
    print(pars)
    if( pars[2] <= 0) return (2^100)
    
    if(any((pars[1] + pars[3]*temporal_covar) <= 0))return (2^100)
    
    print("hi")
    LL <- foreach(i= 1:length(exceedances), .combine = 'c') %dopar% {
      
      my.vario <- function(h){
        if(stat){
          alpha = pars[1]
          beta = pars[2] 
        }else{
          alpha = pars[1] + pars[3]*temporal_covar[i]
          beta = pars[2] 
        }
        
        
        res = rep(NA, length(h))
        res[h == 0] = 0
        res[h != 0] = alpha*(1 - ((((h[h != 0] /beta)^nu) * besselK(x = (h[h != 0] / beta), nu = nu))/(gamma(nu)*2^(nu - 1))))
        res
        
      }
      HRD_ll(lcs = exceedances_locs[i][[1]], dt = exceedances[i], vr = my.vario, zeta = 1, gamma = 1, eta = 1)
    }
    LL = sum(LL)
    print(LL)
    LL
  }

  temp_covars  = read_csv('data/processed/obs_data_mintp_winter.csv') %>% 
    dplyr::select(year, month, loess_temp_anom) %>%
    unique()
  
  data_for_rpareto = readRDS(paste0("data/processed/data_for_rpareto/true/data_for_rpareto_", marg_mod, "_thresh_qnt_", thresh_qnt))
  .GlobalEnv$exceedances <- data_for_rpareto$exceedances
  .GlobalEnv$exceedances_locs <- data_for_rpareto$exceedances_locs
  .GlobalEnv$temporal_covar <- data_for_rpareto$extreme_dates %>%
    mutate(year = lubridate::year(date), month = lubridate::month(date)) %>%
    left_join(temp_covars) %>%
    pull(loess_temp_anom)

  .GlobalEnv$HRD_ll <- HRD_ll
  .GlobalEnv$nu <- nu
  
  .GlobalEnv$stat <- FALSE
  cl <<- parallel::makeForkCluster(num_clusters)
  doParallel::registerDoParallel(cl)
  fit <<- optim(par = c(readRDS(paste0("output/rparp_", 'E2', "_thresh_qnt_", 0.96))$par , -0.2), fn = ngll ,hessian = T) # Variance + range constant
  parallel::stopCluster(cl)
  fit %>% saveRDS(paste0("output/test_nonstat_rparp_", marg_mod, "_thresh_qnt_", thresh_qnt))
}


job::job({fit_rparp_true("E2", 18, 0.1, 0.96)})
# job::job({fit_rparp_true("E1", 1, 0.1, 0.96)})
# job::job({fit_rparp_true("E2", 1, 0.1, 0.96)})
# job::job({fit_rparp_true("E2", 1, 0.1, 0.95)})
# job::job({fit_rparp_true("E2", 1, 0.1, 0.97)})


fit = readRDS(paste0("output/test_nonstat_rparp_", 'E2', "_thresh_qnt_", 0.96))
fisher_info<-solve(fit$hessian)
prop_sigma<-sqrt(diag(fisher_info))
fit$par+1.96*prop_sigma
fit$par-1.96*prop_sigma



