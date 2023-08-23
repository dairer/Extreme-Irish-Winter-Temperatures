# code here is adapted from R package mvPot: Multivariate Peaks-over-Threshold Modelling for Spatial Extreme Events, author: de Fondeville
# https://cran.r-project.org/web/packages/mvPot/index.html


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
    # for some reason when we have 2 sites the matrix is transposed...
    
    if(nlocs == 2){
      omegas = t(omegas)
    }
    
    summ = apply(omegas, MARGIN = 1, function(x){t(x) %*% inversePsi %*% (x)}) %>% sum
    length(dt)*detPsi+ summ # LL
  }
  
  # negative likelihood function
  ngll = function(pars){
    
    print(pars)
    if( pars[2] < 0) return (2^100)
    if( pars[1] < 0) return (2^100)
    
    
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
        
        if(alpha <= 0) return (2^100)
         
        
        if( beta <= 0) return (2^100)

        # nu = 0.2
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

  data_for_rpareto = readRDS(paste0("data/processed/data_for_rpareto/true/data_for_rpareto_", marg_mod, "_thresh_qnt_", thresh_qnt))
  .GlobalEnv$exceedances <- data_for_rpareto$exceedances
  .GlobalEnv$exceedances_locs <- data_for_rpareto$exceedances_locs
  .GlobalEnv$HRD_ll <- HRD_ll
  .GlobalEnv$nu <- nu
  
  .GlobalEnv$stat <- TRUE
  cl <<- parallel::makeForkCluster(num_clusters)
  doParallel::registerDoParallel(cl)
  fit <<- optim(par = c(1,1), fn = ngll ,hessian = F, method = 'Nelder-Mead' ) # Variance + range constant
  parallel::stopCluster(cl)
  
  # method = 'L-BFGS-B'
  fit %>% saveRDS(paste0("output/rparp_", marg_mod, "_thresh_qnt_", thresh_qnt))
}

job::job({fit_rparp_true("E1", 1, 0.1, 0.96)})

job::job({fit_rparp_true("E2", 1, 0.1, 0.96)})
job::job({fit_rparp_true("E2", 1, 0.1, 0.95)})
job::job({fit_rparp_true("E2", 1, 0.1, 0.97)})

readRDS(paste0("output/rparp_", 'E2', "_thresh_qnt_", 0.96))

variogram_model_score <- function(h, vario_params, nu){
  alpha <-  (vario_params[1])
  beta <- (vario_params[2])
  # nu = 0.2
  res = rep(NA, length(h))
  res[h == 0] = 0
  res[h != 0] = alpha*(1 - ((((h[h != 0] /beta)^nu) * besselK(x = (h[h != 0] / beta), nu = nu))/(gamma(nu)*2^(nu - 1))))
  res
}

distances = seq(0.05,4, length.out = 100)


rbind(read_csv(paste0("output/new_chi_obs_bts_25_bins", 0.8,".csv"),
               col_names = c("i", "bts", "dist", "chi")) %>% mutate(p = "p = 0.8"),
      read_csv(paste0("output/new_chi_obs_bts_25_bins", 0.85,".csv"),
               col_names = c("i", "bts", "dist", "chi")) %>% mutate(p = "p = 0.85"),
      read_csv(paste0("output/new_chi_obs_bts_25_bins", 0.9,".csv"),
               col_names = c("i", "bts", "dist", "chi")) %>% mutate(p = "p = 0.9"),
      read_csv(paste0("output/new_chi_obs_bts_25_bins", 0.95,".csv"),
               col_names = c("i", "bts", "dist", "chi")) %>% mutate(p = "p = 0.95"))%>%
  group_by(dist, p) %>%
  summarise(upper = quantile(chi, 0.025, na.rm = T),
            lower = quantile(chi, 0.975, na.rm = T),
            chi = mean(chi, na.rm = T)) %>%
  mutate(d = "Observational") %>%
  ggplot()+
  geom_segment(aes(x = dist*100, xend = dist*100, y = lower, yend = upper, group = d, col = d))+
  # geom_point(aes(x = dist*100, y=chi, shape = d, group = d, col = d))+
  geom_line(data = tibble(distances,
                          mod_4b = 2 * (1 - pnorm(sqrt(variogram_model_score(distances, c(3.750953, 461.775464), 0.1) / 2)))),
            aes(distances*100, mod_4b, col = 'mod_4a'))+
  facet_wrap(~p, nrow = 1)+
  theme_minimal(12)+
  labs(y = expression(chi),
       x = "Distance (km)")+
  theme(legend.position = 'none',
        axis.title.y = element_text(angle=0, vjust=.5))+
  ylim(0,1)
