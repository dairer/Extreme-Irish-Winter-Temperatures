# this script fits rpareto process to event bootstrapped data set
rm(list=ls())
library(tidyverse)
library(doParallel)
library(foreach)


fit_bts_mod = function(marg_mod, bts_seq, change_to_90 = FALSE, thresh_qnt){
  init_vals = readRDS(paste0("output/rparp_", marg_mod, "_thresh_qnt_", thresh_qnt))$par
  
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
    
    tryCatch({
      inversePsi = psi %>% solve()  # inverse of "covarianve" matrix
      
      detPsi = determinant(psi)$modulus[1] # calculates log det by default
      
      omegas = log(t(sapply(dt, function(x) x[-conditioned.site]))/sapply(dt, "[[", conditioned.site)) %>% sweep(2, conditioned.site.lambda[-1], "+")
      # for some reason when we have 2 sites the matrix is transposed...
      
      if(nlocs == 2){
        omegas = t(omegas)
      }
      
      summ = apply(omegas, MARGIN = 1, function(x){t(x) %*% inversePsi %*% (x)}) %>% sum
      length(dt)*detPsi+ summ # LL
      
    }, error=function(e){
      print("omg caugth error")
      return (2^30)
    })
    
  }
  
  ngll = function(pars){
    
    #print(pars)
    if( pars[2] < 0) return (2^100)
    if( pars[1] < 0) return (2^100)
    # if( pars[1] <= 10) return (2^100)
    
    LL <- foreach(i= 1:length(exceedances), .combine = 'c') %dopar% {
      my.vario <- function(h){
        if(stat){
          alpha = pars[1]
          beta = pars[2] 
        }else{
          alpha = pars[1] + pars[3]*temporal_covar[i]
          beta = pars[2] 
        }
        
        if(alpha <= 0){
          return (2^100)
        } 
        
        if( beta <= 0){
          return (2^100)
        } 
        
        nu = 0.1
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
  
  
  for(bs in bts_seq){
    print(bs)

    if(file.exists(paste0('output/bootstrapped_rpareto_data/thresh_qnt_',thresh_qnt,'_bts_',bs ,'_model_',marg_mod))){
      
      data_for_rpareto = readRDS(paste0('output/bootstrapped_rpareto_data/thresh_qnt_',thresh_qnt,'_bts_',bs ,'_model_',marg_mod))
      obs_to_rem = data_for_rpareto$exceedances %>% map(~{sum(.x<=0) <= 0}) %>% unlist  # remove obs with 0
      
      .GlobalEnv$exceedances <- data_for_rpareto$exceedances[obs_to_rem]
      .GlobalEnv$exceedances_locs <- data_for_rpareto$exceedances_locs[obs_to_rem]      
      .GlobalEnv$HRD_ll <- HRD_ll
      .GlobalEnv$stat <- TRUE
      
      tryCatch({
        set.seed(12341234)
        cl <- parallel::makeForkCluster(4)
        doParallel::registerDoParallel(cl)
        
        fit <- optimx::optimx(par = init_vals, fn = ngll ,hessian = F, method = 'L-BFGS-B') # Variance + range constant
        parallel::stopCluster(cl)
        
        
        tibble(bs, fit$p1, fit$p2) %>%
          write_csv(paste0("output/thresh_qnt_",thresh_qnt,"_bts_rpareto_fits_model_", marg_mod, "_L-BFGS-B.csv"), append = T)
        
      }, error=function(e){
        tibble(bs) %>%
          write_csv(paste0("output/thresh_qnt_",thresh_qnt,"_bts_rpareto_fits_failed_model_", marg_mod,".csv"), append = T)
      })
    }
  }
}


# job::job({fit_bts_mod(marg_mod = "E1", bts_seq = seq(200), thresh_qnt = 0.96)})
job::job({fit_bts_mod(marg_mod = "E2", bts_seq = seq(72,200), thresh_qnt = 0.96)})
job::job({fit_bts_mod(marg_mod = "E2", bts_seq = seq(85,200), thresh_qnt = 0.95)})
job::job({fit_bts_mod(marg_mod = "E2", bts_seq = seq(90,200), thresh_qnt = 0.97)})



# visalise fits
variogram_model_score <- function(h, vario_params){
  alpha <-  (vario_params[1])
  beta <- (vario_params[2])
  nu = 0.1
  res = rep(NA, length(h))
  res[h == 0] = 0
  res[h != 0] = alpha*(1 - ((((h[h != 0] /beta)^nu) * besselK(x = (h[h != 0] / beta), nu = nu))/(gamma(nu)*2^(nu - 1))))
  res
}

marg_mod = 'E2'
res = read_csv(paste0("output/bts_rpareto_fits_model_", marg_mod, "_L-BFGS-B.csv"),
         col_names = c('bts', 'alpha', 'beta'))


true = readRDS(paste0("output/rparp_", marg_mod))$par

distances = seq(0.005, 4, length.out = 100)

my_mat = c()
for(i in seq(nrow(res))){
  my_mat = rbind(my_mat, 2 * (1 - pnorm(sqrt(variogram_model_score(distances, c(res[i,]$alpha, res[i,]$beta)) / 2))))
}

tibble(distances,
       true =  2 * (1 - pnorm(sqrt(variogram_model_score(distances, c(true[1]  , true[2])) / 2))),
       lower = my_mat %>% apply(MARGIN = 2, FUN = quantile, 0.025),
       upper = my_mat %>% apply(MARGIN = 2, FUN = quantile, 0.975)
       ) %>%
  ggplot()+
  geom_line(aes(distances*100, true))+
  geom_ribbon(aes(distances*100, ymin = lower, ymax = upper), alpha = 0.25)+
  geom_line(aes(distances*100, upper), alpha = 0.25)+
  geom_line(aes(distances*100, lower), alpha = 0.25)+
  geom_point(data = read_csv(paste0("output/new_chi_obs_bts_25_bins", 0.9,".csv"),
                             col_names = c("i", "bts", "dist", "chi")), aes(dist*100, chi))+
  geom_segment(data = rbind(
    read_csv(paste0("output/new_chi_obs_bts_25_bins", 0.95,".csv"),
                                     col_names = c("i", "bts", "dist", "chi"))#,
    )  %>%
      group_by(dist) %>%
      summarise(upper = quantile(chi, 0.025, na.rm = T),
                lower = quantile(chi, 0.975, na.rm = T),
                chi = mean(chi, na.rm = T)),
    aes(x = dist*100, xend = dist*100, y = lower, yend = upper), col = 'red')+
  geom_segment(data = rbind(
                            read_csv(paste0("output/new_chi_obs_bts_25_bins", 0.9,".csv"),
                               col_names = c("i", "bts", "dist", "chi")))  %>%
                 group_by(dist) %>%
                 summarise(upper = quantile(chi, 0.025, na.rm = T),
                           lower = quantile(chi, 0.975, na.rm = T),
                           chi = mean(chi, na.rm = T)),
               aes(x = dist*100, xend = dist*100, y = lower, yend = upper), position = position_nudge(x = 1))+
  theme_minimal()