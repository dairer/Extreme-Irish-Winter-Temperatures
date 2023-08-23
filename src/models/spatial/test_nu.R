# this function fits r pareto process over range of shape parameter for the matern variogram
# used to select shape vale that gave the best fit
rm(list=ls())
library(tidyverse)
library(doParallel)
library(foreach)

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

marg_mod = 'E2'

for(nu in seq(0.09, 0.2, length.out = 10)){
    
  print(nu)
  data_for_rpareto = readRDS(paste0("data/processed/data_for_rpareto/true/data_for_rpareto_", marg_mod))
  .GlobalEnv$exceedances <- data_for_rpareto$exceedances
  .GlobalEnv$exceedances_locs <- data_for_rpareto$exceedances_locs
  .GlobalEnv$HRD_ll <- HRD_ll
  .GlobalEnv$nu <- nu
  
  .GlobalEnv$stat <- TRUE
  cl <<- parallel::makeForkCluster(5)
  doParallel::registerDoParallel(cl)
  fit <<- optim(par = c(1,1), fn = ngll ,hessian = F, method = 'Nelder-Mead' ) # Variance + range constant
  parallel::stopCluster(cl)
  
  fit %>% saveRDS(paste0("output/nu_",nu,"_rparp_", marg_mod))
}
  

files = list.files("output/")
files = files[grepl(paste0("_rparp_", marg_mod), files) & grepl("nu_", files)]


vals = c()
for(f in files){

  vals = rbind(vals,
               c(str_remove(f, "nu_") %>% str_remove(paste0("_rparp_", marg_mod)) %>% as.numeric(),
                 readRDS(paste0("output/",f))$value))
}

# plot results
readRDS(paste0("output/","nu_0.1_rparp_E2"))
vals %>%
  as_tibble() %>%
  ggplot()+
  geom_line(aes(V1, V2))+
  theme_minimal()+
  labs(x = expression(nu), y = "log likelihood")
