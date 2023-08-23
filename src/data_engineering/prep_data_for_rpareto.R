# get extreme spatial events and standardise them such that they can be modelled with an r pareto process

rm(list = ls())
library(tidyverse)


prep_rpareto_true = function(m, thresh_qnt){
  
  source('src/models/marginal/gpd_models.R')
  
  conditioning_site = 'shannon_airport' # longest running site
  
  obs_data = read_csv("data/processed/obs_data_mintp_winter.csv") %>%
    left_join(read_csv(paste0("data/processed/thresh_exceedance_lambda_mintp.csv")) ) %>%
    left_join(readRDS(paste0("output/quant_models_mintp.csv")))
  
  row_ids = seq(nrow(obs_data))
  
  unif_res = rep(NA, length(row_ids))
  for(i in row_ids){
    unif_res[i] = obs_data$temp_to_tau[i][[1]](obs_data$temp[i])
  }
  
  
  obs_data$unif = unif_res

  # pull date which have 'conditioning' site
  dates_to_keep = obs_data %>%
    group_by(date) %>%
    summarise(check_obs_for_s = conditioning_site %in% Station) %>%
    filter(check_obs_for_s) %>%
    pull(date)
  
  obs_data = obs_data %>% filter(date %in% dates_to_keep)
  
  
  obs_data_standardised = rlang::duplicate(obs_data)
  

  assign(paste0('pars_', m),readRDS(paste0('output/bts_param_est/thresh_qnt_',thresh_qnt,"_model_", m)))
  pred <<- get(paste0('pars_', m)) %>% predict_gpd(model = m, data = obs_data_standardised)
  obs_data_standardised$scale = pred$scale
  obs_data_standardised$shape = pred$shape
  
  
  if(thresh_qnt == 0.95){
    obs_data_standardised$scales_logged = obs_data_standardised$scales_95_logged
    obs_data_standardised$threshold = obs_data_standardised$threshold_tn_l_o_95_w_coast
    obs_data_standardised$threshold_exceedence =  obs_data_standardised$thresh_exceedance_95_coast
    
  }else if(thresh_qnt == 0.96){
    obs_data_standardised$scales_logged = obs_data_standardised$scales_96_logged
    obs_data_standardised$threshold = obs_data_standardised$threshold_tn_l_o_96_w_coast
    obs_data_standardised$threshold_exceedence =  obs_data_standardised$thresh_exceedance_96_coast
    
  }else if(thresh_qnt == 0.97){
    obs_data_standardised$scales_logged = obs_data_standardised$scales_97_logged
    obs_data_standardised$threshold = obs_data_standardised$threshold_tn_l_o_97_w_coast
    obs_data_standardised$threshold_exceedence =  obs_data_standardised$thresh_exceedance_97_coast

  }else{
    adssfgh()
  }
  
  ex_ind = obs_data_standardised$temp > obs_data_standardised$threshold
  obs_data_standardised$unif[ex_ind] = 1 - obs_data_standardised$threshold_exceedence[ex_ind]*(1-evd::pgpd((obs_data_standardised$temp[ex_ind] - obs_data_standardised$threshold[ex_ind]), loc = 0,scale = obs_data_standardised$scale[ex_ind], shape = obs_data_standardised$shape[1]))

  obs_data_standardised$unif[obs_data_standardised$unif<0] = 0
  obs_data_standardised$frechet_marg = -1/log(obs_data_standardised$unif)
  
  
  # -------- MULTIVARIATE EVENTS
  # get cost of each event
  extreme_dates = obs_data_standardised %>%
    dplyr::group_by(date) %>%
    dplyr::summarise(cost = mean(frechet_marg)) %>%
    ungroup() %>%
    arrange(desc(cost))
  
  # get cost threshold
  threshold = quantile(extreme_dates$cost, 0.8) %>% as.numeric
  
  
  threshold %>% saveRDS(paste0("output/r_thres_", m, "_thresh_quant_",thresh_qnt))
  
  # get extreme  dates
  extreme_dates = extreme_dates %>%
    filter(cost > threshold) %>%
    arrange(desc(cost))
  

  
  # temporally decluster events
  my_data_tmp = extreme_dates %>% arrange(date)
  
  min_dif = my_data_tmp %>% arrange(date) %>% pull(date) %>% diff() %>% min
  
  print(my_data_tmp)
  print("declustering events")
  
  while(min_dif<7){
    i<-1
    while(i < nrow(my_data_tmp)){
      if(diff(c(my_data_tmp[i,]$date, my_data_tmp[i+1,]$date))<7){
        # pick the largest
        if(my_data_tmp[i,]$cost > my_data_tmp[i+1,]$cost){
          my_data_tmp <- my_data_tmp[my_data_tmp$date != my_data_tmp[i+1,]$date,]
        }else{
          my_data_tmp <- my_data_tmp[my_data_tmp$date != my_data_tmp[i,]$date,]
        }
      }
      i<-i+1
    }
    min_dif <- my_data_tmp$date %>% diff() %>% min %>% abs()
  }
  extreme_dates = my_data_tmp
  
  # tibble of extreme events
  obs_data_standardised = obs_data_standardised %>%
    filter(date %in% extreme_dates$date) %>%
    arrange(date)
  
  print("getting data in correct format")
  
  # get observations in list
  exceedances = obs_data_standardised %>%
    group_by(date) %>%
    group_map(~{
      
      # our "conditional" site at the top of the list
      if(conditioning_site %in% .x$Station){
        c(.x %>% filter(Station == conditioning_site) %>% pull(frechet_marg),
          .x %>% filter(Station != conditioning_site) %>% pull(frechet_marg))
      }
    })
  
  # remove observatiosn that dont have valentia_observatory
  to_remove = exceedances %>% map(is.null) %>% unlist %>% which()
  
  exceedances_locs = obs_data_standardised %>%
    group_by(date) %>%
    group_map(~{
      if(conditioning_site %in% .x$Station){
        rbind(.x %>% filter(Station == conditioning_site) %>% dplyr::select(Long.projected, Lat.projected) %>% as.matrix(),
              .x %>% filter(Station != conditioning_site) %>% dplyr::select(Long.projected, Lat.projected ) %>% as.matrix())
      }
    })
  
  if(!is.na(to_remove[1])){
    exceedances = exceedances[-to_remove]
    exceedances_locs = exceedances_locs[-to_remove]
  }
  
  dates_to_rem = extreme_dates %>% arrange(date) %>% .[to_remove,] %>% pull(date)
  
  extreme_dates = extreme_dates %>% arrange(date)
  
  extreme_dates = extreme_dates[(exceedances_locs %>% map(nrow) %>% unlist) > 1,]
  exceedances = exceedances[(exceedances_locs %>% map(nrow) %>% unlist) > 1]
  exceedances_locs = exceedances_locs[(exceedances_locs %>% map(nrow) %>% unlist) > 1]
  
  list(exceedances = exceedances,
       exceedances_locs = exceedances_locs,
       thresh = threshold,
       extreme_dates = extreme_dates) %>%
    saveRDS(paste0("data/processed/data_for_rpareto/true/data_for_rpareto_", m, "_thresh_qnt_", thresh_qnt))
  
}

job::job({prep_rpareto_true("E1", 0.96)})
job::job({prep_rpareto_true("E2", 0.96)})
job::job({prep_rpareto_true("E2", 0.97)})
job::job({prep_rpareto_true("E2", 0.95)})


# #### ---- bootstrap


prep_rpareto_bts = function(m, bts_range, thresh_qnt){
  
  
  conditioning_site = 'shannon_airport'
  
  
  
  for(bts_num in bts_range){
    print(bts_num)

    this_dat = readRDS(paste0('output/bootstrapped_data_sets/thresh_qnt_',thresh_qnt,'_model_',m,"_bootstrap_", bts_num)) %>%
      left_join(read_csv('data/processed/sites.csv'))

    extreme_dates = this_dat %>%
      dplyr::group_by(date) %>%
      dplyr::summarise(cost = mean(frechet)) %>%
      ungroup() %>%
      arrange(desc(cost))

    # get cost threshold
    threshold = quantile(extreme_dates$cost, 0.8) %>% as.numeric

    # get extreme  dates
    extreme_dates = extreme_dates %>%
      filter(cost > threshold) %>%
      arrange(desc(cost))


    # temporally decluster events
    my_data_tmp = extreme_dates %>% arrange(date)
    min_dif = my_data_tmp %>% arrange(date) %>% pull(date) %>% diff() %>% min

    print(my_data_tmp)
    print("declustering events")

    while(min_dif<7){
      i<-1
      while(i < nrow(my_data_tmp)){
        if(diff(c(my_data_tmp[i,]$date, my_data_tmp[i+1,]$date))<7){
          # pick the largest
          if(my_data_tmp[i,]$cost > my_data_tmp[i+1,]$cost){
            my_data_tmp <- my_data_tmp[my_data_tmp$date != my_data_tmp[i+1,]$date,]
          }else{
            my_data_tmp <- my_data_tmp[my_data_tmp$date != my_data_tmp[i,]$date,]
          }
        }
        i<-i+1
      }
      min_dif <- my_data_tmp$date %>% diff() %>% min %>% abs()
    }
    extreme_dates = my_data_tmp

    # tibble of extreme events
    this_dat = this_dat %>%
      filter(date %in% extreme_dates$date) %>%
      arrange(date)

    print("getting data in correct format")

    # get observations in list
    exceedances = this_dat %>%
      group_by(date) %>%
      group_map(~{

        # our "conditional" site at the top of the list
        if(conditioning_site %in% .x$Station){
          c(.x %>% filter(Station == conditioning_site) %>% pull(frechet),
            .x %>% filter(Station != conditioning_site) %>% pull(frechet))
        }
      })



    # remove observatiosn that dont have valentia_observatory
    to_remove = exceedances %>% map(is.null) %>% unlist %>% which()

    exceedances_locs = this_dat %>%
      group_by(date) %>%
      group_map(~{
        if(conditioning_site %in% .x$Station){
          rbind(.x %>% filter(Station == conditioning_site) %>% dplyr::select(Long.projected, Lat.projected) %>% as.matrix(),
                .x %>% filter(Station != conditioning_site) %>% dplyr::select(Long.projected, Lat.projected ) %>% as.matrix())
        }
      })


    if(!is.na(to_remove[1])){
      exceedances = exceedances[-to_remove]
      exceedances_locs = exceedances_locs[-to_remove]
      extreme_dates = extreme_dates[-to_remove,]
    }
    extreme_dates = extreme_dates %>% arrange(date)

    extreme_dates = extreme_dates[(exceedances_locs %>% map(nrow) %>% unlist) > 1,]
    exceedances = exceedances[(exceedances_locs %>% map(nrow) %>% unlist) > 1]
    exceedances_locs = exceedances_locs[(exceedances_locs %>% map(nrow) %>% unlist) > 1]

    list(exceedances = exceedances,
         exceedances_locs = exceedances_locs,
         thresh = threshold,
         extreme_dates = extreme_dates) %>%
      saveRDS(paste0('output/bootstrapped_rpareto_data/thresh_qnt_',thresh_qnt,'_bts_',bts_num ,'_model_',m))
  }
}



job::job({prep_rpareto_bts("E1", seq(1, 100), 0.96)})
job::job({prep_rpareto_bts("E1", seq(101, 200), 0.96)})
job::job({prep_rpareto_bts("E2", seq(1, 100), 0.96)})
job::job({prep_rpareto_bts("E2", seq(101, 200), 0.96)})
job::job({prep_rpareto_bts("E2", seq(1, 100), 0.95)})
job::job({prep_rpareto_bts("E2", seq(101, 200), 0.95)})
job::job({prep_rpareto_bts("E2", seq(1, 100), 0.97)})
job::job({prep_rpareto_bts("E2", seq(101, 200), 0.97)})

