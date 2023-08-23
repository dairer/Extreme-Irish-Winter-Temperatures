#------ get full boostrap, run after bootstrap_by_date_blocks.R
gc()
run_bts = function(models_to_run, bts_range, thresh_qnt){

  library(tidyverse)
  source("src/models/marginal/gpd_models.R")

  obs_data = read_csv("data/processed/obs_data_mintp_winter.csv") %>%
    left_join(read_csv(paste0("data/processed/thresh_exceedance_lambda_mintp.csv")) ) %>%
    left_join(readRDS(paste0("output/quant_models_mintp.csv")))


  # transform all data to uniform using bulk, later we transorm tail using gpd and append using thresold exceedence
  row_ids = seq(nrow(obs_data))
  unif_res = rep(NA, length(row_ids))
  for(i in row_ids){
    unif_res[i] = obs_data$temp_to_tau[i][[1]](obs_data$temp[i])
  }

  obs_data$unif = unif_res
  obs_data = obs_data %>%
    dplyr::select(
      tau_to_temp,
      date, 
      Station,
      temp, 
      Long, 
      Lat,
      year,
      dist_sea_logged,
      month,
      threshold_tn_l_o_95,
      threshold_tn_l_o_96,
      threshold_tn_l_o_97,
      thresh_exceedance_95,
      thresh_exceedance_96,
      thresh_exceedance_97,
      threshold_tn_l_o_95_w_coast,
      threshold_tn_l_o_96_w_coast,
      threshold_tn_l_o_97_w_coast,
      thresh_exceedance_95_coast,
      thresh_exceedance_96_coast,
      thresh_exceedance_97_coast,
      scales_95_logged,
      scales_96_logged,
      scales_97_logged,
      temp_anom,
      yearly_residuals,
      loess_temp_anom,
      residuals,
      week,
      unif)
  

  if(thresh_qnt == 0.95){
    obs_data$threshold =  obs_data$threshold_tn_l_o_95_w_coast
    obs_data$threshold_exceedence =  obs_data$thresh_exceedance_95_coast
    obs_data$scales_logged =  obs_data$scales_95_logged
  }else if(thresh_qnt == 0.96){
    obs_data$threshold =  obs_data$threshold_tn_l_o_96_w_coast
    obs_data$threshold_exceedence =  obs_data$thresh_exceedance_96_coast
    obs_data$scales_logged =  obs_data$scales_96_logged
  }else if(thresh_qnt == 0.97){
    obs_data$threshold =  obs_data$threshold_tn_l_o_97_w_coast
    obs_data$threshold_exceedence =  obs_data$thresh_exceedance_97_coast
    obs_data$scales_logged =  obs_data$scales_97_logged
  }else{
    asfdgsafadafdgf() # for debugging
  }
  

  extreme = obs_data %>%
    mutate(excess = temp - threshold) %>%
    filter(excess>0)
  
  ex_ind = (obs_data$temp)>(obs_data$threshold)
  
  # --- transform data to uniform
  for(m in models_to_run){
    
    print(paste0("Running model ", m))
    
    fit = fit_gpd(data = extreme, model = m)
    fit %>% saveRDS(paste0('output/bts_param_est/thresh_qnt_',thresh_qnt,"_model_", m)) # fit on the orig data
    
    pred = fit %>% predict_gpd(data = obs_data, model = m)
    
    obs_data[[paste0("scale_", m)]] = pred$scale
    obs_data[[paste0("shape_", m)]] = pred$shape
    obs_data[[paste0("unif_", m)]] = obs_data$unif
    
    obs_data[[paste0("unif_", m)]][ex_ind] = 1 - obs_data$threshold_exceedence[ex_ind]*(1-evd::pgpd(q = (obs_data$temp[ex_ind]-obs_data$threshold[ex_ind]),
                                                                                                         loc = 0,
                                                                                                         scale = obs_data[[paste0("scale_", m)]][ex_ind],
                                                                                                         shape = obs_data[[paste0("shape_", m)]][1]))
  }

  obs_data = obs_data %>%
    group_by(date) %>%
    arrange(Station, .by_group = TRUE)
  
  # ordering rows
  for(m in models_to_run){
    print(paste0("Running model ", m))
    assign(paste0("unifs_",m), obs_data[c('date', 'Station', paste0("unif_", m))] %>%
             ungroup() %>%
             pivot_wider(names_from = 'Station', values_from = paste0("unif_", m)))
    
    assign(paste0("unifs_",m), get(paste0("unifs_",m))[,order(colnames(get(paste0("unifs_",m))))])
  }
  non_neg_vec = !(get(paste0("unifs_",m)) %>% pivot_longer(-c(date)) %>% pull(value) %>% is.na())

  
  
  for(i in bts_range){
    print(i)
    set.seed(i)
    
    # ---- transform from unif to orig

    
    # read in bootstrap dates
    this_bts = readRDS("data/processed/bootstrap_data/bts_dates")[[i]]

    # taking bulk from any model, all equal
    this_data = tibble(date = this_bts) %>% left_join(get(paste0("unifs_",models_to_run[1]))) # just take first model, bulks equal across models
    this_data = this_data[,order(colnames(this_data))]
    this_data = this_data %>% pivot_longer(-c(date)) %>% pull(value)
    this_bts = rlang::duplicate(obs_data, shallow = FALSE)
    
    this_bts$new_samp = this_data[non_neg_vec]
    this_bts_unifs = this_data[non_neg_vec]
    orig_marg_res = rep(NA, length(this_bts_unifs))
    row_ids = seq(length(this_bts_unifs))
    
    # transforming bulk from unif to orig
    for(k in row_ids){
      if(!is.na(this_bts_unifs[k])){
        orig_marg_res[k] = obs_data$tau_to_temp[k][[1]](this_bts_unifs[k])
      }
    }
    
    temp <<- orig_marg_res
    

    # --- for each marginal model transform tail to orig
    for(m in models_to_run){
      
      
      print(m)
      
      this_bts = readRDS("data/processed/bootstrap_data/bts_dates")[[i]]
      this_data = tibble(date = this_bts) %>% left_join(get(paste0("unifs_",m)))
      this_data = this_data[,order(colnames(this_data))]
      this_data = this_data %>% pivot_longer(-c(date)) %>% pull(value)
      
      this_bts = rlang::duplicate(obs_data, shallow = FALSE)
      this_bts$new_samp = this_data[non_neg_vec]
      
      
      
      this_bts = this_bts %>%
        ungroup() %>%
        dplyr::select(-unif) %>%
        rename(unif = new_samp)
      
      
      # ---- transformed bulk from unif to orig
      this_bts$temp = temp
      
      
      #  check ballincurrig_(peafield)
      this_bts = this_bts %>% drop_na() # should there be no NA?
      this_bts[this_bts$unif<0,]$unif = 0
      
      # doesnt change from one marginal model to the next
      ex_dat_ind = this_bts$unif>(1-this_bts$threshold_exceedence)
      # ex_dat_ind = this_bts$unif>0.96
      
      this_bts$temp[ex_dat_ind] = evd::qgpd(p = (1 + (this_bts$unif[ex_dat_ind]-1)/this_bts$threshold_exceedence[ex_dat_ind]),
                                            loc=this_bts$threshold[ex_dat_ind],
                                            scale=this_bts[[paste0('scale_', m)]][ex_dat_ind],
                                            shape=this_bts[[paste0('shape_', m)]][1])
      
      this_bts = this_bts %>%
        dplyr::select(Station, date, unif, temp) %>%
        mutate(frechet = -1/log(unif)) 
      
      this_bts %>% saveRDS(paste0('output/bootstrapped_data_sets/thresh_qnt_',thresh_qnt,'_model_',m,"_bootstrap_", i))
    }
  }
}


# job::job({run_bts("E2", bts_range = seq(1, 200), thresh_qnt = 0.96)})
# job::job({run_bts("E2", bts_range = seq(1, 200), thresh_qnt = 0.95)})
# job::job({run_bts("E2", bts_range = seq(1, 200), thresh_qnt = 0.97)})
# job::job({run_bts("E1", bts_range = seq(1, 200), thresh_qnt = 0.96)})
