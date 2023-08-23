# script for checking that resampling code works correctly


obs_data = readRDS("output/data_unif_marg")

bts_to_test= sample(seq(200), 25)

for(bts in bts_to_test){
  
  bootstrapped_data = readRDS(paste0('output/bootstrapped_data_sets/thresh_qnt_',0.96,'_model_','m22',"_bootstrap_",bts))
  this_bts = readRDS("data/processed/bootstrap_data/bts_dates")[[bts]]
  dates = tibble(orig_dates = sort(obs_data$date %>% unique),
                 sampled_dates = this_bts)
  
  
  res = 0
  
  
  for(k in seq(20)){
    
    date_to_check = sample_n(dates, 1)
    
    this_sum = bootstrapped_data %>%
      filter(date == date_to_check$orig_dates) %>%
      dplyr::select(Station, unif) %>%
      left_join(obs_data %>%
                  filter(date == date_to_check$sampled_dates) %>%
                  dplyr::select(Station, unif = unif_m22), 
                by = 'Station') %>%
      # filter(unif.x < 0.7, unif.y < 0.7, unif.y > 0) %>%
      filter(unif.y > 0) %>%
      mutate(diff = unif.x - unif.y) %>%
      pull(diff) %>%
      sum()
    
    
    if(this_sum != 0){
      print(this_sum)
      print("!!!!!!!!!!!!!!!!!! FAIL !!!!!!!!!!!!!!!!!!")
      print(paste0("Bootstrap:", bts))
      print(date_to_check)
      skdbasjhg()
    }else{
      paste0("passed :) bootstap: ", bts, ", date: ", date_to_check$sampled_dates) %>%
        print()
    }
  }
}