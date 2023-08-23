# created boostrap dates
rm(list = ls())
library(tidyverse)
library(vroom)
library(evgam)
setwd("~/chapter_6/")
var_modelled = 'mintp'
obs_data = read_csv("data/processed/obs_data_mintp_winter.csv") 

my_data = obs_data %>%
  group_by(date) %>%
  summarise(temp = max(temp)) %>%
  arrange(date)


indices <- which(diff(my_data$date) <= 15)[1]
while (!is.na(indices)) {
  if(my_data$temp[indices] > my_data$temp[indices+1]){
    my_data = my_data[-(indices+1),]
  }else{
    my_data = my_data[-(indices),]
  }
  indices <- which(diff(my_data$date) <= 15)[1]
}

clustered = obs_data %>%
  group_by(date) %>%
  summarise(temp = max(temp)) %>%
  arrange(date)


clustered$pos_blk = NA
clustered$pos_blk[which(clustered$date %in% my_data$date) + c(floor(diff(which(clustered$date %in% my_data$date))/2),0)] = clustered$date[which(clustered$date %in% my_data$date) + c(floor(diff(which(clustered$date %in% my_data$date))/2),0)] 

blk_pos = clustered %>%
  drop_na() %>%
  dplyr::select(pos_blk) %>%
  mutate(pos_blk = lubridate::as_date(as.numeric(pos_blk)))

# group data into its blocks
clustered$block_id = NA
block_id_stamp = 1
for(i in seq(nrow(clustered))){
  clustered[i,]$block_id = block_id_stamp
  if(!is.na(clustered[i,]$pos_blk)){
    block_id_stamp = block_id_stamp + 1
  }
}

clustered = clustered %>%
  dplyr::select(-pos_blk, -temp) 

obs_data = obs_data %>%
  dplyr::select(date, Station, temp) %>%
  left_join(clustered)


obs_data$block_id = as.numeric(obs_data$block_id)

obs_data = obs_data %>%
  arrange(block_id) %>%
  group_by(block_id) %>%
  arrange(Station,.by_group = T) %>%
  ungroup()


# This function determines if set A is a subset of set B
is.subset = function(A, B) all(A %in% B)
num.subset = function(A,B) (sum(A %in% B)/length(A))

obs_data %>%
  dplyr::select(block_id, Station) %>%
  unique() %>%
  group_by(block_id) %>%
  mutate(Station = list(Station)) %>%
  left_join(obs_data %>%
              dplyr::select(block_id, date) %>%
              unique() %>%
              group_by(block_id) %>%
              summarise(num_date = n())) %>%
  ungroup() %>%
  unique() %>%
  saveRDS("data/processed/bootstrap_data/block_sizes")


blocks_and_sites = readRDS("data/processed/bootstrap_data/block_sizes")

blocks_and_resample_choise = c()
loop_to = nrow(blocks_and_sites)

# find out which blocks can resample others
for(i in seq(loop_to)){
  
  print(i)  
  blocks_to_consider  = blocks_and_sites[blocks_and_sites[i,]$num_date <= blocks_and_sites$num_date,]
  num_group_subset = purrr::map(.x = blocks_to_consider$Station,
                         .f = ~num.subset(A = blocks_and_sites[i,]$Station[[1]],
                                          .x)) %>% unlist()
  
  resmple_choices = blocks_to_consider[(num_group_subset>0.9),]$block_id %>%
    unlist %>%
    unique() 
  
  blocks_and_resample_choise = rbind(blocks_and_resample_choise, c(blocks_and_sites[i,]$block_id, blocks_and_sites[i,]$num_date,list(resmple_choices)))
}

blocks_and_resample_choise = as_tibble(blocks_and_resample_choise) 
blocks_and_resample_choise$V1 = unlist(blocks_and_resample_choise$V1)
blocks_and_resample_choise$V2 = unlist(blocks_and_resample_choise$V2)

blocks_and_resample_choise = blocks_and_resample_choise %>%
  rename(block_id = V1,
         num_date = V2,
         resample_choices = V3) %>%
  arrange(block_id)

blocks_and_resample_choise %>% saveRDS('data/processed/bootstrap_data/blocks_and_resample_choise')
blocks_and_resample_choise =readRDS('data/processed/bootstrap_data/blocks_and_resample_choise')

blocks_and_resample_choise$resample_choices %>% map(length) %>% unlist() %>% hist

blocks_and_dates = obs_data %>%
  dplyr::select(date, block_id) %>%
  group_by(block_id) %>%
  unique() %>%
  summarise( date = list(date))


mn_dts = obs_data %>%
  group_by(block_id) %>%
  summarise(mean_date = mean(date)) 


blocks_and_resample_choise = readRDS('data/processed/bootstrap_data/blocks_and_resample_choise') %>%
  left_join(mn_dts)


mn_dts = obs_data %>%
  group_by(block_id) %>%
  summarise(mean_date = mean(date)) 

# ----- only keep obs_within 10 years
for(i in seq(nrow(blocks_and_resample_choise))){
  print(i)
  within_10_yrs = (abs(blocks_and_resample_choise[sort(blocks_and_resample_choise[i,]$resample_choices[[1]]),]$mean_date - blocks_and_resample_choise[i,]$mean_date)/365 <= 2.5)
  blocks_and_resample_choise[i,]$resample_choices = list(blocks_and_resample_choise[i,]$resample_choices[[1]][within_10_yrs])
}

blocks_and_resample_choise %>% saveRDS('data/processed/bootstrap_data/blocks_and_resample_choise')

bts_dates = c()
for(i in seq(200)){
  print(i)

  blocks_and_resample_choise = readRDS('data/processed/bootstrap_data/blocks_and_resample_choise') 

  blocks_and_resample_choise$new_sample = blocks_and_resample_choise$resample_choices %>%
    purrr::map(as.character) %>%
    purrr::map(sample, size = 1) %>%
    purrr::map(as.numeric) %>%
    unlist
  
  
  blocks_and_resample_choise = blocks_and_resample_choise %>%
    left_join(blocks_and_dates, c('new_sample' = 'block_id'))
  
  # new dates sampled:
  # this needs to be an ordered subset
  take_ordered_sample = function(x, n){
    len = length(x)
    ind_start = sample((len-n+1), size = 1)
    x[ind_start:(ind_start+n-1)]
  }
  
  bts_dates = c(bts_dates, list(purrr::map2(blocks_and_resample_choise$date,  blocks_and_resample_choise$num_date,  take_ordered_sample) %>% unlist %>% lubridate::as_date()))
}

saveRDS(bts_dates, paste0("data/processed/bootstrap_data/bts_dates"))

obs_data %>%
  dplyr::select(date, block_id) %>%
  unique() %>%
  write_csv("data/processed/date_and_block_id.csv")