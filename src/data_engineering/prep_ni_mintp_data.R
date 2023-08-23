rm(list=ls())
library(tidyverse)
library(vroom)


files = list.files() 

# only take qcv-1
files <- files[str_detect(files, "qcv-1")]   


NI_data = files %>%
  purrr::map(~{
    
    site_info = str_extract(.x, "202107.*qcv") %>%
      str_remove("202107_") %>%
      str_remove("_qcv") %>%
      str_split("_")
    
    line_to_start = grep("^data$", readLines(.x)) # csv data starts after "data"
    
    my_file_read = read.csv(.x, skip = line_to_start) %>%
      .[-nrow(.),] # remove last line
    
    
    my_file_read %>%
      dplyr::select(ob_end_time,
                    src_id,
                    min_air_temp, 
                    min_air_temp_j,
                    min_air_temp_q) %>%
      mutate(Station = site_info[[1]][3])
  }) %>%
  plyr::rbind.fill()%>%
  as_tibble() %>%
  left_join(read.csv("../NI_meta.csv"))


# create indicator vars
inds = NI_data$min_air_temp_q %>% as.character()

# if 0, has had no processing, can ignore
NI_data$ind_1 = inds %>% substr(start = nchar(inds), stop = nchar(inds)) %>%
  as.numeric()

# should be 0 or NA
NI_data$ind_2 = inds %>% substr(start = nchar(inds)-1, stop = nchar(inds)-1) %>%
  as.numeric()

# remove if ind > 0, its either suspect or estimated (https://dap.ceda.ac.uk/badc/ukmo-midas/metadata/doc/QC_J_flags.html)
NI_data$ind_3 = inds %>% substr(start = nchar(inds)-2, stop = nchar(inds)-2) %>%
  as.numeric()

NI_data$ind_4 = inds %>% substr(start = nchar(inds)-3, stop = nchar(inds)-3) %>%
  as.numeric()


# where there is more than one obs in a day, take the maximum
NI_data = NI_data %>% 
  filter((ind_1 > 0) & !is.na(ind_1)) %>%
  filter((ind_2 == 0) | is.na(ind_2)) %>%
  filter((ind_3 == 0) | is.na(ind_3)) %>%
    dplyr::select(ob_end_time,
                  src_id,
                  min_air_temp,
                  Station,
                  Lat,
                  Long) %>%
  group_by(Station) %>%
  group_map(~{
    .x %>%
      mutate(date = lubridate::as_date(ob_end_time)) %>%
      group_by(date) %>%
      slice(which.min(min_air_temp)) %>%
      ungroup()
  },.keep = T)%>%
  # save as one big tibble
  plyr::rbind.fill()%>%
  as_tibble() %>%
  dplyr::select(date, Station, Long, Lat, id = src_id, mintp = min_air_temp)

# save as csv
NI_data %>%
  write_csv("~/data/processed/NI_mintp_data.csv")
