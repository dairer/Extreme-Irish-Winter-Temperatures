# this script calculates the average length of an extreme event, classified as exceeding the threshold for a consecutive week

obs_data = read_csv("data/processed/obs_data_mintp_winter.csv")

extreme_events_per_site = obs_data %>%
  filter(temp>threshold_tn_l_o_96_w_coast) %>%
  group_by(Station) %>%
  group_map(~{
    
    x = .x %>%
      pull(date) %>%
      sort() %>%
      diff() 

    overall_count = 0
    current_count = 0
    for(i in x){
      
      if(i <= 7){
        current_count = current_count + 1
        if(current_count > overall_count) overall_count = current_count
      }
      else current_count = 0
      
    }
    
    tibble(Station = .x$Station[1], overall_count = overall_count)
  }, .keep = T) %>%
  plyr::rbind.fill() %>%
  as_tibble()



extreme_events_per_site %>%
  ggplot()+
  geom_histogram(aes(overall_count), col = 'black', alpha = 0.5)+
  theme_minimal(10)+
  labs(x = "Length of extreme event",
       y = "Count")+
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

ggsave("output/figs/num_and_len_extreme_events.pdf", width = 4, height = 2.35)



