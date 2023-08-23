# This script plots parameters over a series of high thresholds to investigate what threshold level is appropriate. 
library(tidyverse)

ex_dat = read_csv("data/processed/obs_data_mintp_winter.csv")

shape_res = c()
scale_res = c()
scale_se = c()
tres = c()

# itteratre over a series of high quantiles
for(q in seq(0.8, 0.999, by = 0.001)){
  print(q)

  # get exceedance data
  temporary = ex_dat %>%
    group_by(Station) %>%
    mutate(thresh = quantile(temp, q)) %>%
    mutate(excess = temp - thresh) %>%
    filter(excess > 0) %>%
    ungroup()

  # fit basic gev
  my_fit = evgam::evgam(list(excess ~ 1, ~1),
               data = temporary,
               family = 'gpd')

  # store parameter estimates along with uncertainties
  tres = c(tres, temporary$thresh[1])
  scale_res = rbind(scale_res, summary(my_fit)[[1]]$logscale[,1])
  scale_se = rbind(scale_se, 1.96*summary(my_fit)[[1]]$logscale[,2])
  shape_res = rbind(shape_res,
                    c(summary(my_fit)[[1]]$shape[,1],1.96*summary(my_fit)[[1]]$shape[,2] ))
}

# save estimates and uncertainties
list(threshold = tres,
     shape_res = shape_res,
     scale_res = scale_res,
     scale_se = scale_se) %>%
  saveRDS("output/thresh_stability_res")

# label for x axis of plot
bottom = gridtext::richtext_grob(
  text = 'Site-wise quantile'
)

# Plot paramters estimates 
plt = gridExtra::grid.arrange(tibble(thres = seq(0.8, 0.999, by = 0.001),
                               scale = readRDS("output/thresh_stability_res")$scale_res[,1],
                               scaleu = readRDS("output/thresh_stability_res")$scale_res[,1] + readRDS("output/thresh_stability_res")$scale_se[,1],
                               scalel = readRDS("output/thresh_stability_res")$scale_res[,1] - readRDS("output/thresh_stability_res")$scale_se[,1]) %>%
                          filter(thres < 0.99) %>%
                          ggplot()+
                          geom_vline(xintercept = 0.95, linetype = 'longdash', col = 'magenta')+
                          geom_vline(xintercept = 0.96, linetype = 'longdash', col = 'magenta')+
                          geom_vline(xintercept = 0.97, linetype = 'longdash', col = 'magenta')+
                          geom_ribbon(aes(thres, x = thres, ymin = scalel, ymax = scaleu), alpha = 0.3, fill = 'forestgreen')+
                          geom_line(aes(thres, scale))+
                          geom_line(aes(thres, scaleu), alpha = 0.5)+
                          geom_line(aes(thres, scalel), alpha = 0.5)+
                          scale_x_continuous(breaks = c(0.8, 0.85, 0.9, 0.95, 0.99))+
                          theme_minimal(12)+
                          theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
                                axis.title.x = element_blank())+
                                labs( y = expression(paste("ln", sigma))),
                        tibble(thres = seq(0.8, 0.999, by = 0.001),
                               shape = readRDS("output/thresh_stability_res")$shape_res[,1],
                               shapeu = readRDS("output/thresh_stability_res")$shape_res[,1] + readRDS("output/thresh_stability_res")$shape_res[,2],
                               shapel = readRDS("output/thresh_stability_res")$shape_res[,1] - readRDS("output/thresh_stability_res")$shape_res[,2]) %>%
                          filter(thres < 0.99) %>%
                          ggplot()+
                          geom_vline(xintercept = 0.95, linetype = 'longdash', col = 'magenta')+
                          geom_vline(xintercept = 0.96, linetype = 'longdash', col = 'magenta')+
                          geom_vline(xintercept = 0.97, linetype = 'longdash', col = 'magenta')+
                          geom_ribbon(aes(thres, x = thres, ymin = shapel, ymax = shapeu), alpha = 0.3, fill = 'forestgreen')+
                          # geom_vline(xintercept = 0.98)+
                          scale_x_continuous(breaks = c(0.8, 0.85, 0.9, 0.95, 0.99))+
                          geom_line(aes(thres, shape))+
                          geom_line(aes(thres, shapeu), alpha = 0.5)+
                          geom_line(aes(thres, shapel), alpha = 0.5)+
                          theme_minimal(12)+
                          theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
                                axis.title.x = element_blank())+
                          labs(y = expression(xi)), bottom = bottom, nrow = 1)

ggsave(plt, filename = "output/figs/thresh_stability.pdf", width = 8, height = 3)