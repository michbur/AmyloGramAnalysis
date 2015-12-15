library(ggplot2)
library(grid)
library(gridExtra)

size_mod <- -5

my_theme <- theme(plot.background=element_rect(fill = "transparent",
                                               colour = "transparent"),
                  panel.grid.major = element_line(colour="grey", linetype = "dashed", size = 0.5),
                  panel.grid.major = element_line(colour="lightgrey", linetype = "dashed", size = 0.5),
                  panel.background = element_rect(fill = "transparent",colour = "black"),
                  legend.background = element_rect(fill = "NA"),
                  legend.position = "bottom",
                  axis.text = element_text(size=13 + size_mod),
                  axis.title.x = element_text(size=16 + size_mod, vjust = -1), 
                  axis.title.y = element_text(size=16 + size_mod, vjust = 1),
                  strip.text = element_text(size=17 + size_mod, face = "bold"),
                  legend.text = element_text(size=13 + size_mod), 
                  legend.title = element_text(size=17 + size_mod),
                  plot.title = element_text(size=20 + size_mod),
                  strip.background = element_rect(fill = "NA", colour = "NA"))