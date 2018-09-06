setwd("~/R_stuff/Ku_chips/")
library(data.table)
library(ggplot2)
hist_data <- fread("For_histogram.tsv", sep = "\t", dec = ",")
#hist_data$Cq <- sub(",", ".", hist_data$Cq)
#hist_data$`quantity/100` <- sub(",", ".", hist_data$`quantity/100`)

hist_data$Cq <- as.numeric(hist_data$Cq)
hist_data$`quantity/100` <- as.numeric(hist_data$`quantity/100`)


ggplot(hist_data, aes(x=Target, `quantity/100`))+
  geom_bar(stat = "identity", colour = "black", fill = "green", alpha = "0.6")+
  labs(title = "Top15 TE families",
       x = "TE family name",
       y = "Occupied nucleotides")+
  theme_grey() %+replace% 
  theme(panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "grey20"),
        panel.grid.major = element_line(colour = "grey92"), 
        panel.grid.minor = element_line(colour = "grey92", size = 0.25),
        strip.background = element_rect(fill = "grey85", colour = "grey20"),
        legend.key = element_rect(fill = "white", colour = NA),
        axis.text.x = element_text(angle = 45, hjust = 0.5),
        complete = TRUE)
