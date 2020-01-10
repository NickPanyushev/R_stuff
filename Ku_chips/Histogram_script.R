setwd("~/R_stuff/Ku_chips/")
library(data.table)
library(ggplot2)
hist_data <- fread("For_histogram.tsv", sep = "\t", dec = ",")

hist_data$Cq <- as.numeric(hist_data$Cq)
hist_data$`quantity/100` <- as.numeric(hist_data$`quantity/100`)
hist_data$Target <- as.factor(hist_data$Target)
hist_data$Sample <- as.factor(hist_data$Sample)
hist_data$Cq <- hist_data$`ChIP number` <- NULL
hist_data[, SD := sd(`quantity/100`), by = .(Sample, Target)]
hist_data[, MEAN := mean(`quantity/100`), by = .(Sample, Target)]

mean_data <- hist_data[, sd(`quantity/100`), by = .(Sample, Target)]
mean_data$Means <- hist_data[, mean(`quantity/100`), by = .(Sample, Target)]$V1
mean_data$CI <- hist_data[, qnorm(0.975)*sd(`quantity/100`)/sqrt(.N), by = .(Sample, Target)]$V1


ggplot(mean_data, aes(x = Target, y = Means))+   
  geom_bar(aes(fill = Sample), position = "dodge", stat="identity")+
  
  geom_point(data = hist_data, size = 1.9, position = position_jitterdodge(jitter.width = 0.3,
                                                                          dodge.width = 0.9),  
              aes(x = Target, y = `quantity/100`, 
                  shape = Sample, group = Sample),
    colour="black")+
  
  geom_errorbar(aes(ymin=Means-CI, ymax=Means+CI, group = Sample), 
                colour="black", width=0.1, position=position_dodge(0.9))+
  
  labs(title = "ChIP data",
       x = "Target DNA",
       y = "IP/Input ratio") +
  theme_bw()

