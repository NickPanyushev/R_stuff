library(dplyr)
library(ggplot2)
setwd('~/Рабочий стол/')
table1 <- read.csv('/home/nickolay/Рабочий стол/бредфорд_апрель.csv', header = F)
OD <- (as.numeric(table1[1,]))
curve <- data.frame (Concentrations = c(0, 100, 250, 500, 1000, 2000, 4000, 8000), OD)
curve_trimmed <- data.frame(Concentrations = c(0, 100, 250, 500, 1000, 2000), OD = OD[seq(1,6)])

ggplot(curve_trimmed, aes(curve_trimmed$Concentrations, curve_trimmed$OD)) +
  geom_smooth(method='lm') + 
  geom_point() +
  theme_bw()

curve_trimmed.lm = lm(Concentrations ~ OD, data = curve_trimmed)
coeffs = coefficients(curve_trimmed.lm)

result_table <- table1[2:4,]
result_table <- coeffs[1] + coeffs[2]*result_table
for_WB <- data_frame(Samples = c('H1299', 'H1299+ACTN4', 'H1299+scsh',
                                 'H1299-ACTN4', 'H1299+p65', 'H1299+p65+ACTN4',
                                 'H1299+p65+scsh', 'H1299+p65-ACTN4'),
                     Concentrations = c(result_table[3,1:2], result_table[1,3],
                                        result_table[3,3:5], result_table[3,7],
                                        result_table[3,6]))

for_WB$Concentrations <- (as.numeric(for_WB$Concentrations)/1000)
for_WB$Sample_Volume <- (0.2/for_WB$Concentrations)*50
for_WB$Laemmly <- 10
for_WB$Ripa_Volume <- 50 - for_WB$Sample_Volume - for_WB$Laemmly

