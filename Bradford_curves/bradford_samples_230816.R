library(dplyr)
library(ggplot2)
setwd('~/Рабочий стол/')
# Построим калибровочную кривулину

table1 <- read.csv('/home/nickolay/Рабочий стол/bradford_curve.csv', header = F)
OD <- (as.numeric(table1[1,]))
curve <- data.frame (Concentrations = c(2, 4, 8, 10, 15, 20, 25, 30, 0), OD)

curve <- curve[-5,]

ggplot(curve, aes(curve$Concentrations, curve$OD)) +
  geom_smooth(method='lm') + 
  geom_text(aes(label = curve$Concentrations)) +
  theme_bw()

curve.lm <- lm(Concentrations ~ OD, data = curve)
summary(curve.lm)
coeffs <- coefficients(curve.lm)

#Теперь посмотрим что в образцах:

raw_data <- read.csv('/home/nickolay/Рабочий стол/samples.csv', header = F)
raw_data <- t(raw_data)
concentrations <- coeffs[1] + coeffs[2]*raw_data

for_WB <- data.frame(Names = c('H1299', 'H1299+scsh',
                                 'H1299+p3ex2', 'H1299+ACTN4',
                                 'H1299+shACTN4', 'H1299+p65'),
                     concentrations)

#for_WB$Concentrations <- (as.numeric(for_WB$Concentrations)/1000)
total_volume <- 25 #объем пробы полный
for_WB$Laemmly <- total_volume/5
for_WB$Sample <- min(for_WB$concentrations, na.rm = T)*(total_volume-for_WB$Laemmly)/for_WB$concentrations
for_WB$Water <- total_volume - for_WB$Sample - for_WB$Laemmly
for_WB

#Посчитаем относительно OD:
od_data <- read.csv('/home/nickolay/Рабочий стол/samples.csv', header = F)
od_data <- t(raw_data)

relative_WB <- data.frame(Names = c('H1299', 'H1299+scsh',
                                              'H1299+p3ex2', 'H1299+ACTN4',
                                              'H1299+shACTN4', 'H1299+p65'))
total_volume <- 5 #объем пробы полный
relative_WB$OD <- t(od_data)
#relative_WB$Laemmly <- total_volume/5
relative_WB$Laemmly <- 0
relative_WB$Sample <- min(relative_WB$OD)*(total_volume-relative_WB$Laemmly)/relative_WB$OD
relative_WB$Water <- total_volume - relative_WB$Sample - relative_WB$Laemmly
relative_WB

