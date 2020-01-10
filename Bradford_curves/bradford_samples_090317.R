library(dplyr)
library(ggplot2)
setwd('~/Рабочий стол/')

# Построим калибровочную кривулину

table1 <- read.csv('/home/nickolay/Рабочий стол/bradford_060317.17.txt', header = F) # считали из файла

#сделаем датафрейм для яичного белка
OD <- (as.numeric(table1[1, 2:11]))
curve_EGG <- data.frame (Concentrations = c( 0, 16, 8, 4, 2, 1, 0.5, 0.25, 0.125, 0.063), OD)

#сделаем датафрейм для bsa
OD <- (as.numeric(table1[2, 2:11]))
curve_BSA <- data.frame (Concentrations = c( 0, 16, 8, 4, 2, 1, 0.5, 0.25, 0.125, 0.063), OD)

# Построим график для яичного белка
EGG_iter1 <- ggplot(curve_EGG, aes(curve_EGG$Concentrations, curve_EGG$OD)) +
  geom_point()+
  geom_text(aes(label = curve_EGG$Concentrations)) +
  theme_bw()

# Выкинем последнюю точку 16мг, она явно выпадает
curve_EGG <- curve_EGG [-2,]

#Теперь посмотрим на статы  
curve_EGG.lm <- lm(Concentrations ~ OD, data = curve_EGG)
summary(curve_EGG.lm)
coeffs_EGG <- coefficients(curve_EGG.lm)

# Следующая версия графика:
EGG_iter2 <- ggplot(curve_EGG, aes(curve_EGG$Concentrations, curve_EGG$OD)) +
  geom_point()+
  geom_text(aes(label = curve_EGG$Concentrations)) +
  theme_bw()


# Теперь то же самое для BSA
# Построим график:
BSA_iter1 <- ggplot(curve_BSA, aes(curve_BSA$Concentrations, curve_BSA$OD)) +
              geom_point()+
              geom_text(aes(label = curve_BSA$Concentrations)) +
              theme_bw()

# Все очень плохо
# Выкинем последнюю точки 8 и 16мг, они явно выпадают
curve_BSA <- curve_BSA [c(-2,-3),]

#Теперь посмотрим на статы  
curve_BSA.lm <- lm(Concentrations ~ OD, data = curve_BSA)
summary(curve_BSA.lm)
coeffs_BSA <- coefficients(curve_BSA.lm)

# Построим график2:
BSA_iter2 <- ggplot(curve_BSA, aes(curve_BSA$Concentrations, curve_BSA$OD)) +
  geom_point()+
  geom_text(aes(label = curve_BSA$Concentrations)) +
  theme_bw()


#Теперь посмотрим что в образцах:

shACTN4 <- table1[,12]
concentrations_EGG <- coeffs_EGG[1] + coeffs_EGG[2]*shACTN4
concentrations_BSA <- coeffs_BSA[1] + coeffs_BSA[2]*shACTN4

for_WB <- data.frame(Names = c('sh3', 'sh1'),
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

