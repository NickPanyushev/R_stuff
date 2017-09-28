#!/usr/bin/R

#install.packages(c("dplyr", "data.table", "parallel"), repos = "http://cran.us.r-project.org",lib = "~/R_stuff/Rpackages/")

library(dplyr)
library(data.table)
library(parallel)

all_bins <- fread("all_bins.tsv")
ncells <- length(unique(all_bins$Cell_line))
#Скрипт отжирает очень много оперативы, не стоит его запускать, если оперативы меньше 20гб

#Сделаем лист с векторами, содержащими бины
#Имя каждого элемента - имя линии
bin_vecs <- list()
for (i in unique(all_bins$Cell_line)){
  bin_vecs[[i]] <- all_bins[Cell_line == i, Peak_Bins]
}

#Сделаем датафрейм с именами линий в сочетании в каждой ячейке
#Имя каждого столбца - количество линий в сочетании

comb_df <- matrix(data=NA,
                  nrow = choose(ncells, ncells %/% 2),     
                  ncol = ncells,
                  dimnames = list(NULL, seq(from = 1, to = ncells)))
comb_df <- as.data.table(comb_df)

#Сделали функцию, которая перебирает все сочетания клеточных линий и заменяет имена линий бинами
my_fun1 <- function(x){
  return (combn(unique(all_bins$Cell_line), x, FUN = function(x) bin_vecs[x] , simplify = F))
}
rm(bin_vecs)

#Заполним датафрейм сочетаниями
for (i in 1:ncells){
  length <-  choose(ncells, i)
  comb_df[[i]][1:length] <- my_fun1(i)
}
gc()

#Теперь попробуем что-нибудь посчитать
cl <- makeCluster(detectCores() - 1)

comb_df <- parApply(cl, comb_df, 2, function(x) length(Reduce(intersect, x)))
stopCluster(cl)
#Запишем в файл
fwrite(comb_df, file = "bins_intersect.csv")