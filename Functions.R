library(data.table)
library(dplyr)

read.narrowpeak <- function(filename, filtering=TRUE){
  raw_peaks <-  read.table(filename, header = F, sep = "\t")
  names(raw_peaks) <- c("Chromosome", "Start", "Stop", "Peak_name",
                        "display_int", "dot",
                        "fold_change", "p-value", "q-value",
                        "summit_position_from_start")
  
  if (filtering == TRUE){
    #Выбросим неинформативные колонки
    raw_peaks$dot = raw_peaks$display_int <- NULL
    #удаляем хромосомы, которые странно называются
    raw_peaks <- raw_peaks[grep(pattern = "chr\\w\\d?", raw_peaks$Chromosome), ]
    #Митохондриалку тоже выбросим, там не должно быть ТФ
    raw_peaks <- raw_peaks[!raw_peaks$Chromosome == "chrM",]
    #Переставим имена хромосом в нормальном порядке
    raw_peaks$Chromosome <- factor(raw_peaks$Chromosome, 
                                   levels =c("chr1","chr2","chr3","chr4","chr5",
                                             "chr6","chr7","chr8","chr9","chr10",
                                             "chr11", "chr12", "chr13", "chr14", "chr15",
                                             "chr16", "chr17", "chr18", "chr19", "chr20",
                                             "chr21", "chr22", "chrX", "chrY"))
    droplevels(raw_peaks)
    raw_peaks <- as.data.table(raw_peaks)
  }
  return(raw_peaks)
}

bin.peaks <- function(dataframe, window_size = 50, append = TRUE)
{
  #Эта функция выдает номера бинов в диапазоне [START, STOP]
  fun <- function(Start, Stop)
  { 
    bin_start <- floor(Start/window_size)*window_size
    bin_stop <- ceiling(Stop/window_size)*window_size - window_size
    positions <- seq(bin_start, bin_stop, by = window_size)
    return(positions)
  }
  
  if (append == TRUE){
    dataframe$Peak_Bins <- mapply(fun, dataframe$Start, dataframe$Stop)
    return(dataframe)
  }else{
    bin_vector <- unlist(mapply(fun, dataframe$Start, dataframe$Stop))
    bin_vector <- sort(bin_vector)
    return(bin_vector)
  }
}

#Функция, которая пересекает все бины в сочетаниях от 1 до пересечь все и считает количество общих
#Берет на вход датафрейм с бинами, которые надо пересечь. Возвращает файл .tsv - с таблицей
library(compiler)
bins_intersect <- function(df){
  library(data.table, quietly = T)
  library(parallel)
  
  # Initiate cluster
  cl <- makeCluster(detectCores() - 1)
  
  #Поставим ограничение по памяти в 10гб
  if (.Platform$OS.type == "unix"){
    devtools::install_github("krlmlr/ulimit")
    ulimit::memory_limit(10000) 
  }
  if (.Platform$OS.type == "windows"){
    memory.limit(size = 10000) 
  }
  
  #Сделаем лист с векторами, содержащими бины
  #Имя каждого элемента - имя линии
  bin_vecs <- list()
  for (i in unique(df$Cell_line)){
    bin_vecs[[i]] <- df[Cell_line == i, Peak_Bins]
  }
  
  #Функция, которая перебирает все сочетания клеточных линий, заменяет имена линий бинами, и сразу их пересекает
  #На вход она берет количество линий и выдает вектор с количеством общих бинов
  common_bin_number <- function(x){
    temp <- combn(unique(df$Cell_line), x, 
                  FUN = function(x) bin_vecs[x],
                  simplify = F)
    temp <- sapply(temp, function(x) length(Reduce(intersect, x)))
    
    return(temp)
  }
  
  #Сделаем датафрейм с именами линий в сочетании в каждой ячейке
  #Имя каждого столбца - количество линий в сочетании
  ncells <- length(unique(df$Cell_line))
  comb_df <- matrix(data=NA,
                    nrow = choose(ncells, ncells %/% 2),     
                    ncol = ncells,
                    dimnames = list(NULL, seq(from = 1, to = ncells)))
  comb_df <- as.data.table(comb_df)

  #Заполним датафрейм сочетаниями, с помощью функции common_bin_number
  
  for (i in 1:ncells){
    comb_df[[i]][1:choose(ncells, i)] <- common_bin_number(as.numeric(i))
  }
  rm(i)
  stopCluster(cl)
  registerDoSEQ()
  #Теперь переформуем в 2-столбчатый датасет
  comb_df <- as.data.table(t(comb_df))
  comb_df$Line_number <- c(1:ncells)

  comb_df <- melt(comb_df, id="Line_number")
  comb_df$Bin_numbers <- comb_df$value
  comb_df$variable <- NULL
  comb_df$value <- NULL
  comb_df <- comb_df[!is.na(comb_df$Bin_numbers)]

  #Запишем в файл
  fwrite(comb_df, file = "combinations.tsv", sep = "\t")
  print("Bins combinations were successfully written to combinations.tsv")
  rm(bin_vecs, comb_df, cl)
}
bins_intersect <- cmpfun(bins_intersect) # Bytecode compilation  

rm(list = setdiff(ls(), lsf.str()))