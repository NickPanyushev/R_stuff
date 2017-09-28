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
