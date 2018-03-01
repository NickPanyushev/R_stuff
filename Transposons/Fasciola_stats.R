library(data.table)
setwd("~/R_stuff/Transposons/")
if (!file.exists("RM_rawout.zip")){
  unzip("f.hep_RM.out.zip")
  RM_rawout <- fread("fasciola_hepatica.PRJEB6687.WBPS9.genomic.fa.out" , sep = " ", blank.lines.skip = T, fill = T)
  RM_rawout <- RM_rawout[-1, -16]
  names(RM_rawout) <- c("SW_score", "div.%", "del.%", "ins.%", "scaffold_name", "query_start", "query_stop", 
                      "(left)", "strand", "element_name", "class/family", "repeat_start", "repeat_end", "(left)", "ID")
  RM_rawout$strand <- NULL
  RM_rawout$`(left)` <- NULL
  RM_rawout$`(left)` <- NULL
  RM_rawout$`class/family` <- sapply(RM_rawout$`class/family`, function(x) strsplit(x, "/")[[1]])
  RM_rawout$Class <- lapply(RM_rawout$`class/family`, function(x) x[[1]])
  RM_rawout$Family <- lapply(RM_rawout$`class/family`, function(x) ifelse(length(x) > 1, x[[2]], NA))
  RM_rawout$`class/family` <- NULL
  fwrite(RM_rawout, file = "RM_rawout.tsv", sep = "\t")
  zip("RM_rawout.zip", "RM_rawout.tsv")
  file.remove("RM_rawout.tsv")
  file.remove("fasciola_hepatica.PRJEB6687.WBPS9.genomic.fa.out")
}else{
  unzip("RM_rawout.zip")
  RM_rawout <- fread("RM_rawout.tsv")
  file.remove("RM_rawout.tsv")
}

#Отфильтруем все то, что нам не нужно
RM_rawout <- RM_rawout[Сlass != 'Simple_repeat']
RM_rawout <- RM_rawout[ 'Сlass' != 'Low_complexity']
RM_rawout <- RM_rawout[ 'Сlass' != 'tRNA']
RM_rawout <- RM_rawout[ 'Сlass' != 'rRNA']
RM_rawout <- RM_rawout[ 'Сlass' != 'snRNA']
RM_rawout <- RM_rawout[ 'Сlass' != 'srpRNA']
RM_rawout <- RM_rawout[ 'Сlass' != 'ARTEFACT']

