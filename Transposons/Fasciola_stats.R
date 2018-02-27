library(data.table)
RM_rawout <- fread("fasciola_hepatica.PRJEB6687.WBPS9.genomic.fa.out" , sep = " ", blank.lines.skip = T, fill = T)
RM_rawout <- RM_rawout[-1, -16]
names(RM_rawout) <- c("SW_score", "div.%", "del.%", "ins.%", "scaffold_name", "query_start", "query_stop", 
                      "(left)", "strand", "element_name", "class/family", "repeat_start", "repeat_end", "(left)", "ID")
RM_rawout$strand <- NULL
RM_rawout$`(left)` <- NULL
RM_rawout$`(left)` <- NULL
RM_rawout <- RM_rawout[`class/family`!= "Simple_repeat", ]
RM_rawout <- RM_rawout[`class/family`!= "Low_complexity", ]
RM_rawout$Class <- sapply(RM_rawout$`class/family` function(x) strsplit(x, "/")[1])
RM_rawout$Family <- strsplit(RM_rawout$`class/family`, "/")[2]
unique(RM_rawout$`class/family`)
