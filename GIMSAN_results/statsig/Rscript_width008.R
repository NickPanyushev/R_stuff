#===================================================================

# Reference filename: motif-00000.stdout
#
source("/home/npanyushev/gimsan_cmdline/bin/conf_pval_only.R")
library(MASS)

sample<-scan("/home/npanyushev/gimsan_cmdline/20170820_common_bins/statsig/scores.width008")
getConfPvalLat(21.41, sample, conf=0.1, mins=7, maxs=200)
