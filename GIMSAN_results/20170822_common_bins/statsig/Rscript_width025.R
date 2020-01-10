#===================================================================

# Reference filename: motif-00000.stdout
#
source("/home/nickolay/R_stuff/gimsan_cmdline/bin/conf_pval_only.R")
library(MASS)

sample<-scan("/home/nickolay/R_stuff/GIMSAN_results/20170822_common_bins/statsig/scores.width025")
getConfPvalLat(68.03, sample, conf=0.1, mins=7, maxs=200)
