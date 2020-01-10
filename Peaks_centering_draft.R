Так как у нас есть саммиты для каждого из пиков, то сделаем .bed, где все участки общие будут отцентрованы по позиции саммита.
Заодно удлинним до 500 и 1000 нуклеотидов общие участки, во избежание всяких погрешностей 
```{r 500_centered}
region_length <- 500
for_fasta <- common_peaks[Cell_line == "A549", .(Chromosome, Start, Stop, Summit)]
for_fasta$Start <- for_fasta$Summit - region_length/2
for_fasta$Stop <- for_fasta$Summit + region_length/2
for_fasta$Summit <- NULL
stopifnot(for_fasta$Stop - for_fasta$Start  == region_length) #Проверка, все ли у нас получилось

# Сделаем .bed
fwrite(for_fasta, file = "BEDs/500_centered.bed", sep = "\t", col.names = FALSE)
rm(region_length, for_fasta)
```

```{r 1000_centered}
region_length <- 1000
for_fasta <- common_peaks[Cell_line == "A549", .(Chromosome, Start, Stop, Summit)]
for_fasta$Start <- for_fasta$Summit - region_length/2
for_fasta$Stop <- for_fasta$Summit + region_length/2
for_fasta$Summit <- NULL
stopifnot(for_fasta$Stop - for_fasta$Start  == region_length) #Проверка, все ли у нас получилось

