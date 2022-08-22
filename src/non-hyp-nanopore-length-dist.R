library(data.table)
library(ggplot2)

faidx <- fread("output/nanopore/extract-non-hyp/non-hyp-fasta-faidx.out")


read_lengths <- faidx$V2
read_lengths_dt <- data.table(read_lengths)

min(read_lengths$length) # 124bp
max(read_lengths$length) # 70,440 bp
mean(read_lengths$length) # 7305 bp

ggplot(read_lengths_dt, aes(x=read_lengths))+
  geom_histogram(binwidth=1000, color="black", fill="grey")+
  theme_bw()