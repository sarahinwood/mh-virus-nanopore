library(data.table)
library(tidyverse)

hic_contigs <- fread("data/contig_ids/Mh_hic_contig_ids.txt", header=F)
hic_contigs$label <- paste("HiC")
viral_contigs <- fread("data/contig_ids/Mh_viral_contig_ids.txt", header=F)
viral_contigs$label <- paste("Viral")
contig_ids <- full_join(hic_contigs, viral_contigs)

coverage <- fread("output/long-reads/minimap2/samtools_coverage.out")
hic_viral_cov <- subset(coverage, `#rname` %in% contig_ids$V1)
hic_viral_cov_ids <- merge(hic_viral_cov, contig_ids, by.x="#rname", by.y="V1")

aggregate(list(hic_viral_cov_ids$coverage, hic_viral_cov_ids$meandepth), list(hic_viral_cov_ids$label), FUN=mean)
