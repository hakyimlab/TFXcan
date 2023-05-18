

args <- commandArgs(trailingOnly = TRUE)
print(args[1])

library(GenomicRanges)
library(data.table)
library(tidyverse)
library(glue)

bfile <- data.table::fread(args[1], header=FALSE)[, 1:3]
colnames(bfile) <- c('chr', 'start', 'end')
print(dim(bfile))
rfile_dir <- dirname(args[1])
rfile_bname <- basename(args[1]) %>% gsub('.sorted.bed', '.reduced.bed', .)
print(rfile_bname)

intersection_granges <- with(bfile, GRanges(chr, IRanges(start, end), strand='+', score=0))
intersection_granges <- GenomicRanges::reduce(intersection_granges)
as.data.frame(intersection_granges) %>% data.table::fwrite(., file=glue('{rfile_dir}/{rfile_bname}'), quote=F, sep='\t', row.names = FALSE)