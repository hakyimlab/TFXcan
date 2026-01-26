suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--full_data", help='data to train with enet'),
    make_option("--seed", help='', default = 2025),
    make_option("--num_bound", type="integer", help='How many cv folds?'),
    make_option("--output_basename", help='', type = 'character')
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

library(tidyr)
library(dplyr)
library(glue)
library(data.table)
library(GenomicRanges)
library(tidyverse)

# opt <- list()
# opt$full_data <- "/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2025-04-24/predictor_files/AR_Prostate.info.txt.gz"
# opt$seed <- 2025
# opt$num_bound <- 1000
# opt$output_basename <- ""


# opt$exclude <- 

valid_chromosomes <- paste0('chr', 1:22)

full_data <- data.table::fread(opt$full_data) %>% 
    dplyr::filter(!chr %in% c('chr9', 'chr22')) %>% 
    dplyr::filter(chr %in% valid_chromosomes)

set.seed(opt$seed)
sampled_data <- full_data %>%
    dplyr::group_by(binding_class) %>%
    dplyr::slice_sample(n = opt$num_bound) %>%
    with(., GenomicRanges::GRanges(seqnames = chr, IRanges::IRanges(start = start, end = end), binding_class = binding_class)) %>%
    GenomicRanges::resize(fix = 'center', width = 512) %>%
    as.data.table() %>%
    tidyr::unite('locus', seqnames:end, sep = '_', remove = F) %>%
    dplyr::slice_sample(n = nrow(.))

# print(head(sampled_data))


sampled_data %>%
    data.table::fwrite(glue("{opt$output_basename}.ground_truth.tsv"), sep = '\t', col.names = T, row.names = F, quote = F)
sampled_data %>%
    dplyr::select(locus) %>%
    data.table::fwrite(glue("{opt$output_basename}.loci.tsv"), sep = '\t', col.names = F, row.names = F, quote = F)


# predict the epigenome with Enformer








