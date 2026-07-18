# Author: Temi
# Date: Monday January 26 2026
# Description: inner-join ground truth locus labels with aggregated epigenome features, write merged training table
# Usage: Rscript merge_epigenome_ground_truth.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--epigenome", help='Path to aggregated epigenome feature matrix csv'),
    make_option("--ground_truth", help='Path to ground truth locus/binding_class tsv (from create_training_data.R)'),
    make_option("--output_filename", type="character", help='output file')
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

library(dplyr)
library(data.table)
library(tidyverse)
library(glue)

epigenome <- data.table::fread(opt$epigenome)
ground_truth <- data.table::fread(opt$ground_truth)
merged_data <- dplyr::inner_join(ground_truth, epigenome, by = c('locus' = 'id'))

print(merged_data[1:5, 1:5])
merged_data %>% data.table::fwrite(glue("{opt$output_filename}"), sep = '\t', col.names = T, row.names = F, quote = F)


 #(("/beagle3/haky/users/temi/projects/TFXcan/experiments/training_data/AR_Prostate_50.ground_truth.tsv"))
 #("/beagle3/haky/users/temi/projects/TFXcan/experiments/training_data/AR_Prostate_50.epigenome.csv.gz")#(opt$epigenome)