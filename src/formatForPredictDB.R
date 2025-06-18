# Author: Temi
# Date: Thursday August 10 2023
# Description: script to create predictors, ground truth and info files
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--enpact_scores_file", help='input: A file containing the Enpact scores for the models'),
    make_option("--formatted_escores_file", help='output: The formatted Enpact scores file'),
    make_option("--formatted_annot_file", help='output: The formatted annotation file'),
    make_option("--include_model", help="input: What model should be extracted from the enpact scores file")
)

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)
library(tidyverse)
library(glue)

# opt <- list()
# opt$include_models <- NULL
# opt$blacklist <- '/project/haky/users/temi/projects/TFXcan-snakemake/metadata/HLA_blocks.txt'
# opt$filtered_GWAS_file <- '/project/haky/users/temi/projects/TFXcan-snakemake/data/PCRISK_GWAS_2024-06-07/collection/pc_risk.filteredGWAS.txt.gz'
# opt$enpact_scores_file <- '/project/haky/users/temi/projects/TFXcan-snakemake/data/PCRISK_GWAS_2024-06-07/enpactdb/pc_risk.enpact_scores.txt.gz'
# dt$NAME[startsWith(dt$NAME, 'HSF1_Mammary-Gland')] |> print()
# "/project2/haky/temi/projects/enpact-predict-snakemake/output/baca_2024-10-04/baca_2024-10-04.enpact_scores.array.rds.gz"


rdsfile <- readRDS(opt$enpact_scores_file)
escores <- rdsfile[, , opt$include_model] %>% as.data.frame() %>% tibble::rownames_to_column('NAME') %>% as.data.table()
annot <- escores %>% 
    dplyr::select(NAME) %>%
    tidyr::separate(col = NAME, into=c('chr', 'start', 'end'), sep='_', remove=F) %>%
    dplyr::mutate(gene_id = NAME, gene_name = NAME, gene_type = 'protein_coding', start = as.numeric(start), end = as.numeric(end), chr = gsub('chr', '', chr)) %>%
    dplyr::select(chr, start, end, gene_id, gene_name, gene_type)

escores %>%
    data.table::fwrite(file = opt$formatted_escores_file, quote=F, row.names=F, col.names=T, sep='\t')
annot %>%
    data.table::fwrite(file = opt$formatted_annot_file, quote=F, row.names=F, col.names=T, sep='\t')