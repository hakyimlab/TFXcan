# Author: Temi
# Date: Monday January 26 2026
# Description: Builds the loci x TF/tissue z-ratio matrix from GWAS + TFXcan summary stats, then
#   splits the loci into random-subset batches for repeat_flash.R to run flashier on. First step of
#   the consensus matrix factorization ("tenerife") pipeline, driven by factorize_tfxcan.sbatch.
# Usage: Rscript prepare_summary_matrices.R --gwas_summary_statistic <file> --tfxcan_summary_statistic <file> --output_basename <name> --output_directory <dir> --temporary_directory <dir> --number_resamples <int> --number_batches <int>

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--gwas_summary_statistic", help='[Input] GWAS summary statistic prepare from TFXcan pipeline'),
    make_option("--tfxcan_summary_statistic", help='[Input] TFXcan summary statistics'),
    make_option("--output_basename", help='[Input] The output file basename'),
    make_option("--output_directory", help='[Input] The output directory'),
    make_option("--temporary_directory", help='[Input] Directory for intermediate files (random-subset batch rds/txt files consumed by repeat_flash.R)'),
    make_option("--number_resamples", default = 20, help='[Input] How many random loci subsets (repeats) to generate in total, spread across --number_batches'),
    make_option("--number_batches", default = 5, help='[Input] How many batches to split resamples into?')
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

library(glue)
library(data.table)
library(tidyverse)

# opt <- list()
# opt$gwas_summary_statistic <- "/beagle3/haky/users/temi/projects/TFXcan-snakemake/data/prostate_cancer_risk_2024-09-30/collection/prostate_cancer_risk.filteredGWAS.topSNPs.txt.gz"
# opt$tfxcan_summary_statistic <- "/beagle3/haky/users/temi/projects/Enpact/misc/reruns/prca_tfxcan/prostate_cancer_risk.enpactScores.2025-10-21.complete.spredixcan.tsv"
# opt$output_basename <- 'TEMO'
# opt$output_directory <- '/scratch/beagle3/temi/test'
# opt$temporary_directory <- '/scratch/beagle3/temi/temp'


dt_gwas <- data.table::fread(opt$gwas_summary_statistic) %>% dplyr::mutate(loci = paste('chr', chr, sep='') %>% paste(., pos, sep = ':'))

dt_tfxcan <- data.table::fread(opt$tfxcan_summary_statistic) %>%
    dplyr::filter(!is.na(pvalue), !is.na(zscore)) %>%
    tidyr::separate_wider_delim(tfbs, delim = "_", names = c("tf", 'tissue', 'chrom', 'start', 'end')) %>%
    dplyr::mutate(locus = paste(chrom, start, end, sep = "_"), tf_tissue = paste(tf, tissue, sep = "_")) %>%
    dplyr::select(locus, tf, tissue, tf_tissue, chrom, start, end, zscore, pvalue) %>%
    dplyr::mutate(snp = paste(chrom, start, sep = ':'))

# process the sumstats to make it easier to work with
res_sumstats <- dt_gwas %>%
    dplyr::rename(locus = loci, zscore = zscore, pvalue = pval) %>%
    dplyr::filter(nchar(a0) == 1) %>%
    dplyr::mutate(bind_level = 'GWAS') %>%
    dplyr::select(locus, bind_level, zscore, pvalue) %>% 
    dplyr::filter(locus %in% dt_tfxcan$snp)

if(nrow(res_sumstats) == 0){
    print(glue("ERROR - There are no shared SNPs between GWAS and TFXcan summary statistics."))
    quit(status = 1)
}


mtx <- dt_tfxcan %>% dplyr::select(locus = snp, bind_level = !!rlang::sym('tf_tissue'), zscore, pvalue) 

x_gwas <- res_sumstats %>% 
    dplyr::select(locus, bind_level, zscore, pvalue)
x_tfxcan <- dt_tfxcan %>%
    dplyr::select(locus = snp, bind_level = !!rlang::sym('tf_tissue'), zscore, pvalue, snp) 
mat_result <- dplyr::bind_rows(x_tfxcan, x_gwas) %>%
        dplyr::select(locus, bind_level, zscore) %>%
        tidyr::pivot_wider(id_cols = 'locus', values_from = zscore, names_from = bind_level) %>% 
        tibble::column_to_rownames('locus') %>%
        as.matrix() %>%
        t()

print(mat_result[1:3, 1:3]); print(dim(mat_result))

collect_matrices <- function(input_mat){
    collection_mat <- list()

    # split pdt into sTFXcan and GWAS
    mat_stfxcan <- input_mat[rownames(input_mat) != 'GWAS', ] |> t()
    # arrange by the GWAS pvalues
    mat_gwas <- input_mat['GWAS', ] |> sort(decreasing = T) |> as.matrix() 
    mat_stfxcan <- mat_stfxcan[rownames(mat_gwas), ]
    mat_gwas <- mat_gwas[rownames(mat_gwas), , drop = FALSE]

    collection_mat[["tfxcan"]] <- mat_stfxcan
    collection_mat[["gwas"]] <- mat_gwas
    colnames(collection_mat[["gwas"]]) <- "GWAS"

    return(collection_mat)

}

# collect the Z-scores into two matrices for TFXcan and GWAS results
data_list <- list()
data_list[['zratios']] <- collect_matrices(mat_result)

# divide or get eh Z-ratios
# z-ratio = (TFXcan association z-score)^2 / (GWAS z-score)^2 per locus x TF/tissue pair — this is
# the matrix that gets fed into flashier for factorization downstream.
zratio_matdev <- (data_list[["zratios"]][["tfxcan"]])^2/(data_list[["zratios"]][["gwas"]][,,drop = TRUE])^2

# there are missing values because I removed associations with pvalue at fdr < 0.05
dim(zratio_matdev); zratio_matdev[1:5, 1:5]


saveRDS(data_list, glue("{opt$output_directory}/{opt$output_basename}.zscores.matrices.rds"))
saveRDS(zratio_matdev, glue("{opt$output_directory}/{opt$output_basename}.zratios.matrices.rds"))

print(glue('INFO - Splitting zratios into {opt$number_resamples} random subsets'))

consensus_subsets_directory <- file.path(opt$temporary_directory)
dir.create(consensus_subsets_directory, recursive = TRUE)

# create random subsets; here I use 1000 random subsets, each with 80% of the loci
set.seed(2025)
random_loci <- purrr::map(1:opt$number_resamples, function(k){
    jj <- sample(nrow(zratio_matdev), round(nrow(zratio_matdev) * 0.8))
    return(jj)
})

# split the resamples into --number_batches groups (round-robin by index), so each batch can be
# run as a separate parallel repeat_flash.R job
batches_loci <- purrr::map(1:opt$number_batches, function(k){
    jj <- random_loci[seq(k, opt$number_resamples, by = opt$number_batches) ] #|> unlist() #seq(k, 1, 1) 
    return(jj)
})

names(batches_loci) <- paste0('batch', 1:opt$number_batches)

# save these as a list
out <- lapply(names(batches_loci), function(ni) {
    saveRDS(batches_loci[ni], glue::glue('{consensus_subsets_directory}/{opt$output_basename}.random_subsets.{ni}.rds'))
})

write(names(batches_loci), file = glue::glue('{consensus_subsets_directory}/{opt$output_basename}.random_subsets.txt'))
