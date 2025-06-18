

# option_list <- list(
#     make_option("--data", help='data preferrably in .rds format of a matrix of GWAS loci by TF/tissue paris of ratios of z-scores'),
#     make_option("--splits", help='rds file of the splits'),
#     make_option("--batch", help='batch name'),
#     make_option("--output_basename", help='.rds file to be created as the model'),
#     make_option("--priorL", help='alpha value for enet', default="ebnm_point_exponential"),
#     make_option("--priorF", help='number of cores to use', default="ebnm_point_exponential"),
#     make_option("--transpose", type="logical", action='store_true', help='Should the data be transposed?'),
#     make_option("--greedy_Kmax", type="integer", default=100L, help='How many greedy iterations?'),
#     make_option("--subset", type="numeric", help='Subset of loci to use', default=0.8),
#     make_option("--niterations", type="numeric", default = 500L, help='Number of iterations')
# )

# opt <- parse_args(OptionParser(option_list=option_list))

# print(opt)

library(glue)
library(data.table)
library(flashier)
library(OmicKriging)
library(tidyverse)
library(ggplot2)

library(abind)


source("/beagle3/haky/users/temi/projects/Enpact/notebooks/hsq.R")

Ytrue <- data.table::fread("/beagle3/haky/users/temi/projects/Enpact/experiments/heritability/phenotypes/cwas_arbs.peaks.txt")
Ymat <- Ytrue %>% dplyr::select(-V1) %>% tibble::column_to_rownames('V2') %>% as.matrix()

# Ymat <- Ymat[, 1:5]

# est_mat <- vector("list", length = ncol(Ymat))
# se_mat <- vector("list", length = ncol(Ymat))
# pval_mat <- vector("list", length = ncol(Ymat))

# first read the chromosome grms 
valid_chromosomes <- c(1:22)
out <- purrr::map(valid_chromosomes, function(chr){
    grm_file <- glue("/project2/haky/Data/baca_cwas/vcfs/hg38/grm_files/chr{chr}.dose.grm")
    chrgrm <- OmicKriging::read_GRMBin(grm_file)

    # arrange grm to match the Ymat
    chrgrm <- chrgrm[match(rownames(Ymat), rownames(chrgrm)), match(rownames(Ymat), rownames(chrgrm)), drop = FALSE]

    chrhsq <- apply(Ymat, 2, function(Yarbs){
        calc_mle_from_grm(Yarbs, chrgrm, plotit = FALSE) |> as.matrix()
    })

    # colbind
    xout <- do.call('cbind', chrhsq)
    colnames(xout) <- colnames(Ymat)
    return(xout)
}, .progress = TRUE)

names(out) <- valid_chromosomes

# use abind to create a 3D array
out <- abind::abind(out, along = 0)

saveRDS(out, file = glue("/beagle3/haky/users/temi/projects/Enpact/experiments/heritability/per_arbs_per_chromosome_hsq_estimate.rds.gz"), compress = TRUE)




# Yone <- Yone[match(rownames(xgrm), rownames(Yone)), , drop = FALSE]
# Yone <- Yone[!is.na(Yone), , drop = FALSE]

# #
# xgrm_121 <- xgrm[match(rownames(Yone), rownames(xgrm)), match(rownames(Yone), rownames(xgrm)), drop = FALSE]