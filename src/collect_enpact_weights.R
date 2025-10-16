# Author: Temi
# Date: Friday July 26 2024
# Description: script to create predictors, ground truth and info files
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--models_metadata", help='A transcription factor e.g. AR'),
    make_option("--weights_file_basename", help='file to write weights to')
)

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)
library(stringr)
library(foreach)
library(glmnet)
library(glue)

# models metadata
models_mtdt <- data.table::fread(opt$models_metadata) #data.table::fread('/beagle3/haky/users/temi/projects/TFXcan/files/enpact_models.metadata.tsv') #data.table::fread(opt$models_metadata)

if(!dir.exists(dirname(opt$weights_file_basename))){
    dir.create(dirname(opt$weights_file_basename))
}

# register cores since I want to parallelize
num_clusters <- 8 #- 5 # 12 - 5
doParallel::registerDoParallel(num_clusters)

cat('INFO - Registering', num_clusters, 'clusters for a parallel run\n')

# read in the models and extract the weights

models_weights <- foreach::foreach(i=seq_len(nrow(models_mtdt)), .combine='cbind', .inorder=T) %do% {
    mtdt <- models_mtdt[i, ] |> unlist()
    mpath <- mtdt['path'] #model <- mtdt['model']
    model <- readRDS(mpath)
    weights <- coef(model, s = 'lambda.1se') %>% as.matrix() %>% as.data.table() %>% setNames(mtdt['model']) #|> as.data.frame() |> setNames(model)
    return(weights)
} %>% dplyr::mutate(feature = c('intercept', paste0('f_', seq_len(nrow(.)-1)))) %>% dplyr::relocate(feature)

data.table::fwrite(models_weights, file=glue("{opt$weights_file_basename}.lambda.1se.txt.gz"), sep='\t', row.names=F, col.names=T, quote=F, compress = 'gzip')


models_weights <- foreach::foreach(i=seq_len(nrow(models_mtdt)), .combine='cbind', .inorder=T) %do% {
    mtdt <- models_mtdt[i, ] |> unlist()
    mpath <- mtdt['path']
    model <- readRDS(mpath)
    weights <- coef(model, s = 'lambda.min') %>% as.matrix() %>% as.data.table() %>% setNames(mtdt['model'])
    return(weights)
} %>% dplyr::mutate(feature = c('intercept', paste0('f_', seq_len(nrow(.)-1)))) %>% dplyr::relocate(feature)

data.table::fwrite(models_weights, file=glue("{opt$weights_file_basename}.lambda.min.txt.gz"), sep='\t', row.names=F, col.names=T, quote=F, compress = 'gzip')

doParallel::stopImplicitCluster()

# file.size('/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2024-07-26/models/BARX1_Colon/BARX1_Colon_2024-07-26.logistic.rds')
