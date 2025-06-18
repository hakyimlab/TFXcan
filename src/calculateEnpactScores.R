# Author: Temi
# Date: Thursday August 10 2023
# Description: script to create predictors, ground truth and info files
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--data_directory", help='A transcription factor e.g. AR'),
    make_option("--individuals_list", help='An output file'),
    make_option("--enpact_weights"),
    make_option("--output_file"),
    make_option("--files_pattern"),
    make_option("--for_models", default = NULL),
    make_option("--epifeatures_file", default = NULL)
)

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)
library(tidyverse)
library(glmnet)
library(glue)

print(opt)

# # setwd("/beagle3/haky/users/temi/projects/TFXcan-snakemake/")
# opt <- list()
# opt$enpact_weights <- "/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2025-04-24/statistics/ENPACT_734_2025-04-24.compiled_weights.lambda.1se.txt.gz"
# opt$data_directory <- "/beagle3/haky/users/temi/data/baca_AR_prostate"
# opt$individuals_list <- '/project2/haky/temi/projects/enpact-predict-snakemake/metadata/cwas_individuals.txt'
# opt$output_file <- '/beagle3/haky/users/temi/projects/Enpact/misc/reruns/enpact_predictions/ENPACT_48.BACA.predictions.2025-04-28.rds.gz'
# opt$files_pattern <- '_aggByMeanCenter_AR_Prostate.csv'

# opt$data_directory <- "/beagle3/haky/users/temi/projects/TFXcan-snakemake/data/prostate_cancer_risk_2024-09-30/aggregated_predictions/prostate_cancer_risk"
# opt$individuals <- '/beagle3/haky/users/temi/projects/TFXcan-snakemake/metadata/EUR_individuals.1KG.txt'
# opt$output_file <- '/beagle3/haky/users/temi/projects/Enpact/misc/reruns/enpact_predictions/ENPACT_48.predictions.2025-04-24.rds.gz'

dt_individuals <- data.table::fread(opt$individuals_list, header = FALSE)
individual_enpact_features <- glue::glue("{opt$data_directory}/{dt_individuals$V1}{opt$files_pattern}") |> as.vector()
individual_enpact_features <- sapply(individual_enpact_features, Sys.glob) |> unname()
names(individual_enpact_features) <- dt_individuals$V1

print(length(individual_enpact_features))

individual_enpact_features <- lapply(individual_enpact_features, function(each_file){

    tryCatch({
        if(file.exists(each_file)){
            # print(glue::glue("INFO - {each_file}"))
            return(each_file)
        } 
    }, error = function(e){
        #message(glue::glue("ERROR - {e}"))
        return(NULL)
    })
})
individual_enpact_features <- Filter(Negate(is.null), individual_enpact_features)
individual_enpact_features <- unlist(individual_enpact_features)

print(length(individual_enpact_features))

if(is.null(opt$for_models)){
    # weights 
    print(glue::glue("INFO - predicting for all models supplied in the weights file"))
    weights <- data.table::fread(opt$enpact_weights) %>% 
        dplyr::select(-feature) %>% as.matrix()
} else {
     # weights 
    print(glue::glue("INFO - predicting only for {opt$for_models} models"))
    weights <- data.table::fread(opt$enpact_weights) %>% 
        dplyr::select(-feature) %>% as.matrix()
    weights <- weights[, opt$for_models, drop = FALSE]
}


Y_hats <- purrr:::map(.x=individual_enpact_features, .f = function(each_file){
    dt <- data.table::fread(each_file)
    X <- as.matrix(dt[, -c(1)])

    stopifnot(dim(X)[2] == dim(weights)[1])

    # prediction
    y_hat <- X %*% weights[, , drop = FALSE]
    rownames(y_hat) <- dt$id

    return(y_hat)

}, .progress = TRUE)

library(abind)

# collect the names of the locus and individuals
loci <- purrr::map(.x=Y_hats, .f=rownames) %>%
    base::Reduce(intersect, .) |> unique()
models <- purrr::map(.x=Y_hats, .f=colnames) %>%
    base::Reduce(intersect, .) |> unique()
individuals <- names(individual_enpact_features)

print(length(models))
print(length(individuals))

# ensure the dimensions are matched
out <- purrr::map(.x=Y_hats, function(each_dt){
    X <- as.matrix(each_dt)
    return(X[loci, models, drop = FALSE])
})

names(out) <- individuals

# combine into an array
myarray <- tryCatch({
    abind::abind(out, along=3)
}, error = function(e){
    abind::abind(out, along=2)
})

dimarray <- dim(myarray)
if(length(dimarray) == 2){
    #myarray <- array(myarray, dim=c(dimarray[1], 1, dimarray[2]))
    saveRDS(myarray, file = opt$output_file, compress = "gzip")
} else if(length(dimarray) == 3){
    reshapedarray <- aperm(myarray, c(1, 3, 2), resize=TRUE)
    # print(dimarray)
    # print(dimnames(reshapedarray))
    saveRDS(reshapedarray, file = opt$output_file, compress = "gzip")
}

print("INFO - ENPACT scores calculated and saved to file.")