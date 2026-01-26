# Author: Temi
# Date: Thursday August 10 2023
# Description: script to create predictors, ground truth and info files
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--data_directory", help="The folder with the epigenome data"),
    make_option("--individuals_list", help='The list of individuals for which you want to predict'),
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
library(abind)

print(opt)

# sbatch /beagle3/haky/users/temi/projects/TFXcan/src/calculateEnpactScores.sbatch /beagle3/haky/users/temi/data/1KG_AR_prostate /beagle3/haky/users/temi/projects/TFXcan-snakemake/metadata/EUR_individuals.1KG.txt /beagle3/haky/users/temi/projects/Enpact/data/enpact/weights/ENPACT_734_2025-04-24.compiled_weights.with_intercept.lambda.1se.txt.gz /beagle3/haky/users/temi/projects/TFXcan/data/enpact_scores/ENPACT_734.1KG.CWAS_ARBS.predictions.2025-10-16.rds.gz _aggByCollect_AR_Prostate.csv AR_Prostate

# # setwd("/beagle3/haky/users/temi/projects/TFXcan-snakemake/")
# opt <- list()
# opt$enpact_weights <- "/beagle3/haky/users/temi/projects/Enpact/data/enpact/weights/ENPACT_734_2025-04-24.compiled_weights.with_intercept.lambda.1se.txt.gz"
# opt$data_directory <- "/beagle3/haky/users/temi/data/1KG_AR_prostate"
# opt$individuals_list <- '/beagle3/haky/users/temi/projects/TFXcan-snakemake/metadata/EUR_individuals.1KG.txt'
# opt$files_pattern <- '_aggByCollect_AR_Prostate.csv'
# opt$for_models <- 'AR_Prostate'





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
    weights_intercepts <- weights[1, , drop = F]
    weights_features <- weights[-1, , drop = F]
} else {
     # weights 
    print(glue::glue("INFO - predicting only for {opt$for_models} models"))
    weights <- data.table::fread(opt$enpact_weights) %>% 
        dplyr::select(-feature) %>% as.matrix()
    weights_intercepts <- weights[1, opt$for_models, drop = F]
    weights_features <- weights[-1, opt$for_models, drop = F]
}


Y_hats <- purrr:::map(.x=individual_enpact_features, .f = function(each_file){
    dt <- data.table::fread(each_file)
    X <- as.matrix(dt[, -c(1)])

    stopifnot(dim(X)[2] == dim(weights_features)[1])

    # prediction
    y_hat_noIntercept <- X %*% weights_features[, , drop = FALSE]
    y_hat <- apply(y_hat_noIntercept, 1, function(each_row){
        weights_intercepts + each_row
    }) %>% t()

    if(!is.null(opt$for_models)){
        y_hat <- t(y_hat)
    }

    colnames(y_hat) <- colnames(weights_features)
    rownames(y_hat) <- dt$id
    return(y_hat)

}, .progress = TRUE)



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