# Author: Temi
# Date: Thursday July 27 2023
# Description: script to evaluate TFPred models on train and test data
# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--logistic_model", help='A TFPred model'),
    make_option("--train_data_file", help='training data file'),
    make_option("--test_data_file", help='test data file'),
    make_option("--eval_output", help='evaluation file in without .rds extension')
)

opt <- parse_args(OptionParser(option_list=option_list))

library(glue)
library(glmnet)
library(data.table)
library(parallel)
library(tidyverse)
library(pROC)

models <- list()
models[['logistic']] <- readRDS(opt$logistic_model)

for(i in seq_along(models)){
    # read in the newx data : train or test
    newx_list <- purrr::map(.x=c(opt$train_data_file, opt$test_data_file), function(each_data){
        mat_dt <- data.table::fread(each_data)
        newx <- as.matrix(mat_dt[, -c(1:3)])
        # you only need one : link
        link_pred <- predict(models[[i]], newx, s = "lambda.min", type = 'link') |> as.vector()
        response_pred <- predict(models[[i]], newx, s = "lambda.min", type = 'response') |> as.vector()
        df <- mat_dt[, c(1:3)] |> as.data.frame()
        df$enpact_score <- link_pred
        df$probability <- response_pred
        colnames(df) <- c('locus', 'binding_class', 'binding_count', 'enpact_score', 'probability')
        return(df)
    }, .progress=T)

    data.table::fwrite(newx_list[[1]], glue('{opt$eval_output}.{names(models)[i]}.train_eval.txt.gz'), sep='\t', compress='gzip', quote=F, row.names=F)
    data.table::fwrite(newx_list[[2]], glue('{opt$eval_output}.{names(models)[i]}.test_eval.txt.gz'), sep='\t', compress='gzip', quote=F, row.names=F)

    ev <- purrr::map(newx_list, function(dt){
        tryCatch({
            pp <- with(dt, pROC::roc(binding_class, enpact_score, quiet = TRUE))
            pp_auc <- pp$auc |> as.numeric() |> round(3)
            pp_var <- pROC::var(pp)|> as.numeric()
            pp_ci <- with(dt, pROC::ci.auc(binding_class, enpact_score, conf.level = 0.95, quiet = T)) |> as.numeric() 
            low <- pp_ci[1]
            upp <- pp_ci[3]
            
            return(cbind(pp_auc, pp_var, low, upp))
            }, error=function(e){
            return(cbind(NA, NA, NA, NA))
        })
    })

    ev <- do.call(rbind, ev) |> as.data.table() 
    rownames(ev) <- c('train', 'test')
    print(ev)
    ev %>%
        tibble::rownames_to_column('data') %>%
        data.table::fwrite(glue('{opt$eval_output}.{names(models)[i]}.eval.txt.gz'), sep='\t', compress='gzip', quote=F, row.names=F)
}