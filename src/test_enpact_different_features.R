# Author: Temi
# Date: Thursday July 27 2023
# Description: script to train elastic net TFPred models
# Usage: Rscript train_enet.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--test_data_file", help='data to train with enet'),
    make_option("--features_to_remove", help='data to train with enet', default = NULL),
    make_option("--metadata", help='data to train with enet'),
    make_option("--weights_file", help='data to train with enet'),
    make_option("--predictions_file", help='data to train with enet')
)

opt <- parse_args(OptionParser(option_list=option_list))

library(tidyverse)
library(glue)
library(R.utils)
library(data.table)
library(glmnet)
library(doParallel)
library(parallel)


# rr <- "f_37,f_38,f_127,f_142,f_191,f_654,f_1076,f_1213,f_1530,f_1837,f_1873,f_2061,f_2068,f_2091,f_2140,f_2223,f_2272,f_2305,f_2323,f_2380,f_2414,f_2489,f_2549,f_2645,f_2682,f_2813,f_2867,f_2872,f_2994,f_3201,f_3301,f_3352,f_3527,f_3679,f_3764,f_3807,f_3809,f_3855,f_3890,f_4005,f_4015,f_4062,f_4068,f_4093,f_4096,f_4166,f_4267,f_4271,f_4348,f_4355,f_4478,f_4691,f_4786,f_4787,f_4825,f_4839"
remove_features <- strsplit(opt$features_to_remove, ',')[[1]] #
# data.table::fread(opt$test_data_file) 
seed <- 2023
if(file.exists(opt$test_data_file)){
    dt_test <- data.table::fread(opt$test_data_file) %>% dplyr::select(-all_of(remove_features))
    print(glue('INFO - Read test data, and removed supplied features, new features = {ncol(dt_test) - 3}'))
} else {
    stop(glue('ERROR - Test data cannot be found.'))
}

# remove duplicates
#%>% dplyr::slice_sample(n = 2000)
dt_test <- dplyr::distinct(dt_test)

# remove missing values
dt_test <- dt_test[complete.cases(dt_test), ]

print(glue('INFO - Reading the weights'))
weights <- data.table::fread(glue("{opt$weights_file}")) #

weights.matrix <- weights %>% tibble::column_to_rownames('feature') %>% as.matrix()
weights_intercepts <- weights.matrix[1, , drop = F]
weights_features <- weights.matrix[-1, , drop = F]

print(glue('INFO - Evaluating on test data'))
X_test <- dt_test[, -c(2,3)] %>% tibble::column_to_rownames('locus') %>% as.matrix()
y_test <- dt_test[, c(1,2,3)] |> as.data.frame()

stopifnot(dim(X_test)[2] == dim(weights_features)[1])

predictions <- X_test %*% weights_features[, , drop = FALSE]
y_hat <- apply(predictions, 1, function(each_row){
    weights_intercepts + each_row
}) %>% as.matrix()

colnames(y_hat) <- colnames(weights_features)

y_hat <- y_hat[,,drop = FALSE] %>% as.data.table(keep.rownames = 'locus')
predictions <- dplyr::inner_join(y_hat, y_test, by = 'locus')
print(glue('INFO - Saving the predicted values'))
data.table::fwrite(predictions, file=glue("{opt$predictions_file}"), sep='\t', row.names=F, col.names=T, quote=F, compress = 'gzip')

# Y_hats <- purrr:::map(.x=individual_enpact_features, .f = function(each_file){
#     dt <- data.table::fread(each_file)
#     X <- as.matrix(dt[, -c(1)])

#     stopifnot(dim(X)[2] == dim(weights_features)[1])

#     # prediction
#     y_hat_noIntercept <- X %*% weights_features[, , drop = FALSE]
#     y_hat <- apply(y_hat_noIntercept, 1, function(each_row){
#         weights_intercepts + each_row
#     }) %>% t()

#     if(!is.null(opt$for_models)){
#         y_hat <- t(y_hat)
#     }

#     colnames(y_hat) <- colnames(weights_features)
#     rownames(y_hat) <- dt$id
#     return(y_hat)

# }, .progress = TRUE)