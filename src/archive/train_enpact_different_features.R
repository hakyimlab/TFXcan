# ARCHIVED — superseded by src/train_enpact.R. Kept for reference only.
# Author: Temi
# Date: Monday January 26 2026
# Description: like trainEnpactModel.R, but trains the elastic net model after dropping a supplied list of
#   feature columns (feature-ablation experiment). Note: there's an early `quit()` call below that stops the
#   script right after printing opt — looks unfinished/left mid-debug.
# Usage: Rscript train_enet.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--train_data_file", help='Training data file with binding_class + feature columns'),
    make_option("--features_to_remove", help='Comma-separated list of feature column names to drop before training', default = NULL),
    make_option("--nfolds", type="integer", default=5L, help='How many cv folds?'),
    make_option("--metadata", help='Label to use for the weights column name (e.g. TF_tissue identifier)'),
    make_option("--weights_file", help='Output file to save the model weights (coefficients) to'),
    make_option("--model_rds_file", help='.rds file to be created as the model'),
    make_option("--ncores", default = 12)
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

quit()

library(tidyverse)
library(glue)
library(R.utils)
library(data.table)
library(glmnet)
library(doParallel)
library(parallel)

if(!is.null(opt$features_to_remove)){
    remove_features <- strsplit(opt$features_to_remove, ',')[[1]]
}

seed <- 2023
if(file.exists(opt$train_data_file)){
    dt_train <- data.table::fread(opt$train_data_file) %>% dplyr::select(-all_of(remove_features))
    if(!is.null(opt$features_to_remove)){
        print(glue('INFO - Read training data, and removed supplied features, new features = {ncol(dt_train) - 3}'))
    } else {
        print(glue('INFO - Read training data, number of features = {ncol(dt_train) - 3}'))
    }
} else {
    stop(glue('ERROR - Training data cannot be found.'))
}

# remove duplicates
dt_train <- dplyr::distinct(dt_train) #%>% dplyr::slice_sample(n = 2000)
# remove missing values
cc <- complete.cases(dt_train)
dt_train <- dt_train[cc, ]
# split the data
X_train <- dt_train[, -c(1,2,3)] |> as.matrix()
y_train <- dt_train[, c(1,2,3)] |> as.data.frame()

cl <- opt$ncores #parallel::makeCluster(5)
print(glue('INFO - Found {parallel::detectCores()} cores but using {cl}'))

set.seed(seed)

doParallel::registerDoParallel(cl)
print(glue('INFO - training enet model'))
cv_model <- glmnet::cv.glmnet(x=X_train, y=y_train$binding_class, family = "binomial", type.measure = "auc", alpha = 0.5, keep=T, parallel=T, nfolds=opt$nfolds, trace.it=F)

print(cv_model)
print(glue('INFO - Saving the model to `{opt$model_rds_file}`'))
saveRDS(cv_model, file=opt$model_rds_file)
doParallel::stopImplicitCluster()

#cv_model <- readRDS('/beagle3/haky/users/temi/projects/TFXcan/experiments/auprc/AR_Prostate.no_prostate_features.glmnet_model.rds')

weights <- stats::coef(cv_model, s = 'lambda.1se')
fnames <- rownames(weights)
fnames[1] <- 'intercept'
weights <- weights %>% as.matrix() %>% as.data.table()
colnames(weights) <- opt$metadata  
weights <- weights %>% 
    dplyr::mutate(feature = fnames) %>% 
    dplyr::relocate(feature)

print(glue('INFO - Saving the weights'))
data.table::fwrite(weights, file=glue("{opt$weights_file}"), sep='\t', row.names=F, col.names=T, quote=F, compress = 'gzip')

# weights <- weights %>% tibble::column_to_rownames('feature') %>% as.matrix()
# weights.matrix <- weights %>% dplyr::select(-feature) %>% as.matrix()
# weights_intercepts <- weights.matrix[1, , drop = F]
# weights_features <- weights.matrix[-1, , drop = F]

# print(glue('INFO - Evaluating on test data'))
# X_test <- dt_test[, -c(2,3)] %>% tibble::column_to_rownames('locus') %>% as.matrix()
# y_test <- dt_test[, c(1,2,3)] |> as.data.frame()

# stopifnot(dim(X_test)[2] == dim(weights_features)[1])

# predictions <- as.matrix(X_test %*% weights_features[, , drop = FALSE])
# y_hat <- apply(predictions, 1, function(each_row){
#     weights_intercepts + each_row
# }) %>% t()

# colnames(y_hat) <- colnames(weights_features)

# y_hat <- y_hat[,,drop = FALSE] %>% as.data.table(keep.rownames = 'locus')
# predictions <- dplyr::inner_join(y_hat, y_test, by = 'locus')
# print(glue('INFO - Saving the predicted values'))
# data.table::fwrite(predictions, file=glue("{opt$predictions_file}"), sep='\t', row.names=F, col.names=T, quote=F, compress = 'gzip')



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