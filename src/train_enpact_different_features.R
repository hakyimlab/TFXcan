# Author: Temi
# Date: Thursday July 27 2023
# Description: script to train elastic net TFPred models
# Usage: Rscript train_enet.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--train_data_file", help='data to train with enet'),
    make_option("--test_data_file", help='data to train with enet'),
    make_option("--features_to_remove", help='data to train with enet'),
    make_option("--nfolds", type="integer", default=5L, help='How many cv folds?'),
    make_option("--metadata", help='data to train with enet'),
    make_option("--weights_file", help='data to train with enet'),
    make_option("--predictions_file", help='data to train with enet'),
    make_option("--model_rds_file", help='data to train with enet')
)

opt <- parse_args(OptionParser(option_list=option_list))

library(tidyverse)
library(glue)
library(R.utils)
library(data.table)
library(glmnet)
library(doParallel)
library(parallel)

remove_features <- strsplit(opt$features_to_remove, ',')[[1]]

seed <- 2023
if(file.exists(opt$train_data_file)){
    dt_train <- data.table::fread(opt$train_data_file) %>% dplyr::select(-all_of(remove_features))
    print(glue('INFO - Read training data, and removed supplied features, new features = {ncol(dt_train) - 3}'))
} else {
    stop(glue('ERROR - Training data cannot be found.'))
}

if(file.exists(opt$test_data_file)){
    
    dt_test <- data.table::fread(opt$test_data_file) %>% dplyr::select(-all_of(remove_features))
    print(glue('INFO - Read test data, and removed supplied features, new features = {ncol(dt_test) - 3}'))
} else {
    stop(glue('ERROR - Test data cannot be found.'))
}

# remove duplicates
dt_train <- dplyr::distinct(dt_train) #%>% dplyr::slice_sample(n = 2000)
dt_test <- dplyr::distinct(dt_test)

# remove missing values
cc <- complete.cases(dt_train)
dt_train <- dt_train[cc, ]
cc <- complete.cases(dt_test)
dt_test <- dt_test[cc, ]

# split the data
X_train <- dt_train[, -c(1,2,3)] |> as.matrix()
y_train <- dt_train[, c(1,2,3)] |> as.data.frame()
vbc <- y_train$binding_counts
#nbc <- (vbc - min(vbc))/(max(vbc) - min(vbc)) # min-max normalization
nbc <- log10(1 + y_train$binding_counts)

cl <- 12 #parallel::makeCluster(5)
print(glue('INFO - Found {parallel::detectCores()} cores but using {cl}'))

set.seed(seed)

doParallel::registerDoParallel(cl)
print(glue('INFO - training enet model'))
cv_model <- glmnet::cv.glmnet(x=X_train, y=y_train$binding_class, family = "binomial", type.measure = "auc", alpha = 0.5, keep=T, parallel=T, nfolds=opt$nfolds, trace.it=F)

# save_name <- paste0(opt$rds_file, '.logistic.rds', sep='')

print(cv_model)
print(glue('INFO - Saving the model to `{opt$model_rds_file}`'))
#rds_file <- glue('{model_file_basename}.{each_method}.rds')
saveRDS(cv_model, file=opt$model_rds_file)
doParallel::stopImplicitCluster()

# cv_model <- readRDS('/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2025-04-24/models/AR_Prostate/AR_Prostate_2025-04-24.logistic.rds')
weights <- stats::coef(cv_model, s = 'lambda.1se')[-1]
weights <- weights %>% as.data.frame()
colnames(weights) <- opt$metadata  
weights <- weights %>% 
    dplyr::mutate(feature = paste0('f_', seq_len(nrow(.)))) %>% 
    dplyr::relocate(feature)

print(glue('INFO - Saving the weights'))
data.table::fwrite(weights, file=glue("{opt$weights_file}"), sep='\t', row.names=F, col.names=T, quote=F, compress = 'gzip')

weights <- weights %>% tibble::column_to_rownames('feature') %>% as.matrix()

print(glue('INFO - Evaluating on test data'))
X_test <- dt_test[, -c(2,3)] %>% tibble::column_to_rownames('locus') %>% as.matrix()
y_test <- dt_test[, c(1,2,3)] |> as.data.frame()
predictions <- as.matrix(X_test %*% weights)
predictions <- predictions[,,drop = FALSE] %>% as.data.table(keep.rownames = 'locus')
predictions <- dplyr::inner_join(predictions, y_test, by = 'locus')
print(glue('INFO - Saving the predicted values'))
data.table::fwrite(predictions, file=glue("{opt$predictions_file}"), sep='\t', row.names=F, col.names=T, quote=F, compress = 'gzip')