# ARCHIVED — superseded by src/train_enpact.R. Kept for reference only.
# Author: Temi
# Date: Thursday July 27 2023
# Description: script to train elastic net TFPred models
# Usage: Rscript train_enet.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--train_data_file", help='Training data file with binding_class + feature columns'),
    make_option("--model_rds_file", help='.rds file to be created as the model'),
    make_option("--nfolds", type="integer", default=5L, help='How many cv folds?'),
    make_option("--weights_file", help='Output file to save the model weights (coefficients) to'),
    make_option("--metadata", help='Label to use for the weights column name (e.g. TF_tissue identifier)')
)

opt <- parse_args(OptionParser(option_list=option_list))

library(glue)
library(R.utils)
library(data.table)
library(glmnet)
library(doParallel)
library(parallel)
library(tidyverse)

seed <- 2023
if(file.exists(opt$train_data_file)){
    dt_train <- data.table::fread(opt$train_data_file)
    print(glue('INFO - Read training data, features = {ncol(dt_train) - 2}'))
} else {
    stop(glue('ERROR - Training data cannot be found.'))
}

dt_train <- dplyr::distinct(dt_train)
# remove missing values
cc <- complete.cases(dt_train)
dt_train <- dt_train[cc, ]

# split the data
X_train <- dt_train[, -c(1,2)] |> as.matrix()
y_train <- dt_train[, c(1,2)] |> as.data.frame()

cl <- 6 #parallel::makeCluster(5)
print(glue('INFO - Found {parallel::detectCores()} cores but using {cl}'))
doParallel::registerDoParallel(cl)

print(glue('INFO - training enet model'))
set.seed(seed)
cv_model <- glmnet::cv.glmnet(x=X_train, y=y_train$binding_class, family = "binomial", type.measure = "auc", alpha = 0.5, keep=T, parallel=T, nfolds=opt$nfolds, trace.it=F)
    
print(cv_model)
print(glue('INFO - Saving the model to `{opt$model_rds_file}`'))
saveRDS(cv_model, file=opt$model_rds_file)
doParallel::stopImplicitCluster()

weights <- stats::coef(cv_model, s = 'lambda.1se')
weights <- weights %>% as.data.frame()
colnames(weights) <- opt$metadata  
weights <- weights %>% 
    dplyr::mutate(feature = c('intercept', paste0('f_', seq_len(nrow(.))))) %>% 
    dplyr::relocate(feature)

print(glue('INFO - Saving the weights'))
data.table::fwrite(weights, file=glue("{opt$weights_file}"), sep='\t', row.names=F, col.names=T, quote=F, compress = 'gzip')