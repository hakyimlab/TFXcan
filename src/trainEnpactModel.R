# Author: Temi
# Date: Thursday July 27 2023
# Description: script to train elastic net TFPred models
# Usage: Rscript train_enet.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--train_data_file", help='data to train with enet'),
    make_option("--rds_file", help='.rds file to be created as the model'),
    make_option("--nfolds", type="integer", default=5L, help='How many cv folds?')
)

opt <- parse_args(OptionParser(option_list=option_list))


library(glue)
library(R.utils)
library(data.table)
library(glmnet)
library(doParallel)
library(parallel)

seed <- 2023
if(file.exists(opt$train_data_file)){
    print(glue('INFO - Reading train data...'))
    dt_train <- data.table::fread(opt$train_data_file)
} else {
    stop(glue('ERROR - Training data cannot be found.'))
}

# remove missing values
cc <- complete.cases(dt_train)
dt_train <- dt_train[cc, ]

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

train_methods <- c('logistic')

parallel::mclapply(train_methods, function(each_method){

    cl <- 6 #parallel::makeCluster(5)
    doParallel::registerDoParallel(cl)

    print(glue('INFO - Starting to build {each_method} enet model'))

    if(each_method == 'linear'){

        cv_model <- tryCatch({
            glmnet::cv.glmnet(x=X_train, y=nbc, family = "gaussian", type.measure = "mse", alpha = 0.5, keep=T, parallel=T, nfolds=opt$nfolds)
        }, error = function(e){
            print(glue('ERROR - {e}'))
            return(NULL)
        })
        save_name <- paste0(opt$rds_file, '.linear.rds', sep='') #gsub('.rds', '.linear.rds', opt$rds_file)
    } else if (each_method == 'logistic'){

        cv_model <- tryCatch({
            glmnet::cv.glmnet(x=X_train, y=y_train$binding_class, family = "binomial", type.measure = "auc", alpha = 0.5, keep=T, parallel=T, nfolds=opt$nfolds, trace.it=F)
        }, error = function(e){
            print(glue('ERROR - {e}'))
            return(NULL)
        })
        save_name <- opt$rds_file #gsub('.rds', '.logistic.rds', opt$rds_file)
    }
    print(cv_model)
    print(glue('INFO - Saving `{save_name}`'))
    #rds_file <- glue('{model_file_basename}.{each_method}.rds')
    saveRDS(cv_model, file=save_name)
    doParallel::stopImplicitCluster()

}, mc.cores=1)

print(glue('INFO - Finished with model training and saving'))
doParallel::stopImplicitCluster()