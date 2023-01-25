# Predict on test, train

arguments <- commandArgs(trailingOnly=TRUE)
#print(arguments)

library(glue)
library(glmnet)
library(data.table)
library(parallel)
library(tidyverse)

model_dir <- arguments[1] # ./models
model_id <- arguments[2] # e.g. cistrome, kawakami
model_type <- arguments[3] # e.g. linear or logistic
predict_on <- arguments[4] # e.g. train, test
data_dir <- arguments[5] # 
TF <- arguments[6]
run_date <- arguments[7]
output_dir <- arguments[8] # ./model_evaluation

print(glue('model dir is {model_dir}\nTF is {TF}\nmodel_id is {model_id}\noutput directory is {output_dir}\nrun date is {run_date}\nmodel type is {model_type}\npredict_on is {predict_on}\n\n'))

# aggregation models
agg_center <- 'aggByCenter'
agg_precenter <- 'aggByPreCenter'
agg_postcenter <- 'aggByPostCenter'
agg_meancenter <- 'aggByMeanCenter'
agg_mean <- 'aggByMean'
agg_upstream <- 'aggByUpstream'
agg_downstream <- 'aggByDownstream'
agg_upstream_downstream <- 'aggByUpstreamDownstream'

agg_methods <- c(agg_postcenter, agg_meancenter, agg_mean, agg_upstream, agg_center, agg_downstream, agg_upstream_downstream, agg_precenter)

# read in the models
# load all the models coefficients
models_list <- purrr::map(.x=agg_methods, function(each_method){
    agg_rds <- glue('{model_dir}/{model_id}_{TF}_{each_method}_{model_type}_{run_date}.rds')
    model <- readRDS(agg_rds)
    return(model)
}, .progress=T)

names(models_list) <- agg_methods

# read in the newx data : train or test
newx_list <- purrr::map(.x=agg_methods, function(each_method){
    mat_file <- glue('{data_dir}/{predict_on}_{each_method}.csv.gz')
    mat_dt <- data.table::fread(mat_file)
    return(mat_dt)
}, .progress=T)
names(newx_list) <- agg_methods

# predict
predictions_list <- parallel::mclapply(agg_methods, function(each_method){

    newx <- as.matrix(newx_list[[each_method]][, -c(1:4)])
    if(model_type == 'linear'){

        # you only need one : link
        link_pred <- predict(models_list[[each_method]], newx, s = "lambda.1se", type = 'link') |> as.vector()
        df <- newx_list[[each_method]][, c(1:4)] |> as.data.frame()
        df$prediction_link <- link_pred

    } else if (model_type == 'logistic'){
        # return two things

        link_pred <- predict(models_list[[each_method]], newx, s = "lambda.1se", type = 'link') |> as.vector()
        response_pred <- predict(models_list[[each_method]], newx, s = "lambda.1se", type = 'response') |> as.vector()

        df <- newx_list[[each_method]][, c(1:4)] |> as.data.frame()
        df$prediction_link <- link_pred
        df$prediction_response <- response_pred
    }

    return(df)

}, mc.cores=length(agg_methods))

names(predictions_list) <- agg_methods

# save the object to be read later
print(glue('[INFO] Saving `{model_id}_{TF}_{predict_on}_{model_type}_{run_date}.rds` to {output_dir}')) # cistrome_FOXA1_train_linear_2023-12-12.rds
rds_file <- glue('{output_dir}/{model_id}_{TF}_{predict_on}_{model_type}_{run_date}.rds')
saveRDS(predictions_list, file=rds_file)