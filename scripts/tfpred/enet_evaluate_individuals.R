arguments <- commandArgs(trailingOnly=TRUE)
#print(arguments)

library(glue)
library(glmnet)
library(data.table)
library(parallel)
library(tidyverse)

<<<<<<< HEAD
model_dir <- arguments[1] # /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/TFPred_pipeline/models/cistrome_AR_aggByMeanCenter_linear_2023-01-24.rds
predict_on <- arguments[2] # '/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/compare_predictors/metadata/cwas_individuals.txt'
individuals_data_dir <- arguments[3] #'/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/predictions_folder/cwas_AR_Prostate/predictions_2023-05-16/aggregated_predictions'
output_dir <- arguments[4] # /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/output
agg_method <- arguments[5] # aggByMeanCenter
metainfo <- arguments[6] # e.g. AR_Prostate_1KG'
id_num <- arguments[7] # e.g. 00, 01, 02 e.t.c


# model_dir <- '/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFPred_models/models/cistrome_AR_Prostate_2023-05-16'
# predict_on <- '/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/compare_predictors/metadata/cwas_individuals.txt'
# individuals_data_dir <- '/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/predictions_folder/cwas_AR_Prostate/predictions_2023-05-16/aggregated_predictions'
# output_dir <- '/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/output'
# agg_method <- 'aggByMeanCenter'
# metainfo <- 'AR_Prostate_1KG'
# id_num <- 00
=======
model_rds <- arguments[1] # /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/TFPred_pipeline/models/cistrome_AR_aggByMeanCenter_linear_2023-01-24.rds
predict_on <- arguments[2] # /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/info/cwas_individuals.txt
model_type <- arguments[3] # linear
individuals_data_dir <- arguments[4] #'/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/predictions_folder/cwas_AR_Prostate/predictions_2023-05-16/aggregated_predictions'
output_rds <- arguments[5] # /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/output

# model_rds='/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/TFPred_pipeline/models/cistrome_AR_aggByMeanCenter_linear_2023-01-24.rds'
# predict_on='/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/metadata/cwas_individuals.txt'
# model_type='linear'
# individuals_data_dir='/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/predictions_folder/cwas_AR_Prostate/predictions_2023-05-16/aggregated_predictions'
# output_dir='/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/output'

# model_rds <- '/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFPred_models/models/cistrome_AR_Prostate_2023-05-16/aggByMeanCenter_AR_Prostate.logistic.rds' 
# predict_on <- '/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/metadata/cwas_individuals.txt' 
# model_type <- 'logistic' 
# individuals_data_dir <- '/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/predictions_folder/cwas_imputed_AR_Prostate/predictions_2023-05-30/aggregated_predictions' 
#output_rds <- '/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/output/imputed_cwas_predictions/aggByMeanCenter_AR_Prostate_cwas.logistic.rds'
>>>>>>> 8d1e50b0eb3dd8bdfc1db9105ea1984253c9762d

#print(glue('model dir is {model_dir}\nTF is {TF}\nmodel_id is {model_id}\noutput directory is {output_dir}\nrun date is {run_date}\nmodel type is {model_type}\n\n'))

# individuals
individuals <- data.table::fread(predict_on, header=F)
individuals <- individuals$V1 |> unique()#[1:5]
print(individuals[1:5])
<<<<<<< HEAD
#agg_method <- 'aggByMeanCenter'

# read in the models
# load all the models coefficients
models_list <- purrr::map(.x=agg_method, function(each_method){
    linear_model <- readRDS(glue('{model_dir}/{each_method}_{metainfo}.linear.rds'))
    logistic_model <- readRDS(glue('{model_dir}/{each_method}_{metainfo}.logistic.rds'))
    model <- list(linear=linear_model, logistic=logistic_model)
    return(model)
}, .progress=T)

names(models_list) <- agg_method

model_types <- c('linear', 'logistic')
=======
agg_methods <- 'aggByMeanCenter'

# read in the models
# load all the models coefficients
models_list <- purrr::map(.x=agg_methods, function(each_method){
    model <- readRDS(model_rds)
    return(model)
}, .progress=T)

names(models_list) <- agg_methods
>>>>>>> 8d1e50b0eb3dd8bdfc1db9105ea1984253c9762d

# read in the individual's training data & merge with the ground truth and get the pscores
#ground_truth_file <- data.table::fread(individuals_ground_truth_file) |> as.data.frame()

# first gather valid
valid_names <- c()
for(name in individuals){
<<<<<<< HEAD
    ind_gt <- glue('{individuals_data_dir}/{name}_{agg_method}_{metainfo}.csv')
    if(file.exists(ind_gt)){
        print(name)
=======
    ind_gt <- glue('{individuals_data_dir}/{name}_aggByMeanCenter_AR_Prostate.csv')
    if(file.exists(ind_gt)){
>>>>>>> 8d1e50b0eb3dd8bdfc1db9105ea1984253c9762d
        valid_names <- append(valid_names, name)
    }
}
print(glue('INFO - Found {length(valid_names)} valid individuals'))

<<<<<<< HEAD
# here I am using just the linear model
for(model_type in model_types[1]){
    #model_type <- 'linear'

    predictions_list <- parallel::mclapply(valid_names, function(each_ind){
        #each_ind <- 'DFCI_1636'
        ind_gt <- glue('{individuals_data_dir}/{each_ind}_{agg_method}_{metainfo}.csv')
        if(file.exists(ind_gt)){
            mat_dt <- data.table::fread(ind_gt)
        }
        #out <- lapply(agg_method, function(each_method){
=======
predictions_list <- parallel::mclapply(valid_names, function(each_ind){

    #ind_gt <- ground_truth_file[, c('region', each_ind)]
   # print(each_ind)

    ind_gt <- glue('{individuals_data_dir}/{each_ind}_aggByMeanCenter_AR_Prostate.csv')
    if(file.exists(ind_gt)){
        mat_dt <- data.table::fread(ind_gt)
    }

    out <- lapply(agg_methods, function(each_method){

>>>>>>> 8d1e50b0eb3dd8bdfc1db9105ea1984253c9762d
        # compute prediction scores
        newx <- as.matrix(mat_dt[, -c(1)])
        regions <- mat_dt[, 1]
        if(model_type == 'linear'){
<<<<<<< HEAD
            # you only need one : link
            link_pred <- predict(models_list[[agg_method]][[model_type]], newx, s = "lambda.1se", type = 'link') |> as.vector()
=======

            # you only need one : link
            link_pred <- predict(models_list[[each_method]], newx, s = "lambda.1se", type = 'link') |> as.vector()
>>>>>>> 8d1e50b0eb3dd8bdfc1db9105ea1984253c9762d
            df <- cbind(regions, link_pred) |> as.data.frame() #new_dt[, c(1:2)] |> as.data.frame()
            colnames(df) <- c('regions', 'prediction_link')

        } else if (model_type == 'logistic'){
            # return two things
<<<<<<< HEAD
            link_pred <- predict(models_list[[agg_method]][[model_type]], newx, s = "lambda.1se", type = 'link') |> as.vector()
            response_pred <- predict(models_list[[agg_method]][[model_type]], newx, s = "lambda.1se", type = 'response') |> as.vector()
            df <- cbind(regions, link_pred, response_pred) |> as.data.frame()
            colnames(df) <- c('regions', 'prediction_link', 'prediction_response')
        }

        # merge by region, and class => they should all be the same but this provides a sanity check
        return(df)

    }, mc.cores=16)

    #predictions_list <- list(a=1, b=2)
    names(predictions_list) <- valid_names
    out_dt <- purrr::reduce(predictions_list, dplyr::left_join, by=c('regions'))
    colnames(out_dt) <- c('locus', valid_names)

    # save the object to be read later
    output_file <- glue('{output_dir}/{agg_method}_{metainfo}_1KG.{model_type}.{id_num}.txt')
    if(!dir.exists(dirname(output_file))){
        dir.create(dirname(output_file))
    } else {
        print(glue('INFO - {dirname(output_file)} exists'))
    }

    print(glue('INFO - Saving to {output_file}'))
    data.table::fwrite(out_dt, file=output_file, quote=FALSE, row.names=FALSE)
}

=======

            link_pred <- predict(models_list[[each_method]], newx, s = "lambda.1se", type = 'link') |> as.vector()
            response_pred <- predict(models_list[[each_method]], newx, s = "lambda.1se", type = 'response') |> as.vector()
            df <- cbind(regions, link_pred, response_pred) |> as.data.frame()
            colnames(df) <- c('regions', 'prediction_link', 'prediction_response')
        }
        return(df)
    })

    # merge by region, and class => they should all be the same but this provides a sanity check
    names(out) <- agg_methods
    return(out)

}, mc.cores=16)

names(predictions_list) <- valid_names

print(glue('INFO - Length of predictions list is {length(predictions_list)}'))

# save the object to be read later

if(!dir.exists(dirname(output_rds))){
    dir.create(dirname(output_rds))
} else {
    print(glue('INFO - {dirname(output_rds)} already exists'))
}

print(glue('INFO - Saving to {output_rds}'))
#rds_file <- glue('{output_dir}/cwas_AR_Prostate_linear.rds')
saveRDS(predictions_list, file=output_rds)
>>>>>>> 8d1e50b0eb3dd8bdfc1db9105ea1984253c9762d
