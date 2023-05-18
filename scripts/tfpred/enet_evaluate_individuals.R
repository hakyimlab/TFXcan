arguments <- commandArgs(trailingOnly=TRUE)
#print(arguments)

library(glue)
library(glmnet)
library(data.table)
library(parallel)
library(tidyverse)

model_rds <- arguments[1] # /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/TFPred_pipeline/models/cistrome_AR_aggByMeanCenter_linear_2023-01-24.rds
predict_on <- arguments[2] # /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/info/cwas_individuals.txt
model_type <- arguments[3] # linear
individuals_data_dir <- arguments[4] #'/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/predictions_folder/cwas_AR_Prostate/predictions_2023-05-16/aggregated_predictions'
output_dir <- arguments[5] # /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/output

# model_rds='/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/TFPred_pipeline/models/cistrome_AR_aggByMeanCenter_linear_2023-01-24.rds'
# predict_on='/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/metadata/cwas_individuals.txt'
# model_type='linear'
# individuals_data_dir='/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/predictions_folder/cwas_AR_Prostate/predictions_2023-05-16/aggregated_predictions'
# output_dir='/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/AR_prostate/output'

#print(glue('model dir is {model_dir}\nTF is {TF}\nmodel_id is {model_id}\noutput directory is {output_dir}\nrun date is {run_date}\nmodel type is {model_type}\n\n'))

# individuals
individuals <- data.table::fread(predict_on, header=F)
individuals <- individuals$V1[-1] |> unique()#[1:5]
print(individuals[1:5])
agg_methods <- 'aggByMeanCenter'

# read in the models
# load all the models coefficients
models_list <- purrr::map(.x=agg_methods, function(each_method){
    model <- readRDS(model_rds)
    return(model)
}, .progress=T)

names(models_list) <- agg_methods

# read in the individual's training data & merge with the ground truth and get the pscores
#ground_truth_file <- data.table::fread(individuals_ground_truth_file) |> as.data.frame()

predictions_list <- parallel::mclapply(individuals, function(each_ind){

    #ind_gt <- ground_truth_file[, c('region', each_ind)]
    print(each_ind)

    ind_gt <- glue('{individuals_data_dir}/{each_ind}_aggByMeanCenter_AR_Prostate.csv')
    if(file.exists(ind_gt)){
        mat_dt <- data.table::fread(ind_gt)
    }

    out <- lapply(agg_methods, function(each_method){

        # compute prediction scores
        newx <- as.matrix(mat_dt[, -c(1)])
        regions <- mat_dt[, 1]
        if(model_type == 'linear'){

            # you only need one : link
            link_pred <- predict(models_list[[each_method]], newx, s = "lambda.1se", type = 'link') |> as.vector()
            df <- cbind(regions, link_pred) |> as.data.frame() #new_dt[, c(1:2)] |> as.data.frame()
            colnames(df) <- c('regions', 'prediction_link')

        } else if (model_type == 'logistic'){
            # return two things

            link_pred <- predict(models_list[[each_method]], newx, s = "lambda.1se", type = 'link') |> as.vector()
            response_pred <- predict(models_list[[each_method]], newx, s = "lambda.1se", type = 'response') |> as.vector()

            df <- new_dt[, c(1:2)] |> as.data.frame()
            df <- cbind(regions, link_pred, response_pred) |> as.data.frame()
            colnames(df) <- c('regions', 'prediction_link', 'prediction_response')
        }

        # pdt <- cbind(new_dt[, c(1:2)], prediction_scores) |> as.data.frame()
        # colnames(pdt) <- c('region', 'class', each_method)
        return(df)

    })

    # merge by region, and class => they should all be the same but this provides a sanity check
    #out <- out %>% purrr::reduce(dplyr::full_join, by = c('region', 'class'))]
    names(out) <- agg_methods
    return(out)

}, mc.cores=16)

names(predictions_list) <- individuals

# save the object to be read later
print(glue('INFO - Saving `cwas_AR_Prostate_linear.rds` to {output_dir}'))
rds_file <- glue('{output_dir}/cwas_AR_Prostate_linear.rds')
saveRDS(predictions_list, file=rds_file)