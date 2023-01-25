arguments <- commandArgs(trailingOnly=TRUE)
#print(arguments)

library(glue)
library(glmnet)
library(data.table)
library(parallel)
library(tidyverse)

model_dir <- arguments[1]
model_id <- arguments[2]
model_type <- arguments[3] # e.g. linear or binary
predict_on <- arguments[4]
individuals_data_dir <- arguments[5]
individuals_ground_truth_file <- arguments[6]
TF <- arguments[7]
run_date <- arguments[8]
output_dir <- arguments[9]

print(glue('model dir is {model_dir}\nTF is {TF}\nmodel_id is {model_id}\noutput directory is {output_dir}\nrun date is {run_date}\nmodel type is {model_type}\n\n'))

# individuals
#individuals <- c('LuCaP_136', 'LuCaP_141', 'LuCaP_167', 'LuCaP_145')
#individuals <- data.table::fread('/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFpred_pipeline/metadata/individuals.txt', header=F)
individuals <- data.table::fread(predict_on, header=F)
individuals <- individuals$V1#[-1]#[1:5]

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

# read in the individual's training data & merge with the ground truth and get the pscores
ground_truth_file <- data.table::fread(individuals_ground_truth_file) |> as.data.frame()

predictions_list <- parallel::mclapply(individuals, function(each_ind){

    ind_gt <- ground_truth_file[, c('region', each_ind)]

    out <- lapply(agg_methods, function(each_method){

        # merge with the ground truth
        mat_file <- glue('{individuals_data_dir}/{each_ind}_{each_method}_{TF}.csv.gz')
        mat_dt <- data.table::fread(mat_file)

        gt <- ind_gt[ind_gt$region %in% mat_dt$id, ] # match the available regions
        gt_dedup <- gt[!duplicated(gt[['region']]),] # remove duplicates (there should be none anyway)
        new_dt <- base::merge(gt_dedup, mat_dt, by.x='region', by.y='id') # merge by region
        colnames(new_dt) <- c('region', 'class', paste('f_', 1:(ncol(new_dt)-2), sep=''))

        # compute prediction scores
        newx <- as.matrix(new_dt[, -c(1:2)])
        if(model_type == 'linear'){

            # you only need one : link
            link_pred <- predict(models_list[[each_method]], newx, s = "lambda.1se", type = 'link') |> as.vector()
            df <- new_dt[, c(1:2)] |> as.data.frame()
            df$prediction_link <- link_pred

        } else if (model_type == 'logistic'){
            # return two things

            link_pred <- predict(models_list[[each_method]], newx, s = "lambda.1se", type = 'link') |> as.vector()
            response_pred <- predict(models_list[[each_method]], newx, s = "lambda.1se", type = 'response') |> as.vector()

            df <- new_dt[, c(1:2)] |> as.data.frame()
            df$prediction_link <- link_pred
            df$prediction_response <- response_pred
        }

        # pdt <- cbind(new_dt[, c(1:2)], prediction_scores) |> as.data.frame()
        # colnames(pdt) <- c('region', 'class', each_method)
        return(df)

    })

    # merge by region, and class => they should all be the same but this provides a sanity check
    #out <- out %>% purrr::reduce(dplyr::full_join, by = c('region', 'class'))]
    names(out) <- agg_methods
    return(out)

}, mc.cores=length(individuals))

names(predictions_list) <- individuals

# save the object to be read later
print(glue('[INFO] Saving `freedman_{TF}_{model_type}_evaluation_{run_date}.rds` to {output_dir}'))
rds_file <- glue('{output_dir}/freedman_{TF}_{model_type}_evaluation_{run_date}.rds')
saveRDS(predictions_list, file=rds_file)