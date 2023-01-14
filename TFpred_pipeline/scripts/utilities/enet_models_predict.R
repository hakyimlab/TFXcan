arguments <- commandArgs(trailingOnly=TRUE)
#print(arguments)

library(glue)
library(glmnet)
library(data.table)
library(parallel)
library(tidyverse)

model_dir <- arguments[1]
model_id <- arguments[2]
individual_data_dir <- arguments[3]
TF <- arguments[4]
model_type <- arguments[5] # e.g. linear or binary
run_date <- arguments[6]
output_dir <- arguments[7]
ground_truth_file <- arguments[8]

print(glue('model dir is {model_dir}\nTF is {TF}\nmodel_id is {model_id}\noutput directory is {output_dir}\nrun date is {run_date}\nmodel type is {model_type}\n\n'))

# individuals
#individuals <- c('LuCaP_136', 'LuCaP_141', 'LuCaP_167', 'LuCaP_145')
individuals <- data.table::fread('/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFpred_pipeline/metadata/individuals.txt', header=F)
individuals <- individuals$V1[-1]#[1:5]

# aggregation models
agg_center <- 'aggByCenter'
agg_precenter <- 'aggByPreCenter'
agg_postcenter <- 'aggByPostCenter'
agg_mean <- 'aggByMean'
agg_upstream <- 'aggByUpstream'
agg_downstream <- 'aggByDownstream'
agg_upstream_downstream <- 'aggByUpstreamDownstream'

agg_methods <- c(agg_postcenter, agg_mean, agg_upstream, agg_center, agg_downstream, agg_upstream_downstream, agg_precenter)

# read in the models
# load all the models coefficients
linear_models_betas <- purrr::map(.x=agg_methods, function(each_method){
    agg_rds <- glue('{model_dir}/{model_id}_{TF}_{each_method}_{model_type}_{run_date}.rds')
    model <- readRDS(agg_rds)
    whlm <- which(model$lambda == model[['lambda.1se']])
    model_beta <- model$glmnet.fit$beta |> as.matrix()
    bst_beta <- model_beta[, whlm]
    return(bst_beta)
}, .progress=T)

names(linear_models_betas) <- agg_methods

# read in the individual's training data & merge with the ground truth and get the pscores
freedman_ground_truth <- data.table::fread(ground_truth_file) |> as.data.frame()

pscores_list <- parallel::mclapply(individuals, function(each_ind){

    ind_gt <- freedman_ground_truth[, c('region', each_ind)]

    out <- lapply(agg_methods, function(each_method){
        mat_file <- glue('{individual_data_dir}/{each_ind}_{each_method}_{TF}.csv.gz')
        mat_dt <- data.table::fread(mat_file)

        gt <- ind_gt[ind_gt$region %in% mat_dt$id, ] # match the available regions
        gt_dedup <- gt[!duplicated(gt[['region']]),] # remove duplicates (there should be none anyway)
        new_dt <- base::merge(gt_dedup, mat_dt, by.x='region', by.y='id') # merge by region
        colnames(new_dt) <- c('region', 'class', paste('f_', 1:(ncol(new_dt)-2), sep=''))

        # compute prediction scores
        model_betas <- linear_models_betas[[each_method]] |> as.matrix()
        new_x <- new_dt[, -c(1:2)] |> as.matrix()
        prediction_scores <- new_x %*% model_betas # e.g. 40000 by 5313 %*% 5313 by 1 => 40000 by 1

        pdt <- cbind(new_dt[, c(1:2)], prediction_scores) |> as.data.frame()
        colnames(pdt) <- c('region', 'class', each_method)
        return(pdt)

    })

    # merge by region, and class => they should all be the same but this provides a sanity check
    out <- out %>% purrr::reduce(dplyr::full_join, by = c('region', 'class'))
    return(out)

}, mc.cores=length(individuals))

names(pscores_list) <- individuals

# save the object to be read later
print(glue('[INFO] Saving `freedman_{TF}_prediction_scores_{run_date}.rds` to {output_dir}'))
rds_file <- glue('{output_dir}/freedman_{TF}_prediction_scores_{run_date}.rds')
saveRDS(pscores_list, file=rds_file)