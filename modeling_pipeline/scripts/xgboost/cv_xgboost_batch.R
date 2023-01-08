arguments <- commandArgs(trailingOnly=TRUE) #|> as.numeric()

library(glue)
library(xgboost)
library(data.table)
library(foreach)
library(parallel)
library(doParallel)

#print(glue('[INFO] Arguments are {arguments}'))

project_dir <- '/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline'
# date
run_date <- Sys.Date()

id <- 'kawakami'
TF <- 'FOXA1'

xgb_params_file <- arguments[1]
data_file <- arguments[2] #glue('{project_dir}/data/enet_data')
unique_id <- arguments[3]
batch_num <- arguments[4]

# read in the training data
data_file <- "/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/data/enet_data/data_2022-12-16/train_enet_2022-12-16.csv.gz"
xgb_params_file <- '/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/metadata/cv_split/xgb_gridcv_parameters_1.txt'

# do some checkings
if(!file.exists(xgb_params_file)){
    stop(glue('[ERROR] Parameter file: {xgb_params_file} does not exist.'))
} else {
    params_data <- data.table::fread(xgb_params_file, header=F)[1:5, ]
    print(head(params_data, 3))
}

if(!file.exists(data_file)){
    stop('[ERROR] Training data file does not exist.')
} else {
    train_data <- data.table::fread(data_file)
    X <- train_data[, -c(1,2,3)][1:10000, ] |> as.matrix()
    print(X[1:5, ])
    y <- train_data$class[1:10000]
}

#print(glue('[INFO] Parameter file is {xgb_params_file}'))
print(glue('[INFO] Batch number is {batch_num}'))

output_dir <- glue('{project_dir}/models/xgboost_cv')
if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive=T)
}


 #|> unname()

batch_num <- 0
# register workers/cores
hcores <- parallel::detectCores()
print(glue('[INFO] Detected {hcores} cores.'))

# # how many cores to register
# njobs <- nrow(params_data)
# if(njobs >= (hcores - 4)){
#     ncores <- (hcores - 4) - 4
# } else {
#     ncores <- njobs
# }

ncores <- 10
cl <- parallel::makeCluster(ncores, outfile = glue('{project_dir}/log/batch_{batch_num}.out'))
doParallel::registerDoParallel(cl)
print(glue('[INFO] Registered {ncores} cores.'))

showConnections()

#print(glue('[INFO] {ls(globalenv())}'))

#objects_to_export <- c()
packages_to_export <- c('xgboost', 'data.table', 'base')
xgb_grid_result <- foreach::foreach(i=1:nrow(params_data), .combine='rbind', .packages=packages_to_export) %dopar% {
    print(glue::glue('[INFO] Currently on row {i}'))

    dtrain <- xgboost::xgb.DMatrix(data=X, label=y)
    pparams <- params_data[i, ] |> unlist() |> unname()
    xgb_params <- list(colsample_bytree=pparams[1], max_depth=pparams[2], eta=pparams[3], gamma=pparams[4], lambda=pparams[5], alpha=pparams[6])
    params <- c(xgb_params, list(objective = "binary:logistic", eval_metric='auc', verbosity=1))

    xgb_model <- xgboost::xgb.cv(data=dtrain, params=params, nrounds=100, nfold=3, metrics='auc', early_stopping_rounds = 20, print_every_n=20)

    scores <- as.data.frame(xgb_model$evaluation_log)
    idmax_test_auc <- which.max(scores$test_auc_mean)
    test_auc <- scores$test_auc_mean[idmax_test_auc]
    train_auc <- scores$train_auc_mean[idmax_test_auc]

    return(c(test_auc=test_auc, train_auc=train_auc, unlist(xgb_params)))
}


# dtrain <- xgboost::xgb.DMatrix(data=X, label=y)
# i <- 1
# pparams <- params_data[i, ] |> unlist() |> unname()
# xgb_params <- list(colsample_bytree=pparams[1], max_depth=pparams[2], eta=pparams[3], gamma=pparams[4], lambda=pparams[5], alpha=pparams[6])
# params <- c(xgb_params, list(objective = "binary:logistic", eval_metric='auc', verbosity=1))
# xgb_model <- xgboost::xgb.cv(data=dtrain, params=params, nrounds=200, nfold=3, metrics='auc', early_stopping_rounds = 20, print_every_n=49)

parallel::stopCluster(cl)
doParallel::stopImplicitCluster()
foreach::registerDoSEQ()
showConnections()

unique_id <- 'test'

print(glue('[INFO] Saving cv result for batch {batch_num}'))
write.table(xgb_grid_result, file=glue('{output_dir}/cv_results_{unique_id}_{batch_num}_{run_date}.txt'), append=T, row.names=F, quote=F, col.names=F)
print(glue('[INFO] Job complete for row {batch_num}'))

# dtrain <- xgboost::xgb.DMatrix(data=X, label=y)
# xgb_grid_res <- parallel::mclapply(seq(1:nrow(params_data)), function(i){
#     print(glue('[INFO] Currently on row {i}'))

#     pparams <- params_data[i, ] |> unlist() |> unname()
#     xgb_params <- list(colsample_bytree=pparams[1], max_depth=pparams[2], eta=pparams[3], gamma=pparams[4], lambda=pparams[5], alpha=pparams[6])
#     params <- c(xgb_params, list(objective = "binary:logistic", eval_metric='auc', verbosity=1))

#     xgb_model <- xgboost::xgb.cv(data=dtrain, params=params, nrounds=10, nfold=3, metrics='auc', early_stopping_rounds = 5, print_every_n=5)

#     # scores <- as.data.frame(xgb_model$evaluation_log)
#     # idmax_test_auc <- which.max(scores$test_auc_mean)
#     # test_auc <- scores$test_auc_mean[idmax_test_auc]
#     # train_auc <- scores$train_auc_mean[idmax_test_auc]

#     # return(c(test_auc=test_auc, train_auc=train_auc, unlist(xgb_params)))
# }, mc.cores=nrow(params_data))