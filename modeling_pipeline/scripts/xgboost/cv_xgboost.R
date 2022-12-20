arguments <- commandArgs(trailingOnly=TRUE) #|> as.numeric()

library(glue)
library(xgboost)
library(data.table)

xgb_arguments <- arguments[1:6]

#print(glue('[INFO] Arguments are {arguments}'))

project_dir <- '/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline'
# date
run_date <- Sys.Date()

id <- 'kawakami'
TF <- 'FOXA1'

data_file <- arguments[7] #glue('{project_dir}/data/enet_data')
unique_id <- arguments[8]
output_dir <- glue('{project_dir}/models/xgboost_cv')
if(!dir.exists(output_dir)){
    dir.create(output_dir, recursive=T)
}

# read in the training data
train_data <- data.table::fread(data_file)

X <- train_data[, -c(1,2,3)] |> as.matrix()
y <- train_data$class

# set up param list
params <- list(colsample_bytree=xgb_arguments[1], max_depth=xgb_arguments[2], eta=xgb_arguments[3], gamma=xgb_arguments[4], lambda=xgb_arguments[5], alpha=xgb_arguments[6], verbosity=0, objective = "binary:logistic", eval_metric='auc', booster='gbtree')
dtrain <- xgboost::xgb.DMatrix(data=X, label=y)

# # uses 5-folds cross validation
print(glue('[INFO] Starting cross-validation for row {arguments[9]}'))
xgb_model <- xgboost::xgb.cv(data=dtrain, params=params, nrounds=100, nfold=5, metrics='auc', early_stopping_rounds = 20, verbose=0)
#print('[INFO] Cross-validation done')

scores <- as.data.frame(xgb_model$evaluation_log)
idmax_test_auc <- which.max(scores$test_auc_mean)
test_auc <- scores$test_auc_mean[idmax_test_auc]
train_auc <- scores$train_auc_mean[idmax_test_auc]

results <- c(test_auc, train_auc, xgb_arguments) |> as.data.frame() |> t()
print(glue('[INFO] Saving cv result for row {arguments[9]}'))
write.table(results, file=glue('{output_dir}/cv_results_{unique_id}_{argument[10]}_{run_date}.txt'), append=T, row.names=F, quote=F, col.names=F)
print(glue('[INFO] Job complete for row {arguments[9]}'))