
arguments <- commandArgs(trailingOnly=TRUE) |> as.numeric()

library(glue)
library(xgboost)
library(data.table)

project_dir <- '/grand/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline'
# date
run_date <- Sys.Date()

id <- 'kawakami'
TF <- 'FOXA1'

data_dir <- glue('{project_dir}/data/enet_data')
output_dir <- glue('{project_dir}/models/xgboost_cv')

# read in the training data
train_data <- data.table::fread(paste0(data_dir, '/train_enet.csv.gz'))

X <- train_data[, -c(1,2,3)] |> as.matrix()
y <- train_data$class
# set up param list
params <- list(colsample_bytree=arguments[1], max_depth=arguments[2], eta=arguments[3], gamma=arguments[4], lambda=arguments[5], alpha=arguments[6], verbosity=0, objective = "binary:logistic", eval_metric='auc', booster='gbtree')
dtrain <- xgboost::xgb.DMatrix(data=X, label=y)

# # uses 5-folds cross validation
cat('[INFO] Starting cross-validation\n')
xgb_model <- xgboost::xgb.cv(data=dtrain, params=params, nrounds=99, nfold=5, metrics='auc', early_stopping_rounds = 20, print_every_n = 10, nthread=2, verbose=T)
cat('[INFO] Cross-validation done\n')

scores <- as.data.frame(xgb_model$evaluation_log)
test_auc <- tail(scores$test_auc_mean, 1)
train_auc <- tail(scores$train_auc_mean, 1)

results <- c(test_auc, train_auc, arguments) |> as.data.frame() |> t()
cat('[INFO] Saving output\n')
write.table(results, file=glue('{output_dir}/cv_results_{run_date}.txt'), append=T, row.names=F, quote=F, col.names=F)
cat('[INFO] Job complete\n')