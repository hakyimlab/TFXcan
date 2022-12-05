


library(glue)
library(R.utils)
library(data.table)
library(glmnet)
#library(bit64)
library(foreach)
library(doParallel)
library(parallel)
library(xgboost)

project_dir <- '/projects/covid-ct/imlab/users/temi/projects/TFXcan'
kawakami_data_dir <- glue('{project_dir}/train-test-val/kawakami-human')
TF <- 'FOXA1'

# using the aggByCenter 
kawakami_center_dt <- data.table::fread(glue('{kawakami_data_dir}/aggByCenter_{TF}_40000.csv'))

# split the data
X_train <- kawakami_center_dt[, -c(1,2)] |> as.matrix()
y_train <- kawakami_center_dt[, c(1,2)] |> as.data.frame()

cat('[SPLITTING]\n')

# for optimization, I need a test set
set.seed(27102022)

# first bind the data 
dt_all <- cbind(y_train[, 2], X_train) |> as.data.frame()

# split and do your magic
split_dt <- lapply(split(dt_all, f=dt_all$V1), function(each_dt){

  dt_indices <- sample(1:nrow(each_dt), 5000)
  a <- each_dt[dt_indices, ] 
  b <- each_dt[-dt_indices, ]

  return(list(train=a, test=b))

})

train_dt <- rbind(split_dt[[1]][[1]], split_dt[[2]][[1]])
test_dt <- rbind(split_dt[[1]][[2]], split_dt[[2]][[2]])

train_dt <- train_dt[sample(1:nrow(train_dt)), ] |> as.matrix()
test_dt <- test_dt[sample(1:nrow(test_dt)), ] |> as.matrix()

cat('[STARTING TO TRAIN] \n')

# here we use a grid search to find the optimal hyperparameters
xgb_grid_search <- expand.grid(colsample_bytree = 0.7, max_depth=20, eta=c(0.01, 0.1, 1), gamma=c(5, 10, 20), lambda=c(1, 3, 5), alpha=c(0, 3, 7))

cat(glue('[DETECTING CORES] {detectCores()} \n'))

fk_cl <- parallel::makeForkCluster(detectCores() - 200, outfile=glue("{project_dir}/log/xgboost-progress.log"))
registerDoParallel(fk_cl)

cat(glue('[REGISTERING & USING] {foreach::getDoParWorkers()} workers/cores \n'))

#set.seed(27102022)
output_cv_result <- foreach(i=1:nrow(xgb_grid_search), .combine='rbind', .packages=c('xgboost'), .inorder=F) %dopar% {

  # this is important - you need to initialize this within the loop !!!
  dtrain <- xgboost::xgb.DMatrix(data=train_dt[, -1], label=train_dt[, 1])

  params <- list(max_depth=xgb_grid_search[i, 'max_depth'], eta=xgb_grid_search[i, 'eta'], colsample_bytree=xgb_grid_search[i, 'colsample_bytree'], gamma=xgb_grid_search[i, 'gamma'], lambda=xgb_grid_search[i, 'lambda'], alpha=xgb_grid_search[i, 'alpha'], nthreads=, verbosity=0, objective = "binary:logistic", eval_metric='auc', booster='gbtree', tree_method='gpu_hist')

  cat('Running', i, '\n')

  # # uses 5-folds cross validation
  xgb_model <- xgboost::xgb.cv(data=dtrain, params=params, nrounds=50, nfold=5, metrics='auc', early_stopping_rounds = 20, print_every_n = 10, nthread=24, verbose=T)

  scores <- as.data.frame(xgb_model$evaluation_log)
  test_auc <- tail(scores$test_auc_mean, 1)
  train_auc <- tail(scores$train_auc_mean, 1)

  return(
    c(test_auc, train_auc, 
    xgb_grid_search[i, 'max_depth'], xgb_grid_search[i, 'eta'], 
    xgb_grid_search[i, 'colsample_bytree'], xgb_grid_search[i, 'gamma'], 
    xgb_grid_search[i, 'lambda'], xgb_grid_search[i, 'alpha'])
  )
}

stopCluster(fk_cl)

cat('[ENDING TRAINING]... \n')

cat('[DONE] Now saving... \n')

# save the model
output <- list(cv_result=output_cv_result, train_data=train_dt, test_data=test_dt)

saveRDS(output, file=glue('{project_dir}/models/kawakami-human/xgboost-model-cv.rds'))

cat('[JOB COMPLETE]')