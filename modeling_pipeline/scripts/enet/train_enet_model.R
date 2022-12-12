#args = commandArgs(trailingOnly=TRUE)

#setwd('/projects/covid-ct/imlab/users/temi/projects/TFXcan/scripts/')


library(glue)
library(R.utils)
library(data.table)
library(glmnet)
library(foreach)
library(doParallel)

project_dir <- '/grand/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline'
setwd(project_dir)

TF <- 'FOXA1'
id <- 'kawakami'

# date
run_date <- Sys.Date()

data_dir <- glue('{project_dir}/data/enet_data')

dt_train <- data.table::fread(paste0(data_dir, '/train_enet.csv.gz'))
#dt_test <- data.table::fread(paste0(data_dir, '/test_enet.csv.gz'))

# split the data
X_train <- dt_train[, -c(1,2,3)] |> as.matrix()
y_train <- dt_train[, c(1,2,3)] |> as.data.frame()

# X_test <- dt_test[, -c(1,2,3)] |> as.matrix()
# y_test <- dt_test[, c(1,2,3)] |> as.data.frame()

cat('[INFO] Starting to build enet model\n')

# register a parallel backend
cl <- 24
doParallel::registerDoParallel(cl)

cat(glue('[INFO] Registering {foreach::getDoParWorkers()} workers/cores\n'))

set.seed(2022)
enet_center_binary_model <- glmnet::cv.glmnet(x=X_train, y=y_train$class, family = "binomial", type.measure = "auc", alpha = 0.5, keep=T, parallel=T, nfolds=10)

cat('[INFO] Saving the model\n')

# save the model
saveRDS(enet_center_binary_model, file=glue('{project_dir}/models/enet_models/{id}_center_binary_{run_date}.rds'))

registerDoSEQ()
cat('[INFO] Job completed')