#args = commandArgs(trailingOnly=TRUE)

#setwd('/projects/covid-ct/imlab/users/temi/projects/TFXcan/scripts/')

library(glue)
library(R.utils)
library(data.table)
library(glmnet)
#library(bit64)
library(foreach)
library(doParallel)

project_dir <- '/projects/covid-ct/imlab/users/temi/projects/TFXcan'
kawakami_data_dir <- glue('{project_dir}/train-test-val/kawakami-human')
TF <- 'FOXA1'

# using the aggByCenter 
kawakami_center_dt <- data.table::fread(glue('{kawakami_data_dir}/aggByCenter_{TF}_40000.csv.gz'))

# split the data
X_train <- kawakami_center_dt[, -c(1,2)] |> as.matrix()
y_train <- kawakami_center_dt[, c(1,2)] |> as.data.frame()

cat('Starting to build enet model\n')

# register a parallel backend
cl <- 512
doParallel::registerDoParallel(cl)

cat(glue('[REGISTERING] {foreach::getDoParWorkers()} workers/cores\n'))

set.seed(26102022)
enet_center_binary_model <- glmnet::cv.glmnet(x=X_train, y=y_train$class, family = "binomial", type.measure = "auc", alpha = 0.5, keep=T, parallel=T)

#registerDoSEQ()

#stopCluster(cl)

cat('[DONE] Now saving...\n')

# save the model
saveRDS(enet_center_binary_model, file=glue('{project_dir}/models/kawakami-human/enet-center-binary-40000.rds'))

cat('[JOB COMPLETE]')