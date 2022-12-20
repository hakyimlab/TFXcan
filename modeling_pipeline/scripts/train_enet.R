library(glue)
library(R.utils)
library(data.table)
library(glmnet)
library(doParallel)
library(parallel)

project_dir <- '/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline'
setwd(project_dir)

TF <- 'FOXA1'
id <- 'kawakami'

# date
run_date <- Sys.Date()

#data_dir <- glue('{project_dir}/data/enet_data')
data_file <- "/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/data/enet_data/data_2022-12-16/train_enet_2022-12-16.csv.gz" #arguments[1]
print(data_file)

dt_train <- data.table::fread(data_file)
print(dim(dt_train))

# split the data
X_train <- dt_train[, -c(1,2,3)] |> as.matrix()
y_train <- dt_train[, c(1,2,3)] |> as.data.frame()
print(head(y_train))

## metadata
tracks_info <- data.table::fread(glue('{project_dir}/metadata/enformer_tracks_annotated-resaved.txt'))
tracks_info$feature_names <- paste0('f_', 1:nrow(tracks_info))

categories <- tracks_info$assay |> unique()
categories <- categories[!is.na(categories)]


print(glue('[INFO] Found {parallel::detectCores()} cores\n'))

# register a parallel backend
# cl <- 24
# doParallel::registerDoParallel(cl)

#print(glue('[INFO] Registering {foreach::getDoParWorkers()} workers/cores\n'))

# build two models
print('[INFO] Starting to build enet model\n')
set.seed(2022)


grouped_models <- parallel::mclapply(categories, function(each_category){
    X_fnames <- X_train[, tracks_info[tracks_info$assay == each_category, ]$feature_names]

    mixing_parameters <- c(0, 0.25, 0.5, 0.75, 1)
    enet_center_binary_models_list <- parallel::mclapply(mixing_parameters, function(mix){

        cl <- 5 #parallel::makeCluster(5)
        doParallel::registerDoParallel(cl)
        #print(glue('[INFO] Registering {len(foreach::getDoParWorkers())} workers/cores for {mix} mixing parameter\n'))

        model <- glmnet::cv.glmnet(x=X_fnames, y=y_train, family = "binomial", type.measure = "auc", alpha = mix, keep=T, parallel=T, nfolds=3)
        print(model)
        doParallel::stopImplicitCluster()
        #registerDoSEQ()
        return(model)

    }, mc.cores=30)

    names(enet_center_binary_models_list) <- mixing_parameters

    return(enet_center_binary_models_list)
})
names(grouped_models) <- categories


# mixing_parameters <- c(0, 0.25, 0.5, 0.75, 1)
# enet_center_binary_models_list <- parallel::mclapply(mixing_parameters, function(mix){

#     cl <- 5 #parallel::makeCluster(5)
#     doParallel::registerDoParallel(cl)
#     #print(glue('[INFO] Registering {len(foreach::getDoParWorkers())} workers/cores for {mix} mixing parameter\n'))

#     model <- glmnet::cv.glmnet(x=X_train, y=y_train$class, family = "binomial", type.measure = "auc", alpha = mix, keep=T, parallel=T, nfolds=3)
#     print(model)
#     doParallel::stopImplicitCluster()
#     #registerDoSEQ()
#     return(model)

# }, mc.cores=48)

# names(enet_center_binary_models_list) <- mixing_parameters

print('[INFO] Finished with binary, mixing models\n')


# vbc <- y_train$binding_counts
# nbc <- (vbc - min(vbc))/(max(vbc) - min(vbc))
# enet_center_linear_model <- glmnet::cv.glmnet(x=X_train, y=nbc, family = "gaussian", type.measure = "mse", alpha = 0.5, keep=T, parallel=T, nfolds=5)
# cat('[INFO] Finished with linear model\n')

# enet_center_multinomial_model <- glmnet::cv.glmnet(x=X_train, y=y_train$binding_counts, family = "multinomial", type.multinomial = "ungrouped", alpha=0.5, keep=T, parallel=T, nfolds=5)
# cat('[INFO] Finished with multinomial model\n')

print('[INFO] Saving the models')

# save the model
saveRDS(enet_center_binary_models_list, file=glue('{project_dir}/models/enet_models/{id}_center_binary_{unique_id}_{run_date}.rds'))
#saveRDS(enet_center_linear_model, file=glue('{project_dir}/models/enet_models/{id}_center_linear_{unique_id}_{run_date}.rds'))
#saveRDS(enet_center_multinomial_model, file=glue('{project_dir}/models/enet_models/{id}_center_multinomial_{unique_id}_{run_date}.rds'))

#registerDoSEQ()
#cat('[INFO] Job completed')