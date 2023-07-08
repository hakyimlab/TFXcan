


library(data.table)
library(glmnet)
library(glue)


ff <- '/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/compare_predictors/prediction_matrices'

training_loci <- data.table::fread(glue('/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/compare_predictors/metadata/cwas_intervals.txt'),header=FALSE)$V1

tloci <- 'chr4:102674243-102674943'
which_chr <- strsplit(tloci, split=':|-')[[1]]

training_matrix <- data.table::fread(glue('{ff}/{which_chr[1]}/{tloci}_1KG_prediction_matrix.csv'))

X <- training_matrix[, -c(1:2)] |> as.matrix()
y <- training_matrix[, 2] |> unlist() |> unname()

#test_point <- 'chr1:99730260'
test_loci <- colnames(training_matrix)[-c(1:2)][1:100]
fit <- glmnet::cv.glmnet(X, y, alpha = 0.5, family = 'gaussian', type.measure = "mse", keep=T, nfolds=5) # may add : parallel=T, 
betas <- coef(fit, s = "lambda.min") |> as.matrix() # or lambda.1se
betas <- betas[!betas == 0, ][-1] |> as.matrix()
colnames(betas) <- 'beta'

#plot(fit, xvar = "lambda", label = TRUE)



# training_result <- lapply(test_loci, function(each_locus){
#     dt <- cbind(y = y, x = X[, each_locus]) |> as.data.frame()

#     out <- tryCatch({
#         fit <- summary(glm(y ~ x, data=dt)) |> coef()
#         return(fit[2, c(1,2,4)])
#     }, error = function(e){
#         return(matrix(data=NA, nrow = 1, ncol = 3))
#     })
    
# })

training_result <- do.call('rbind', training_result)
rownames(training_result) <- test_loci
colnames(training_result) <- c('beta', 'se', 'pvalue')
training_result

# 

save_dir <- '/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/experiments/compare_predictors/lweights'
data.table::fwrite(training_result, file = glue('{save_dir}/{tloci}_weights.txt'), sep = '\t', row.names=TRUE, quote=FALSE)

##
library(RSQLite)

library(DBI)
# Create an ephemeral in-memory RSQLite database
con <- dbConnect(RSQLite::SQLite(), ":memory:")

dbListTables(con)

load("/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/baca_cwas/data/AR/chr14:94862150-94862750.wgt.RDat", ex <- new.env())


