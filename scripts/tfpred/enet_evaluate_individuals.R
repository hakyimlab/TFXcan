suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--model_rds", help='The command to run.'),
    make_option("--predict_on", default=NULL, help = "Prefix of files before the extension."),
    make_option("--individuals_data_dir", default=NULL, help = "Prefix of files before the extension."),
    make_option("--output_file", default=NULL, help="Full path plus prefix of output file."),
    make_option("--agg_method", default=NULL, help='Full path to output folder'),
    make_option("--metadata", default=NULL, help='Full path to output folder')
)

opt <- parse_args(OptionParser(option_list=option_list))

#print(arguments)

library(glue)
library(glmnet)
library(data.table)
library(parallel)
library(tidyverse)


#print(glue('model dir is {model_dir}\nTF is {TF}\nmodel_id is {model_id}\noutput directory is {output_dir}\nrun date is {run_date}\nmodel type is {model_type}\n\n'))

# individuals
individuals <- data.table::fread(opt$predict_on, header=F)
individuals <- individuals$V1 |> unique()#[1:5]
print(individuals[1:5])

# read in the model
if(file.exists(opt$model_rds)){
    #print(glue('INFO - Found model at {opt$model_rds}.'))
    model <- readRDS(opt$model_rds)
} else {
    stop('ERROR - model rds file does not exist.')
}

# first gather valid : check that the files exist
valid_names <- c()
for(name in individuals){
    ind_gt <- glue('{opt$individuals_data_dir}/{name}_{opt$agg_method}_{opt$metadata}.csv')
    if(file.exists(ind_gt)){
        valid_names <- append(valid_names, name)
    }
}
print(glue('INFO - Found {length(valid_names)} valid individuals'))

start_time <- proc.time()
ncores <- 32
model_predictions <- parallel::mclapply(valid_names, function(each_ind){
    ind_gt <- glue('{opt$individuals_data_dir}/{each_ind}_{opt$agg_method}_{opt$metadata}.csv')
    if(file.exists(ind_gt)){
        mat_dt <- data.table::fread(ind_gt)
        #print(glue('INFO - Dimension of data is {dim(mat_dt)}'))
    } else {
        #print(glue('INFO - {each_ind} predictors does not exist.'))
        next
    }
    #out <- lapply(agg_method, function(each_method){
    # compute prediction scores
    newx <- as.matrix(mat_dt[, -c(1)])
    locus <- mat_dt[, 1]
    link_pred <- predict(model, newx, s = "lambda.1se", type = 'link') |> as.vector()
    df <- cbind(locus, link_pred) |> as.data.frame() #new_dt[, c(1:2)] |> as.data.frame()
    colnames(df) <- c('locus', 'tfpred_score')
    print(glue('INFO - Dimension of predictions is {dim(df)}'))

    # merge by locus, and class => they should all be the same but this provides a sanity check
    return(df)

}, mc.cores=ncores)

#predictions_list <- list(a=1, b=2)
#names(model_predictions) <- valid_names
out_dt <- purrr::reduce(model_predictions, dplyr::left_join, by=c('locus'))
colnames(out_dt) <- c('locus', valid_names)
#out_dt <- out_dt %>% tibble::column_to_rownames('locus') %>% as.matrix()

end_time <- proc.time()

# save the object to be read later
if(!dir.exists(dirname(opt$output_file))){
    dir.create(dirname(opt$output_file))
} else {
    print(glue('INFO - {dirname(opt$output_file)} exists'))
}

print(glue('INFO - Saving to {opt$output_file}'))
data.table::fwrite(out_dt, file=opt$output_file, quote=FALSE, row.names=FALSE, compress = 'gzip')

summary_dt <- rbind(as.matrix(end_time - start_time), ncores = ncores, nindividuals=length(valid_names))
# save summary
summary_file <- gsub('.txt.gz', '.summary.txt', opt$output_file)
sink(summary_file)
summary_dt
sink()
