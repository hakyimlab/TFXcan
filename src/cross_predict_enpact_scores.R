suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--tf_tissue", help=''),
    make_option("--weights_file", help='.txt file to be created'),
    make_option("--output_basename", help='.txt file to be created'),
    make_option("--input_file", help='.txt file to be created')
)

opt <- parse_args(OptionParser(option_list=option_list))

library(glue)
library(data.table)
library(dplyr)
library(magrittr)
library(glmnet)

# opt <- list()
# opt$tf_tissue <- 'AR_Prostate'
# opt$input_directory <- "/beagle3/haky/users/temi/projects/TFXcan/data/enpact/training/test"
# opt$weights_file <- "/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_CHIPPEAKS_2026-01-02/statistics/ENPACT_CHIPPEAKS_2026-01-02.compiled_weights.lambda.1se.txt.gz"

# tf_tissue <- glue('{opt$transcription_factor}_{opt$tissue}')

if(!dir.exists(dirname(opt$output_basename))){
    dir.create(dirname(opt$output_basename), recursive = TRUE)
}

# read in the test data
#test_file <- Sys.glob(glue::glue("{opt$input_directory}/test_ENPACT*.{opt$tf_tissue}*.csv.gz"))
if(file.exists(opt$input_file)){
    dt.test_file <- data.table::fread(opt$input_file)
    dt.test_file <- dplyr::distinct(dt.test_file)
    # remove missing values
    dt.test_file <- dt.test_file[complete.cases(dt.test_file), ]
    X_test <- dt.test_file %>% dplyr::select(starts_with("f_")) |> as.matrix()
    y_test <- dt.test_file %>% dplyr::select(!starts_with("f_")) |> as.data.frame()
    rownames(X_test) <- y_test$locus
} else {
    print(glue::glue("ERROR - {test_file} does not exist."))
}


# read in the weights
weights <- data.table::fread(opt$weights_file) %>% dplyr::select(all_of(c('feature', opt$tf_tissue)))
weights.matrix <- weights %>% tibble::column_to_rownames('feature') %>% as.matrix()
weights_intercepts <- weights.matrix[1, , drop = F]
weights_features <- weights.matrix[-1, , drop = F]

print(dim(weights))
print(dim(dt.test_file))
# weights[1:5, ]
# dt.test_file[1:5, 1:5]

predictions <- X_test %*% weights_features[, , drop = FALSE]
y_hat <- apply(predictions, 1, function(each_row){
    weights_intercepts + each_row
}) %>% as.matrix()

colnames(y_hat) <- colnames(weights_features)

y_hat <- y_hat[,,drop = FALSE] %>% as.data.table(keep.rownames = 'locus')
predictions <- dplyr::inner_join(y_hat, y_test, by = 'locus')
outfile <- glue("{opt$output_basename}.test.{opt$tf_tissue}.CHIPPEAKS.tsv.gz")
print(glue('INFO - Saving the predicted values to {outfile}'))
data.table::fwrite(predictions, file=glue("{opt$output_basename}.test.{opt$tf_tissue}.CHIPPEAKS.tsv.gz"), sep='\t', row.names=F, col.names=T, quote=F, compress = 'gzip')