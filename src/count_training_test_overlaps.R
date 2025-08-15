# Usage: Rscript create_training_sets.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--enpact_metadata_file", help='A csv or tsv: transcription_factor, tissue'),
    make_option("--files_directory", help='pattern to match path e.g. /beagle3/haky/users/temi/projects/Enpact/data/enpact/training/train_test'),
    make_option("--output_file")
)

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)
library(foreach)
library(glue)
library(parallel)


# opt <- list()
# opt$enpact_metadata_file <- '/beagle3/haky/users/temi/projects/TFXcan/experiments/auprc/enpact_metadata.734.tsv'
# opt$files_directory <- '/beagle3/haky/users/temi/projects/Enpact/data/enpact/training/train_test'

metadata <- data.table::fread(opt$enpact_metadata_file)

num_clusters <- 2 #- 5 # 12 - 5
cl <- parallel::makeCluster(num_clusters, outfile="/beagle3/haky/users/temi/projects/TFXcan/logs/count_training_test_overlaps.log",type = 'FORK')
doParallel::registerDoParallel(cl)


calculate_overlaps <- function(dt, row_number, files_directory){

    tf <- dt[[row_number, 'transcription_factor']]
    tissue <- dt[[row_number, 'tissue']]
    train_file <- list.files(files_directory, pattern = glue::glue('train_.*.{tf}_{tissue}\\..*.csv.gz'), full.names = TRUE)
    test_file <- list.files(files_directory, pattern = glue::glue('test_.*.{tf}_{tissue}\\..*.csv.gz'),
    full.names = TRUE)

    if(file.exists(train_file) && file.exists(test_file)){
        train_file_size <- file.info(train_file)$size
        test_file_size <- file.info(test_file)$size
        # check that they are not
        condition <- (is.na(train_file_size) || train_file_size == 0L) || ((is.na(test_file_size) || test_file_size == 0L))
        if(isTRUE(condition)){
            return(cbind(transcription_factor = tf, tissue = tissue, overlap = NA))
        } else {
            dt_train <- data.table::fread(train_file)
            dt_test <- data.table::fread(test_file)
            res <- sum(dt_train$locus %in% dt_test$locus)
            return(cbind(transcription_factor = tf, tissue = tissue, overlap = res))
        }
    } else {
        return(cbind(transcription_factor = tf, tissue = tissue, overlap = NA))
    }
}

train_test_overlaps <- foreach::foreach(i=seq_len(nrow(metadata)), .combine='rbind', .inorder=T, .packages = c("glue", "data.table")) %dopar%{
    tryCatch({
        tr_te_overlaps <- calculate_overlaps(dt = metadata, row_number = i, files_directory = opt$files_directory)
        rtf <- tr_te_overlaps[1, 'transcription_factor']
        rtissue <- tr_te_overlaps[1, 'tissue']
        rres <- tr_te_overlaps[1, 'overlap']
        cat(glue('INFO: {rtf}_{rtissue}: {rres}\n\n'))
        return(tr_te_overlaps)
    }, error = function(msg){
            message(glue("ERROR: At index {i}"))
            traceback()
            return(cbind(transcription_factor = metadata[[i, 'transcription_factor']], tissue = metadata[[i, 'tissue']], overlap = NA))
        }
    )
}

data.table::fwrite(train_test_overlaps, file=opt$output_file, sep='\t', row.names=F, col.names=T, quote=F)
parallel::stopCluster(cl)
doParallel::stopImplicitCluster()