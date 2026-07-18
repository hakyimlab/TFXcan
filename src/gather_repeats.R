# Author: Temi
# Date: Wednesday June 18 2025
# Description: Merges the per-batch flashier .rds outputs from repeat_flash.R into one combined
#   list of flashier runs, for cluster_programs.R to consensus-cluster.
# Usage: Rscript gather_repeats.R --input_basename <path used for repeat_flash.R output_basename> --batch_list <file, one batch name per line> --output_basename <path> --priorL <ebnm_...> --priorF <ebnm_...>

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--input_basename", help='output_basename that was passed to repeat_flash.R; used to glob for its per-batch .rds.gz result files'),
    make_option("--batch_list", help='file listing batch names (one per line) to gather, same file used to drive the repeat_flash.R parallel run'),
    make_option("--output_basename", help='output path prefix for the merged result; written to {output_basename}.{priorL}-{priorF}.{N}Iters.rds.gz'),
    make_option("--priorL", help='flashier EBNM prior name used for the loadings (L) matrix in the upstream repeat_flash.R run; only used here to build the input/output filenames', default="ebnm_point_exponential"),
    make_option("--priorF", help='flashier EBNM prior name used for the factors (F) matrix in the upstream repeat_flash.R run; only used here to build the input/output filenames', default="ebnm_point_exponential")
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

# opt <- list()
# opt$output_basename <- "/beagle3/haky/users/temi/projects/Enpact/data/tenerife/repeat_flash/prca_risk.retrain_repeats"
# opt$batch_list <- "/beagle3/haky/users/temi/projects/Enpact/data/tenerife/repeat_flash/prca_risk.random_subsets.txt"
# opt$priorL <- "ebnm_point_exponential"
# opt$priorF <- "ebnm_point_exponential"

library(glue)
library(data.table)
library(flashier)
library(magrittr)

# read in the batch list
batch_list <- data.table::fread(opt$batch_list, header = FALSE)$V1
# create regex pattern to match the batch names


batch_files <- paste0(opt$input_basename, ".", batch_list, ".*", '.rds.gz') 

batch_files <- Sys.glob(batch_files) 


if(all(file.exists(batch_files))){
    print(glue('INFO - Found {length(batch_files)} files'))
} else {
    stop(glue('ERROR - Not all files exist.'))
}

rds_list <- purrr::map(batch_files, function(x){
    print(glue('INFO - Reading {x}'))
    readRDS(x)
}, .progress = TRUE) %>% unlist(recursive = FALSE)

# save file
output_file <- glue("{opt$output_basename}.{opt$priorL}-{opt$priorF}.{length(rds_list)}Iters.rds.gz")
saveRDS(rds_list, file = output_file, compress = TRUE)
print(glue('INFO - Saved {length(rds_list)} repeats to {output_file}'))