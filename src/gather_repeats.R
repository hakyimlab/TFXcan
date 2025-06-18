
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--output_basename", help='data preferrably in .rds format of a matrix of GWAS loci by TF/tissue paris of ratios of z-scores'),
    make_option("--batch_list", help='batch name'),
    make_option("--priorL", help='alpha value for enet', default="ebnm_point_exponential"),
    make_option("--priorF", help='number of cores to use', default="ebnm_point_exponential")
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

batch_files <- paste0(opt$output_basename, ".", batch_list, ".*", '.rds.gz') 
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
output_file <- glue("{opt$output_basename}.{opt$priorL}-{opt$priorL}.{length(rds_list)}Iters.rds.gz")
saveRDS(rds_list, file = output_file, compress = TRUE)
print(glue('INFO - Saved {length(rds_list)} repeats to {output_file}'))