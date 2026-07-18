# Author: Temi
# Date: Tuesday February 24 2026
# Description: apply Enpact weights to one individual's HDF5 features, output per-locus Enpact scores
# Usage: Rscript enpact_predict_hdf5.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--loci_list", help='a list of loci'),
    make_option("--weights", help="Path to Enpact weights file (tsv, with a 'feature' column)"),
    make_option("--individual", help='Individual ID (used to find <input_directory>/<individual>.h5)'),
    make_option("--input_directory", help='Directory containing per-individual HDF5 feature files'),
    make_option("--output_basename", help='Output basename; writes <basename>.enpact_scores.tsv.gz')
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

library(glue)
library(data.table)
library(tidyverse)
library(rhdf5)

# read in the list of loci e.g tss
dt_loci <- data.table::fread(opt$loci_list, header = F)
valid_tss <- dt_loci$V1
# print(length(valid_tss))

# read in the weights 
wgts <- data.table::fread(opt$weights) %>% 
    tibble::column_to_rownames('feature') %>% as.matrix()
# print(dim(wgts))

# collect the tss into epigenome matrix
iFile <- glue::glue("{opt$input_directory}/{opt$individual}.h5")
if(file.exists(iFile)){
    tempout <- purrr::map(valid_tss, function(tss_locus){
        tryCatch({
            exc <- rhdf5::h5read(iFile, glue('{tss_locus}_predictions'))
            epigenome <- t(matrix(rowMeans(exc)))
            return(epigenome)
        }, error = function(err){
            return(matrix(NA, ncol = 5313, nrow = 1)) # 5313 = number of Enformer output tracks
        })
    }, .progress = TRUE) %>% do.call('rbind', .)
    rownames(tempout) <- valid_tss
} else {    
    glue("INFO - {iFile} does not exist")
    base::quit(-1)
}

# print(dim(tempout))
# print(tempout[1:5, 1:5])


# predict using weights: features x weight matrix, plus intercept (row 1 of wgts)
Xpred <- as.matrix(tempout %*% wgts[-1,,drop = FALSE]) + wgts[1,]
# print(Xpred[1:5, 1:5])
dtpred <- Xpred %>% as.data.table(keep.rownames = 'locus')
# print(dtpred[1:5, 1:5])
# save
data.table::fwrite(dtpred, glue("{opt$output_basename}.enpact_scores.tsv.gz"), sep = '\t', quote = F, col.names = T, row.names = F, compress = 'gzip')
