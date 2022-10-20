# This pipeline works with Cistrome data
# Author: Temi
# Date: Circa July 23 2022

setwd('/projects/covid-ct/imlab/users/temi/projects/TFXcan/scripts/')

rm(list=ls())

library(glue)
library(GenomicRanges)
library(reticulate)
library(R.utils)
library(data.table)


# I very much prefer to extract the bed files this way because of the way the files have been archived i.e. zip:gzip...
# and using the python code I wrote makes it easy - also a good way to learn how to use reticulate ^^
reticulate::use_python("~/miniconda3/envs/r-env/bin/python3.10", required = T)
reticulate::py_config()
reticulate::source_python('./load-cistrome-data.py')

py$a <- 'FOXA1' # our TF of interest
pbfile <- py$read_kawakami_data(py$a)

# pick a transcription factor
cell_line <- 'T47D'

# where is the chip-seq data, as well as other directories
cistrome_dir <- '/projects/covid-ct/imlab/data/cistrome/compressed'
homer_dir <- '~/miniconda3/share/homer'
project_dir <- '../'
output_dir <- glue('../processed-data/{TF}')

list.files(cistrome_dir)

# make sure that the TF's Chip data and motif data exist
# check for the motif in Homer
if(!file.exists(glue('{homer_dir}/data/knownTFs/motifs/{tolower(TF)}.motif'))){
    print(glue('No motif information for {TF} available'))
} else {
    print(glue('Motif information for {TF} available. Moving on.'))
}

# check that the chip-seq data is avaialable
# hc_info <- data.table::fread(glue('{cistrome_dir}/human_ca_full_QC.txt'))
# hf_info <- data.table::fread(glue('{cistrome_dir}/human_factor_full_QC.txt'))
# hm_info <- data.table::fread(glue('{cistrome_dir}/human_hm_full_QC.txt'))
hf_info <- data.table::fread(glue('{cistrome_dir}/human_factor_full_QC.txt'))

mtdata_names <- c('Cell_line', 'Cell_type', 'Tissue_type')

# https://genome.ucsc.edu/FAQ/FAQformat.html#format12

if(!TF %in% hf_info$Factor){
    print(glue('No ChIP information for {TF} available'))
} else {
    print(glue('ChIP information for {TF} available. Loading the data...'))

    TF_INFO <- hf_info[hf_info$Factor == TF, c('DCid', 'Cell_line', 'Cell_type', 'Tissue_type')] # across different cell lines and tissues
    py$TF_DCID <- as.integer(TF_INFO$DCid)
    py$b <- 'TF'

    TF_CHIP <- py$get_bed_file(py$TF_DCID, py$b)
    # convert each object to an R object and retain only the first 3 columns
    TF_CHIP <- lapply(seq_along(TF_CHIP), function(i){
        each_chip <- TF_CHIP[[i]]
        each_chip_dcid <- as.numeric(names(TF_CHIP)[i])

        temp <- reticulate::py_to_r(each_chip)
        each_chip_extra_info <- TF_INFO[TF_INFO$DCid == each_chip_dcid, c('Cell_line', 'Cell_type', 'Tissue_type')]
        cbind(temp[, 1:3], each_chip_extra_info)
    })
    TF_CHIP <- do.call(rbind, TF_CHIP)
}

