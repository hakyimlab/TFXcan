# Author: Temi
# Date: Thursday July 27 2023
# Description: script to train elastic net TFPred models
# Usage: Rscript train_enet.R [options]

# suppressPackageStartupMessages(library("optparse"))

# option_list <- list(
#     make_option("--train_data_file", help='data to train with enet'),
#     make_option("--features_to_remove", help='', default = NULL),
#     make_option("--nfolds", type="integer", default=5L, help='How many cv folds?'),
#     make_option("--metadata", help='data to train with enet'),
#     make_option("--weights_file", help='data to train with enet'),
#     make_option("--model_rds_file", help='data to train with enet'),
#     make_option("--ncores", default = 12)
# )

# opt <- parse_args(OptionParser(option_list=option_list))

# print(opt)

library(tidyverse)
library(glue)
library(R.utils)
library(data.table)
library(glmnet)
library(doParallel)
library(parallel)
library(foreach)
library(GenomicRanges)


opt <- list()

opt$input_dataframe <- "/beagle3/haky/users/temi/projects/TFXcan/experiments/atac/enpact_tf_tissues.test_path.tsv"
opt$db_metadata <- "/project2/haky/Data/TFXcan/cistrome/raw/human_ca_full_QC.txt"
opt$db_directory <- "/project2/haky/Data/TFXcan/cistrome/raw/human_ca"
opt$output_file <- "/beagle3/haky/users/temi/projects/TFXcan/experiments/atac/atac.bound_vs_unbound.tsv.gz"
# opt$transcription_factor <- 'AR'
# opt$tissue <- 'Prostate'
# opt$test_data <- "/beagle3/haky/users/temi/projects/TFXcan/data/enpact/evaluations/AR_Prostate_2025-04-24.logistic.test_eval.txt.gz"

enpact_metadata <- data.table::fread(opt$input_dataframe)

enpact_metadata <- enpact_metadata %>%
    dplyr::filter(file.info(test_data)$size > 0)


db_metadata <- data.table::fread(opt$db_metadata) 
valid_chromosomes <- paste0('chr', 1:22)

# tissue_case <- stringi::stri_replace_all_regex(
#     qi[['tissue']],
#     "(?<=[a-z])(?=[A-Z])",
#     " "
# )
db_metadata.split <- db_metadata %>% dplyr::filter(Factor == 'ATAC-seq') %>% split(f = .$Tissue_type)
names(db_metadata.split) <- gsub(' ', '', names(db_metadata.split))

find_atac_peaks <- function(db, db_directory, valid_chromosomes){
    use_dcids <- glue::glue("{db_directory}/{db$DCid}_sort_peaks.narrowPeak.bed")
    use_dcids <- use_dcids[file.info(use_dcids)$size > 0]
    atac_cistrome.tissue <- purrr::map(use_dcids, function(ff){
        tryCatch({
            xt <- data.table::fread(ff) %>%
                dplyr::select(chr=V1, start=V2, end=V3, strand = V6) %>% 
                dplyr::filter(chr %in% valid_chromosomes) %>%
                with(., GRanges(chr, IRanges(start, end), strand="*")) %>% GenomicRanges::reduce()
            return(xt)
        }, error = function(err){
            return(err)
        })
    }) 
    atac_cistrome.tissue <- Filter(Negate(is.null), atac_cistrome.tissue)
    atac_cistrome.tissue <- atac_cistrome.tissue %>% as(., "GRangesList") %>% unlist(., recursive = TRUE, use.names = TRUE)
    background.atac <- atac_cistrome.tissue[seqnames(atac_cistrome.tissue) %in% valid_chromosomes]
    seqlevels(background.atac) <- valid_chromosomes
    background.atac <- sort(background.atac) %>% GenomicRanges::reduce()
    return(background.atac)
}




# xx <- db_metadata %>% dplyr::filter(Factor == 'ATAC-seq' & Tissue_type == 'Prostate')
# find_atac_peaks(xx, opt$db_directory, valid_chromosomes)

# glue::glue("{opt$db_directory}/{db_metadata$DCid}_sort_peaks.narrowPeak.bed")

estimate <- function(transcription_factor, tissue, test_data, background){
    atac_cistrome.tissue <- background
    motif_eval <- data.table::fread(test_data) %>% 
        dplyr::select(locus, score = TFPred_score, binding_class) %>% 
        dplyr::mutate(method = case_when(
            binding_class == 1 ~ "bound",
            binding_class == 0 ~ 'unbound'
        )) %>% 
        tidyr::separate_wider_delim(col = locus, delim = '_', names = c('chrom', 'start', 'end'), cols_remove = F) %>%
        dplyr::mutate(across(c(start, end), as.numeric)) %>%
        with(., GRanges(chrom, IRanges(start, end), strand="*", locus, score, binding_class, method))

        # those in accessible regions
    ovalaps <- GenomicRanges::findOverlaps(motif_eval, atac_cistrome.tissue)
    motif.accessible <- motif_eval[queryHits(ovalaps)] %>% as.data.table() %>% dplyr::mutate(feature = 'open')
    motif.inaccessible <- motif_eval[-queryHits(ovalaps)] %>% as.data.table() %>% dplyr::mutate(feature = 'closed')

    cond <- c(nrow(motif.accessible), nrow(motif.inaccessible))
    if(any(cond == 0)){
        return(NULL)
    }

    data_eval <- bind_rows(motif.accessible, motif.inaccessible) %>% dplyr::rename(status = method) %>% 
        dplyr::mutate(group = paste0(feature, status))

    out <- coef(summary(glm(score ~ status + feature, data = data_eval))) %>% 
        as.data.table(keep.rownames = 'term') %>% 
        dplyr::mutate(tf = transcription_factor, tissue = tissue) %>%
        dplyr::relocate(tf, tissue) %>%
        setnames(c('transcription_factor', 'tissue', 'term', 'estimate', 'se', 'tvalue', 'pvalue'))

    return(out)
    
}

num_clusters <- 16 #- 5 # 12 - 5
doParallel::registerDoParallel(num_clusters)

cat('INFO - Registering', num_clusters, 'clusters for a parallel run\n')

evaluations <- foreach::foreach(i=seq_len(nrow(enpact_metadata)), .combine = 'rbind', .inorder=T) %dopar% {
    null_out <- data.table('transcription_factor' = NA, 'tissue' = NA, 'term' = NA, 'estimate' = NA, 'se' = NA, 'tvalue' = NA, 'pvalue' = NA)
    res <- tryCatch({
        qi <- enpact_metadata[i, ]
        dbu <- db_metadata.split[[qi[['tissue']]]]
        bck <- find_atac_peaks(dbu, opt$db_directory, valid_chromosomes)
        out <- estimate(transcription_factor = qi[['transcription_factor']], tissue = qi[['tissue']], test_data = qi[['test_data']], background = bck)
        return(out)
        # if(is.null(out)){
        #     return(null_out)
        # } else {
        #     return(out)
        # }
    }, error = function(err){
        return(null_out)
    })
    return(res)
    
}

doParallel::stopImplicitCluster()

data.table::fwrite(evaluations, opt$output_file, quote = F, compress = 'gzip', col.names = T, row.names = F, sep = '\t')

# db_metadata %>% dplyr::filter(Tissue_type == "Lymph Node")
# db_metadata$Tissue_type |> unique() |> sort()

# enpact_metadata$tissue |> unique()

# stringi::stri_replace_all_regex(
#         "LymphNode",
#         "(?<=[a-z])(?=[A-Z])",
#         " "
#     )

# find_atac_peaks(db_metadata.split[['MammaryGland']], opt$db_directory, valid_chromosomes)


    # if(nrow(dbu) == 0){
    #     return(null_out)
    # } else {
    #     tryCatch({
    #         bck <- find_atac_peaks(dbu, opt$db_directory, valid_chromosomes)
    #         out <- estimate(transcription_factor = qi[['transcription_factor']], tissue = qi[['tissue']], test_data = qi[['test_data']], background = bck)
    #         if(is.null(out)){
    #             return(null_out)
    #         }
    #         return(out)
    #     }, error = function(err){
    #         return(null_out)
    #     })
    # }