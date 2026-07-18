# Author: Temi
# Date: Monday January 26 2026
# Description: like accessibility_effect_on_binding.R, but compares bound vs unbound predicted-binding scores with a
#   Wilcoxon rank-sum test (not a literal t-test, despite the filename) instead of a linear model; run per TF-tissue,
#   once on all sites ("motif" group) and once restricted to accessible ("open") sites, reporting the median score
#   difference alongside the test
# Usage: Rscript accessibility_effect_on_binding.ttest.R (no CLI args -- all inputs/outputs are hardcoded below)

library(tidyverse)
library(glue)
library(R.utils)
library(data.table)
library(glmnet)
library(doParallel)
library(parallel)
library(foreach)
library(GenomicRanges)


# hardcoded inputs/outputs -- not exposed as flags, edit here to point at different data
opt <- list()

opt$input_dataframe <- "/beagle3/haky/users/temi/projects/TFXcan/experiments/atac/enpact_tf_tissues.test_path.tsv"
opt$db_metadata <- "/project2/haky/Data/TFXcan/cistrome/raw/human_ca_full_QC.txt"
opt$db_directory <- "/project2/haky/Data/TFXcan/cistrome/raw/human_ca"
opt$output_file <- "/beagle3/haky/users/temi/projects/TFXcan/experiments/atac/atac.bound_vs_unbound.wilcoxon_tests.dnase_atac.reduced.tsv.gz"

# opt$transcription_factor <- 'AR'
# opt$tissue <- 'Prostate'
# opt$test_data <- "/beagle3/haky/users/temi/projects/TFXcan/data/enpact/evaluations/AR_Prostate_2025-04-24.logistic.test_eval.txt.gz"
# xx <- db_metadata %>% dplyr::filter(Factor == 'ATAC-seq' & Tissue_type == 'Prostate')
# find_atac_peaks(xx, opt$db_directory, valid_chromosomes)
# glue::glue("{opt$db_directory}/{db_metadata$DCid}_sort_peaks.narrowPeak.bed")

# opt$transcription_factor <- 'ATF4'
# opt$tissue <- 'Liver'
# opt$test_data <- "/beagle3/haky/users/temi/projects/TFXcan/data/enpact/evaluations/ATF4_Liver_2024-07-26.logistic.test_eval.txt.gz"
# xx <- db_metadata %>% dplyr::filter(Factor == 'ATAC-seq' & Tissue_type == 'Prostate')
# find_atac_peaks(xx, opt$db_directory, valid_chromosomes)
# glue::glue("{opt$db_directory}/{db_metadata$DCid}_sort_peaks.narrowPeak.bed")


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
# keep only ATAC-seq entries, then split by tissue to build a per-tissue accessible-region background
db_metadata.split <- db_metadata %>% dplyr::filter(Factor == 'ATAC-seq') %>% split(f = .$Tissue_type)
names(db_metadata.split) <- gsub(' ', '', names(db_metadata.split))

find_atac_peaks <- function(db, db_directory, valid_chromosomes){
    use_dcids <- glue::glue("{db_directory}/{db$DCid}_sort_peaks.narrowPeak.bed")
    use_dcids <- use_dcids[file.info(use_dcids)$size > 0]
    atac_cistrome.tissue <- purrr::map(use_dcids, function(ff){
        tryCatch({
            xt <- data.table::fread(ff) %>%
                dplyr::select(chr=V1, start=V2, end=V3, strand = V6) %>% 
                dplyr::filter(chr %in% valid_chromosomes)
            if(nrow(xt) > 0){
                xt <- with(xt, GRanges(chr, IRanges(start, end), strand="*")) %>% GenomicRanges::reduce()
                return(xt)
            } else {
                return(NULL)
            }
                
            
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

# data.table::fread('/project2/haky/Data/TFXcan/cistrome/raw/human_ca/5027_sort_peaks.narrowPeak.bed')

# db <- db_metadata.split[['Prostate']]
# db_directory <- opt$db_directory
# valid_chromosomes



# runs a wilcoxon rank-sum test (bound vs unbound) plus the median score difference, once over all sites
# ("motif") and once restricted to accessible sites ("open")
calculate_differences <- function(input_data){
    motif.pwc <- input_data %>%
        rstatix::wilcox_test(score ~ status) %>%
        rstatix::add_xy_position(x = "status") %>%
        dplyr::mutate(group = 'motif')

    motif.diff <- input_data %>%
        dplyr::group_by(status) %>% 
        dplyr::summarise(med = median(score, na.rm = TRUE)) %>% 
        dplyr::ungroup() %>% 
        tidyr::pivot_wider(names_from = 'status', values_from = med) %>% 
        dplyr::mutate(median_diff = round(bound - unbound, 3)) %>% dplyr::mutate(group = 'motif')

    open.pwc <- input_data %>% 
        dplyr::filter(feature == 'open') %>%
        rstatix::wilcox_test(score ~ status) %>%
        rstatix::add_xy_position(x = "status") %>%
        dplyr::mutate(group = 'open')

    open.diff <- input_data %>% 
        dplyr::filter(feature == 'open') %>%
        dplyr::group_by(status) %>% 
        dplyr::summarise(med = median(score, na.rm = TRUE)) %>% 
        dplyr::ungroup() %>% 
        tidyr::pivot_wider(names_from = 'status', values_from = med) %>% 
        dplyr::mutate(median_diff = round(bound - unbound, 3)) %>% dplyr::mutate(group = 'open')


    open.metrics <- dplyr::inner_join(open.diff, open.pwc, by = 'group')
    motif.metrics <- dplyr::inner_join(motif.diff, motif.pwc, by = 'group')
    data.diff <- dplyr::bind_rows(motif.metrics, open.metrics)
    return(data.diff)
}

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
        dplyr::mutate(group = paste0(feature, status)) %>% distinct()
    out <- calculate_differences(data_eval) %>% 
        dplyr::mutate(transcription_factor = transcription_factor, tissue = tissue) %>%
        dplyr::relocate(transcription_factor, tissue) 

    return(out)
    
}


# bck <- find_atac_peaks(db_metadata.split[['Prostate']], opt$db_directory, valid_chromosomes)
# rest <- estimate('AR', 'Prostate', '/beagle3/haky/users/temi/projects/TFXcan/data/enpact/evaluations/AR_Prostate_2025-04-24.logistic.test_eval.txt.gz', bck)
# calculate_differences(rest)

# rest %>% 
#   group_by(locus) %>% 
#   filter(n()>1)

num_clusters <- 24 #- 5 # 12 - 5
doParallel::registerDoParallel(num_clusters)

cat('INFO - Registering', num_clusters, 'clusters for a parallel run\n')

evaluations <- foreach::foreach(i=seq_len(nrow(enpact_metadata)), .combine = 'rbind', .inorder=T) %dopar% {
    null_out <- data.table('transcription_factor' = NA, 'tissue' = NA, bound = NA, unbound = NA, median_diff = NA, group = NA, '.y.' = NA,  group1 = NA, group2 = NA, n1 = NA, n2 = NA,  statistic= NA, p=NA, y.position =NA, groups = NA, xmin =NA, xmax = NA)
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