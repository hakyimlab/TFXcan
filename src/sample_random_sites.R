# Author: Temi
# Description: builds a bootstrap null distribution for AUPRC by repeatedly resampling (stratified by
#   binding_class) a TF/tissue's real test set, re-extracting Enformer features for each resampled batch,
#   applying the compiled Enpact weights, and computing AUPRC per replicate -- gives a sense of how much AUPRC
#   varies just from which test sites happen to get sampled
# Usage: Rscript sample_random_sites.R (no CLI args -- all inputs/outputs are hardcoded below)

# real inputs/outputs -- hardcoded here (no CLI parsing; a stale optparse block left over from a
# training-data-generation script used to sit here and get silently overwritten, removed since it had no effect)
opt <- list()
opt$transcription_factor <- "AR"
opt$tissue <- "Prostate"
opt$data_file <- glue::glue("/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_CHIPPEAKS_2026-01-02/predictor_files/{opt$transcription_factor}_{opt$tissue}.info.txt.gz")
opt$num_test_observations <- 1000
opt$outputs_directory <- "/scratch/beagle3/temi/random_tests"

print(opt)

seed <- 2026

library(glue)
library(data.table)
library(tidyverse)
library(GenomicRanges)
library(plyranges)
library(purrr)
library(readr)

valid_chromosomes <- paste0('chr', 1:22)

# save 
test_data.dt <- data.table::fread(opt$data_file) %>%
    dplyr::filter(split == 'test') %>%
    dplyr::group_by(binding_class)

set.seed(seed)

# 1000 replicate resamples of the real test set, each stratified by binding_class (num_test_observations sites
# per class) then shuffled -- these are NOT random genomic sites, they're subsamples of the true test labels
permuted_test_set <- replicate(n = 1000, expr = {
    test_set <- test_data.dt %>%
        dplyr::slice_sample(n = opt$num_test_observations) %>%
        dplyr::ungroup() %>%
        dplyr::slice_sample(n = nrow(.)) %>%
        dplyr::mutate(locus = paste(chrom, start, end, sep = '_'))
    return(test_set)
}, simplify = FALSE)


# make sure that they are not all the same

# save in

saveRDS(permuted_test_set, file = glue("{opt$outputs_directory}/{opt$transcription_factor}_{opt$tissue}.data.rds.gz"), compress = TRUE)

# Save each dataframe using its index in the name
out <- purrr::iwalk(permuted_test_set, .f = function(.x, .y){

    ofile <- glue("{opt$outputs_directory}/{opt$transcription_factor}_{opt$tissue}.batch_{.y}.predictors.tsv")

    .x %>%
        dplyr::select(locus) %>%
        data.table::fwrite(glue("{ofile}"), sep = '\t', quote = F, col.names = F, row.names = F)
})


# predict the epigenome

## gather the files and gtools sort
pfiles <- list.files(opt$outputs_directory, pattern = "*.predictors.tsv") %>%
gtools::mixedsort()

out <- purrr::map(seq_along(pfiles), function(i) {

    pfile <- glue("{opt$outputs_directory}/{opt$transcription_factor}_{opt$tissue}.batch_{i}.predictors.tsv")

    efile <- glue("{opt$outputs_directory}/{opt$transcription_factor}_{opt$tissue}.batch_{i}.epigenome.csv.gz")

    # fires off one Enformer-feature-aggregation job per replicate batch and moves on (doesn't wait for it)
    cmd <- glue::glue("sbatch /beagle3/haky/users/temi/projects/TFXcan/src/gather_epigenomes.sbatch {pfile} /beagle3/haky/data/enformer-reference-epigenome {efile}")

    system(cmd)

}, .progress = TRUE)


# merge with ground truth

ground_truth_list <- readRDS(glue("{opt$outputs_directory}/{opt$transcription_factor}_{opt$tissue}.data.rds.gz"))

efiles <- list.files(opt$outputs_directory, pattern = "*.epigenome.csv.gz") %>%
gtools::mixedsort()

merged_data_list <- purrr::map(1:1000, function(i){
    efile <- glue("{opt$outputs_directory}/{opt$transcription_factor}_{opt$tissue}.batch_{i}.epigenome.csv.gz")
    if(file.exists(efile)){
        epigenome <- data.table::fread(efile)
        ground_truth <- ground_truth_list[[i]] %>% dplyr::select(locus, binding_class)
        merged_data <- dplyr::inner_join(ground_truth, epigenome, by = c('locus' = 'id')) %>% dplyr::rename(binding_class = binding_class)
        return(merged_data)
    } else {
        return(NULL)
    }
    
}, .progress = TRUE)

merged_data_list <- Filter(Negate(is.null), merged_data_list)


# X_values_list <- purrr::map(merged_data_list, function(dt){
#     random.test <- dt %>% select(-binding_counts)
#     random.test <- dplyr::distinct(random.test)
#     cc <- complete.cases(random.test)
#     random.test <- random.test[cc, ] %>% dplyr::slice_sample(n = nrow(.))


# })


tfxcan_directory <- '/beagle3/haky/users/temi/projects/TFXcan'

wgts.motifs <- data.table::fread(file.path(tfxcan_directory, 'data/enpact/weights/ENPACT_734_2025-04-24.compiled_weights.with_intercept.lambda.1se.txt.gz')) %>% dplyr::select(feature, motifs_wgts = AR_Prostate) %>% tibble::column_to_rownames('feature') %>% as.matrix()
wgts.motifs[1:5, , drop = F]

predictions <- purrr::map(merged_data_list, function(X_test){

    X_test <- dplyr::distinct(X_test)
    cc <- complete.cases(X_test)
    X_test <- X_test[cc, ] %>% dplyr::slice_sample(n = nrow(.))
    X <- X_test %>% dplyr::select(starts_with("f_")) |> as.matrix()
    y <- X_test %>% dplyr::select(!starts_with("f_")) |> as.data.frame()
    rownames(X) <- y$locus
    # apply the compiled Enpact elastic-net weights (intercept in row 1) to this replicate's features
    Xpred <- as.matrix(X %*% wgts.motifs[-1,,drop = FALSE]) + wgts.motifs[1,]
    Xpred <- Xpred %>% as.data.table(keep.rownames = 'locus') %>% dplyr::inner_join(y, by = 'locus')
    Xpred
}, .progress = TRUE)

calculate_auprc <- function(dt, ycolumn = 'AR_Prostate', labelcolumn = "binding_class", ...){
    dt <- dt[complete.cases(dt), ]
    tryCatch({
        pr <- PRROC::pr.curve(scores.class0 = dt[[ycolumn]], weights.class0=dt[[labelcolumn]], curve = FALSE)
        return(pr$auc.integral)
    }, error = function(e){
        return(NA)
    })
}

auprcs <- sapply(predictions, function(pred){
    auprc <- calculate_auprc(pred, ycolumn = 'motifs_wgts', labelcolumn = "binding_class") %>% round(3)
})

hist(auprcs, ylab = 'Number of replicates')
box()

lower_bound <- quantile(auprcs, 0.025)
upper_bound <- quantile(auprcs, 0.975)
mean_value <- mean(auprcs)
median_value <- median(auprcs)
















# --- everything below is commented-out/dead: an older approach that built the null set from truly random genomic
# sites (peak vs no-peak regions), superseded by the resampling-of-real-test-set approach above. Kept for reference.

# nopeaks.granges <- data.table::fread(glue::glue("{opt$output_directory}/nullpeaks.{opt$transcription_factor}_{opt$tissue}.sorted.bed"), col.names = c('chrom', 'start', 'end', 'id_num', 'width', 'strand', 'score', 'neglog10p', 'neglog10q', 'value')) %>% dplyr::filter(chrom %in% valid_chromosomes) %>%
#     with(GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges(start, end), strand = '*')) %>% GenomicRanges::reduce()

# peaks.granges <- data.table::fread(glue::glue("{opt$output_directory}/peaks.{opt$transcription_factor}_{opt$tissue}.sorted.bed"), col.names = c('chrom', 'start', 'end', 'id_num', 'width', 'strand', 'score', 'neglog10p', 'neglog10q', 'value')) %>% dplyr::filter(chrom %in% valid_chromosomes) %>%
#     with(GenomicRanges::GRanges(seqnames = chrom, ranges = IRanges(start, end), strand = '*')) %>% GenomicRanges::reduce()  

# # find and remove overlaps
# ovlaps <- findOverlaps(query = peaks.granges, subject = nopeaks.granges) 

# if(length(ovlaps) > 0){
#     nopeaks.granges <- nopeaks.granges[-subjectHits(ovlaps)]
#     peaks.granges <- peaks.granges[-queryHits(ovlaps)]
# }


# # add metadata column of bound and unbound
# mcols(nopeaks.granges)$binding_class <- 0
# mcols(peaks.granges)$binding_class <- 1

# training_data.granges <- unique(c(peaks.granges, nopeaks.granges)) %>% GenomicRanges::resize(fix = 'center', width = 512)

# avl_chroms <- seqnames(training_data.granges)@values %>% as.character() %>% unique()

# if(!any(avl_chroms %in% c('chr9', 'chr22'))){
#     print(glue('INFO - There are no training or test observations for {opt$transcription_factor} in {opt$tissue}'))

#     # create empty files
#     columns <- c('chrom', 'start', 'end', 'locus', 'binding_class', 'split')

#     # Create an empty data.table with 0 rows but specified columns
#     empty_dt <- data.table(matrix(nrow = 0, ncol = length(columns)))
#     colnames(empty_dt) <- columns

#     empty_dt %>% data.table::fwrite(glue("{opt$info_file}"), sep = '\t', quote = F, col.names = T, row.names = F)

#     empty_dt %>%
#         dplyr::select(locus, binding_class, split) %>%
#         data.table::fwrite(glue("{opt$ground_truth_file}"), sep = '\t', quote = F, col.names = T, row.names = F)

#     empty_dt %>%
#         dplyr::select(locus) %>% 
#         data.table::fwrite(glue("{opt$predictors_file}"), sep = '\t', quote = F, col.names = F, row.names = F)

#     quit(status = 0)
# }


# training_data.dt <- training_data.granges %>% plyranges::mutate(
#     chrom = as.character(seqnames)
#   ) %>% plyranges::mutate(split =
#     case_when(
#       !chrom %in% c('chr9', 'chr22') ~ "train",
#       chrom %in% c('chr9', 'chr22') ~ "test" # Default case
#     )) %>% as.data.table() %>% dplyr::select(-seqnames) %>% dplyr::relocate(chrom)

# training_data.dt <- training_data.dt %>% dplyr::filter(chrom %in% valid_chromosomes)

# set.seed(seed)
# training_set <- training_data.dt %>% 
#     dplyr::filter(split == 'train') %>%
#     dplyr::group_by(binding_class) %>%
#     dplyr::slice_sample(n = opt$num_train_observations) %>%
#     dplyr::ungroup() %>%
#     dplyr::slice_sample(n = nrow(.)) %>% 
#     dplyr::mutate(locus = paste(chrom, start, end, sep = '_'))

# test_set <- training_data.dt %>% 
#     dplyr::filter(split == 'test') %>%
#     dplyr::group_by(binding_class) %>%
#     dplyr::slice_sample(n = opt$num_test_observations) %>%
#     dplyr::ungroup() %>%
#     dplyr::slice_sample(n = nrow(.)) %>% 
#     dplyr::mutate(locus = paste(chrom, start, end, sep = '_'))

# selected_train <- dplyr::bind_rows(training_set, test_set)
