

args <- commandArgs(trailingOnly = TRUE)
TF <- args[1]
tissue <- args[2]
data_dir <- args[3]

valid_chromosomes <- c(paste('chr', 1:22, sep=''), "chrX")

library(yaml)
library(data.table)
library(glue)
library(tidyverse)
library(GenomicRanges)

output_file <- glue('{data_dir}/predictor_files/{TF}_{tissue}_predictors.txt')

genome_wide_predicted_motifs <- data.table::fread(glue('{data_dir}/homer_files/{TF}/merged_motif_file.txt'))
genome_wide_predicted_motifs <- genome_wide_predicted_motifs 
    %>% dplyr::select(chr=V2, start=V3, end=V4, strand=V5, score=V6) 
    %>% dplyr::filter(chr %in% valid_chromosomes)

# get threshold
threshold <- genome_wide_predicted_motifs 
    %>% dplyr::pull(score) 
    %>% quantile(0.05)
genome_wide_predicted_motifs <- genome_wide_predicted_motifs 
    %>% dplyr::filter(score >= threshold)

# turn into GRanges
tf_motifs_granges <- with(genome_wide_predicted_motifs, GRanges(chr, IRanges(start,end), strand, score))
tf_motifs_granges <- tf_motifs_granges[seqnames(tf_motifs_granges) %in% valid_chromosomes]

# peak files
peak_files_paths <- list.files(glue('{data_dir}/bed_files/{TF}_{tissue}/'), pattern = '*.bed', full.names = TRUE)
peak_files_paths <- peak_files_paths[file.info(peak_files_paths)$size != 0]
peak_files_list <- purrr::map(.x=peak_files_paths, .f=data.table::fread, .progress=T)

pmi_dt_list <- purrr::map(.x=peak_files_paths, function(each_file){
    dt <- data.table::fread(each_file) %>%
        distinct(V1, V2, .keep_all=T) %>%
        dplyr::select(chr=V1, start=V2, end=V3) %>% # select the chr, start and end columns
        with(., GRanges(chr, IRanges(start, end), strand='+', score=0))
    
    dt <- dt[seqnames(dt) %in% valid_chromosomes]

    overlaps <- GenomicRanges::findOverlaps(query=dt, subject=tf_motifs_granges, type='any')

    positive_dt <- tf_motifs_granges[subjectHits(overlaps), ] %>% # because I only want the motifs
        as.data.frame() %>%
        dplyr::select(chr=seqnames, start, end) %>%
        dplyr::mutate(class = 1)

    negative_dt <- tf_motifs_granges[-subjectHits(overlaps), ] %>% # because I only want the motifs
        as.data.frame() %>%
        dplyr::select(chr=seqnames, start, end) %>%
        dplyr::mutate(class = 0)

    return(rbind(positive_dt, negative_dt) |> as.data.frame())

}, .progress=T)

pmi_dt_list <- lapply(seq_along(pmi_dt_list), function(i){
    colnames(pmi_dt_list[[i]])[4] <- paste('class_', i, sep='')
    return(pmi_dt_list[[i]])
})

dt_merged <- pmi_dt_list %>% purrr::reduce(full_join, by = c('chr', 'start', 'end')) 
dt_merged$binding_counts <- rowSums(dt_merged[, -c(1:3)], na.rm=T)
dt_merged$binding_class <- ifelse(dt_merged$binding_counts > 0, 1, 0)
dt_merged <- dt_merged %>%
    dplyr::relocate(c('binding_class', 'binding_counts'), .after=end)

# shuffle the data
set.seed(2023)
dt_merged <- dt_merged[sample(nrow(dt_merged)), ]
dt_merged$chr <- as.character(dt_merged$chr)

# write to file
write.table(dt_merged, file=output_file, sep='\t', quote=F, row.names=F)