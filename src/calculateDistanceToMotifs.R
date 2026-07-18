# Author: Temi
# Date: Wednesday June 18 2025
# Description: for a given TF, finds the distance from each GWAS locus to the nearest HOMER motif instance for that TF
# Usage: Rscript calculateDistanceToMotifs.R --transcription_factor <TF> --output_basename <path>

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--transcription_factor", help='transcription factor name (e.g. AR); used to look up its HOMER motif file'),
    make_option("--output_basename", help='output file prefix; ".distance_to_<TF>_motif.txt.gz" is appended')
)

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table) |> suppressPackageStartupMessages()
library(tidyverse) |> suppressPackageStartupMessages()
library(glue) |> suppressPackageStartupMessages()
library(GenomicRanges)|> suppressPackageStartupMessages()
library(IRanges) |> suppressPackageStartupMessages()
library(plyranges) |> suppressPackageStartupMessages()

glue::glue('INFO - Running for {opt$transcription_factor}...')

# get the loci list
# NOTE: hardcoded input path, not a flag -- this script only works for this one PrCa z-ratio matrix as-is
Y_mat <- readRDS("/beagle3/haky/users/temi/projects/Enpact/data/tenerife/PrCa.tfxcan.zratios.matrix.rds")

gwas_loci <- rownames(Y_mat) %>% strsplit(":") %>% do.call('rbind', .) %>% data.table::as.data.table()

# homer motif files
homer_path <- '/beagle3/haky/users/temi/projects/TFPred-snakemake/data/homer_instances'
tf_motifs <- data.table::fread(file.path(homer_path, opt$transcription_factor, 'merged_motif_file.txt'))

# valid chromosomes
valid_chr <- paste('chr', 1:22, sep = '')

# create the granges objects
granges.tf_motifs <- with(tf_motifs, GenomicRanges::GRanges(seqnames = V2, IRanges::IRanges(start = V3, end = V4), score = V6)) %>%
  plyranges::filter(seqnames %in% valid_chr)

granges.gwas_loci <- with(gwas_loci, GenomicRanges::GRanges(seqnames = V1, IRanges::IRanges(start = as.numeric(V2), end = as.numeric(V2)))) %>%
  plyranges::filter(seqnames %in% valid_chr)

# find the distance to the nearest motif
dd <- IRanges::distanceToNearest(granges.gwas_loci, granges.tf_motifs)

# process and merge the results
aa <- granges.tf_motifs[subjectHits(dd)] %>% as.data.table() %>% dplyr::mutate(motif_site = paste0(seqnames, '_', start, '_', end)) %>% dplyr::select(motif_site, score)

bb <- granges.gwas_loci[queryHits(dd)] %>% as.data.table() %>%
  dplyr::mutate(gwas_loci = paste0(seqnames, '_', start), distance = elementMetadata(dd)$distance) %>% dplyr::select(gwas_loci, distance)

dt <- bind_cols(aa, bb) %>% dplyr::mutate(transcription_factor = opt$transcription_factor)

# save the results
data.table::fwrite(dt, glue::glue("{opt$output_basename}.distance_to_{opt$transcription_factor}_motif.txt.gz"), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE, compress = 'gzip')