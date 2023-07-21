args <- commandArgs(trailingOnly = TRUE)

output_folder <- dt_dir <- args[1] # '/lus/grand/projects/TFXcan/imlab/data/GEUVADIS/formatted_geno'

#output_folder <- '/lus/grand/projects/TFXcan/imlab/data/GEUVADIS/formatted_geno'

library(glue)
library(tidyverse)
library(data.table)

valid_chr <- paste0('chr', c(1:22), sep='')

#dt_dir <- '/lus/grand/projects/TFXcan/imlab/data/GEUVADIS/formatted_geno'
all_geno_files <- purrr::map(valid_chr, function(chrom){
    ff <- glue('{dt_dir}/ALL.{chrom}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.geno.txt')
    if(file.exists(ff)){
        dt <- data.table::fread(ff) %>% dplyr::mutate(varID=as.character(varID))
        return(dt)
    }
}, .progress=TRUE)

all_snpannot_files <- purrr::map(valid_chr, function(chrom){
    ff <- glue('{dt_dir}/ALL.{chrom}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.snp_annot.txt')
    if(file.exists(ff)){
        dt <- data.table::fread(ff) %>% dplyr::mutate(chr=as.character(chr))
        return(dt)
    }
}, .progress=TRUE)

geno_file <- dplyr::bind_rows(all_geno_files)
snpannot_file <- dplyr::bind_rows(all_snpannot_files)

data.table::fwrite(snpannot_file, file=glue("{output_folder}/all_chrs.snp_annot.txt"), sep='\t', quote=F, row.names=F)
data.table::fwrite(geno_file, file=glue("{output_folder}/all_chrs.geno.txt"), sep='\t', quote=F, row.names=F)