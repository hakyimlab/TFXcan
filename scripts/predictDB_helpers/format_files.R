args <- commandArgs(trailingOnly = TRUE)

file_prefix <- args[1]
output_prefix <- args[2]

library(glue)
library(tidyverse)
library(data.table)

snp_annot <- fread(glue::glue("{file_prefix}.bim")) %>% 
    setnames(.,names(.), c("chr", "snp", "CM", "pos", "alt_vcf", "ref_vcf")) %>%
    dplyr::mutate(rsid = paste(chr, pos, ref_vcf, alt_vcf, 'b37', sep=':')) %>%
    dplyr::mutate(maf = 0.01) %>%  
    dplyr::mutate(varID = str_replace_all(rsid,":","_")) %>%
    dplyr::select(chr, pos, varID, ref_vcf, alt_vcf, maf, rsid)

genotype <- readr::read_table(glue::glue("{file_prefix}.traw"))

genotype <- genotype %>% 
    tidyr::unite('varID', CHR, POS, COUNTED, ALT, sep = '_', remove=FALSE) %>%
    dplyr::mutate(varID = paste(varID, 'b37', sep='_')) %>%
    dplyr::select(-c(CHR,`(C)M`,POS, COUNTED, ALT, SNP)) %>% 
    setnames(.,names(.),gsub("0_", "", colnames(.)))

data.table::fwrite(snp_annot, file=glue("{output_prefix}.snp_annot.txt"), sep='\t', quote=F, row.names=F)
data.table::fwrite(genotype, file=glue("{output_prefix}.geno.txt"), sep='\t', quote=F, row.names=F)

# 1_51479_T_A_b37 1       2       2       2       2       2       1       1       1       2       1       1       2       1       2       1       1       2       1       2    >
# 1_54490_G_A_b37 1       2       2       2       2       2       1       1       1       2       1       1       2       1       2       1       1       2       1       2    >
# 1_54708_G_C_b37 1       0       1       1       1       0       2       1       1       1       1       2       2       2       0       1       1       1       1       2    >
# 1_54716_C_T_b37 1       0       1       1       1       1       2       1       1       1       1       2       2       2       1       1       1       1       1       2    >
# 1_54753_T_G_b37 2