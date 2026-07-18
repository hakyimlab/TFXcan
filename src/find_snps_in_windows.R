# Author: Temi
# Date: Monday January 26 2026
# Description: finds SNPs (from per-chromosome dbSNP bed files) that fall inside the significant ARBS/loci windows
#   from the Enpact paper's supplementary table
# Usage: Rscript find_snps_in_windows.R (no CLI args -- all inputs/outputs are hardcoded below)

library(tidyverse)
library(data.table)
library(glue)

tfxcan_directory <- '/beagle3/haky/users/temi/projects/TFXcan'

# NOTE: hardcoded input -- always reads sheet 9 of this specific supplementary-table workbook
tfxcan_significant_results <- readxl::read_excel(file.path(tfxcan_directory, 'results/files', 'enpact-paper.supplementary_tables.xlsx'), sheet = 9, skip = 1)

tfxcan_dt <- tfxcan_significant_results %>% tidyr::separate_wider_delim(arbs, delim = '_',names = c('chr', 'start', 'end'), cols_remove = FALSE) %>% dplyr::mutate(across(c(start, end), as.numeric), chrom = gsub('chr', '', chr) %>% as.numeric) %>% dplyr::select(chrom, chr, start, end, arbs) #%>% dplyr::filter(!is.na(TFXcan.P))

tfxcan_dt.split <- tfxcan_dt %>% base::split(.$chrom) 

out_x <- purrr::map(names(tfxcan_dt.split), function(ch){
    # NOTE: hardcoded SNP bed directory, one gzipped bed file per chromosome
    snp_dataframe <- data.table::fread(glue::glue('/beagle3/haky/data/human_variation/bed/bed_chr_{ch}.bed.gz'), skip = 1) %>% setnames(c('chrom', 'start', 'end', 'rsid', 'score', 'strand'))
    # define the function each time
    find_snps_in_window <- function(wstart, wend){
        found_snps <- snp_dataframe %>% dplyr::filter(data.table::between(start, wstart, wend)) %>% dplyr::select(rsid, pos = start)
        if(nrow(found_snps) > 0){
            return(found_snps)
        } else {return(NULL)}
    }

    data_x <- tfxcan_dt.split[[ch]][, ]
    setDT(data_x)
    out <- data_x[, coefs := Map(find_snps_in_window, start, end)]
    out <- out[, coefs[[1]], by=list(start, end, arbs)]
    return(out)
}, .progress = TRUE)

print(length(out_x))
out_x <- Filter(Negate(is.null), out_x)
print(length(out_x))
out_x_dt <- dplyr::bind_rows(out_x)
print(dim(out_x_dt))
data.table::fwrite(out_x_dt, file.path(tfxcan_directory, 'experiments/coloc', 'prca_cwas.arbs.signif.snps-in-arbs.tsv.gz'), col.names = TRUE, row.names = F, sep = '\t', quote = F, compress = 'gzip')