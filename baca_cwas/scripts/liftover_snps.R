arguments <- commandArgs(trailingOnly=TRUE)

library(biomaRt)
library(glue)

snps_file <- arguments[1]
lifted_snps <- arguments[2]

snp_mart <- biomaRt::useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
snp_ids <- scan(input_snps, character(), quote = "", sep=',')
snp_attributes <- c("refsnp_id", 'allele', "chr_name", "chrom_start", "chrom_end")
snp_locations <- biomaRt::getBM(attributes=snp_attributes, filters="snp_filter", values=snp_ids, mart=snp_mart)
write.table(snp_locations, file=lifted_snps, quote=F, row.names=F)