# Author: Temi
# Date: Tuesday February 24 2026
# Description: collect per-individual Enpact score files into one loci x individual x TF/tissue array (.rds)
# Usage: Rscript collect_predictions_into_rds.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--individuals", help='Path to a list of individual IDs (no header)'),
    make_option("--gene_expression", help="Path to gene expression matrix (genes x individuals, with a 'gene' column)"),
    make_option("--predictions_directory", help='Directory containing per-individual <id>.enpact_scores.tsv.gz files'),
    make_option("--output_file", help='Output path for the combined loci x individual x TF array (.rds)')
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

library(glue)
library(data.table)
library(tidyverse)
library(abind)


dt_ids <- data.table::fread(opt$individuals, header = F)
mat_gexpr <- data.table::fread(opt$gene_expression) %>%
    tibble::column_to_rownames('gene') %>% as.matrix()

print(dim(mat_gexpr))
print(length(dt_ids$V1))

predicted_binding <- purrr::map(dt_ids$V1, function(idd){
    ff <- file.path(opt$predictions_directory, glue("{idd}.enpact_scores.tsv.gz"))
    if(file.exists(ff)){
        mat_tf <- data.table::fread(ff) %>% tibble::column_to_rownames('locus') %>% as.matrix()
        # select from gexpr
        mat_gexpr.id <- mat_gexpr[, idd, drop = F]
        # gene_expression is only used here, to restrict output to TSS with expression data (expression values themselves aren't kept)
        common_tss <- intersect(rownames(mat_gexpr.id), rownames(mat_tf))

        # arrange
        #mat.GEXPR <- mat_gexpr.id[common_tss, , drop = F]
        mat.TFB <- mat_tf[common_tss, ]
        #print(dim(mat.GEXPR)); print(dim(mat.TFB))
        return(mat.TFB)
    } else {
        return(NULL)
    }  
}, .progress = T) 

names(predicted_binding) <- dt_ids$V1
predicted_binding <- Filter(Negate(is.null), predicted_binding)

print(length(predicted_binding))

# combine into an array
myarray <- abind::abind(predicted_binding, along=3)
dimnames(myarray)[[3]] <- names(predicted_binding)
reshapedarray <- aperm(myarray, c(1, 3, 2), resize=TRUE)
lapply(dimnames(reshapedarray), length)
print(reshapedarray[1:3, 1:3, 1:3])
print(dim(reshapedarray))
saveRDS(reshapedarray, file = opt$output_file, compress = "gzip")