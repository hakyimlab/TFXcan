# Author: Temi
# Description: all-pairs spearman correlation between every row of one gene-level matrix and every row of another
#   (e.g. TF-gene vs non-TF-gene expression, for the ACR/gene-module analyses); writes one temp file per row of
#   --binding_matrix, then gathers them all into --final_output_file
# Usage: Rscript correlate_binding_expression.acrs.R --binding_matrix <tsv> --expression_matrix <tsv> --temp_output_basename <prefix> --final_output_file <path>

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--binding_matrix", help='tsv(.gz) matrix, genes x samples, first column "gene_name" -- the row-side of the pairwise correlation'),
    make_option("--expression_matrix", help='tsv(.gz) matrix, genes x samples, first column "gene_name" -- correlated against every row of --binding_matrix'),
    make_option("--temp_output_basename", help='prefix for per-row intermediate correlation files (one per row of --binding_matrix); not deleted afterwards'),
    make_option("--final_output_file", help='path to write the final merged tsv.gz of all pairwise correlation results')
)

opt <- parse_args(OptionParser(option_list=option_list))

# opt <- list()
# tfxcan_directory <- '/beagle3/haky/users/temi/projects/TFXcan'
# opt$binding_matrix <- file.path(tfxcan_directory, "experiments/gene_modules/files", "gene_expression.tf_genes.matrix.tsv.gz")
# opt$expression_matrix <- file.path(tfxcan_directory, "experiments/gene_modules/files", "gene_expression.nontf_genes.matrix.tsv.gz")
# opt$output_basename <- '/scratch/beagle3/temi/gene_modules/correlations/temp'

print(opt)

library(glue)
library(data.table)
library(tidyverse)
library(foreach)

bmat <- data.table::fread(opt$binding_matrix) %>% tibble::column_to_rownames('gene_name') %>% as.matrix()
emat <- data.table::fread(opt$expression_matrix) %>% tibble::column_to_rownames('gene_name') %>% as.matrix()

coln <- intersect(colnames(bmat), colnames(emat))

bmat <- bmat[, coln]
emat <- emat[, coln]

print(dim(bmat))
print(dim(emat))

rn_bmat <- rownames(bmat)
rn_emat <- rownames(emat)

print('Running correlation')

num_clusters <- 16 #- 5 # 12 - 5
doParallel::registerDoParallel(num_clusters)

results <- foreach::foreach(i = seq_len(nrow(bmat)), .inorder = T, .combine = 'rbind', .packages = c("data.table")) %dopar% {
    out <- purrr::map(seq_len(nrow(emat)), function(j){
        #print(glue::glue("Inner: {j}"))
        res <- tryCatch({
            test <- cor.test(bmat[i, ], emat[j, ], method = 'spearman', exact=FALSE)
            out <- data.table(tf_gene = rn_bmat[i], target_gene = rn_emat[j], pvalue = test$p.value, corr_coef = test$estimate)
            return(out)
        }, error = function(err){
            return(data.table(tf_gene = NA, target_gene = NA, pvalue = NA, corr_coef = NA))
        })
    }, .progress = TRUE) %>%
        dplyr::bind_rows()
    #print(glue::glue("Outer: {i}"))
    if(nrow(out) > 0){
        # one file per row i of bmat (all its correlations against every row of emat); gathered back below
        miniFile <- glue("{opt$temp_output_basename}.outer_{i}.correlation.tsv.gz")
        tk <- data.table::fwrite(out, file = miniFile, sep = '\t', col.names = T, row.names = F, compress = 'gzip', quote = F)
        return(miniFile)
    } else {
        return(NULL)
    }
    
} 

doParallel::stopImplicitCluster()
print(glue::glue("INFO - Finished running correlations"))


# gather all
print(glue::glue("INFO - Gathering correlations results into {opt$final_output_file}"))
jumble <- purrr::map(results[, 1], function(ff){
    if(file.exists(ff)){
        data.table::fread(ff)
    } else {
        return(NULL)
    }
}) %>% Filter(Negate(is.null), .) %>% dplyr::bind_rows()

print(head(jumble))
data.table::fwrite(jumble, file = glue("{opt$final_output_file}"), sep = '\t', col.names = T, row.names = F, compress = 'gzip', quote = F)

print(glue::glue("INFO - Gathered correlations results into {opt$final_output_file}"))





# print(head(results))

# data.table::fwrite(results, file = glue("{opt$output_basename}.correlation.tsv.gz"), sep = '\t', col.names = T, row.names = F, compress = 'gzip', quote = F)


# out.pvalues <- apply(bmat, 1, function(col_b){
#   apply(emat, 1, function(col2, col1) {
#     tryCatch({
#       return(cor.test(col2, col1)$p.value)
#     }, error = function(err){
#       return(NA)
#     })
#   }, col1=col_b)
# })

# print('Running correlation')

# out.corr <- apply(bmat, 1, function(col_b){
#   apply(emat, 1, function(col2, col1) {
#     tryCatch({
#       return(cor.test(col2, col1)$estimate)
#     }, error = function(err){
#       return(NA)
#     })
#   }, col1=col_b)
# })

# data.table::fwrite(out.corr, file = glue("{opt$output_basename}.corrcoef.tsv.gz"), sep = '\t', col.names = T, row.names = T, compress = 'gzip', quote = F)


# results <- foreach::foreach(i = seq_len(nrow(bmat))[1:3], j = seq_len(nrow(emat))[1:6], .inorder = T, .combine = 'rbind') %dopar% {
#         res <- tryCatch({
#             test <- cor.test(bmat[i, ], emat[j, ], method = 'spearman', exact=FALSE)
#             return(data.table(tf_gene = rn_bmat[i], target_gene = rn_emat[j], pvalue = test$p.value, corr_coef = test$estimate))
#         }, error = function(err){
#             return(c(NA, NA, NA, NA))
#         })

#         return(res)
#     }




# results <- foreach::foreach(i = seq_len(nrow(bmat))[1], .inorder = T, .combine = 'rbind', .packages = c("data.table")) %:% 
#   foreach::foreach(j = seq_len(nrow(emat))[1:4], .combine = 'rbind') %dopar% {
#         res <- tryCatch({
#             test <- cor.test(bmat[i, ], emat[j, ], method = 'spearman', exact=FALSE)
#             out <- data.table(tf_gene = rn_bmat[i], target_gene = rn_emat[j], pvalue = test$p.value, corr_coef = test$estimate)
#             return(out)
#         }, error = function(err){
#             return(data.table(tf_gene = NA, target_gene = NA, pvalue = NA, corr_coef = NA))
#         })

#         data.table::fwrite(out, file = glue("{opt$output_basename}.correlation.tsv.gz"), sep = '\t', col.names = T, row.names = F, compress = 'gzip', quote = F)
        
#         #return(res)
#     }




# results <- foreach::foreach(i = seq_len(nrow(bmat))[1:3], .inorder = T, .combine = 'rbind', .packages = c("data.table")) %:% 
#   foreach::foreach(j = seq_len(nrow(emat))[1:4], .combine = 'rbind') %dopar% {
#         res <- tryCatch({
#             test <- cor.test(bmat[i, ], emat[j, ], method = 'spearman', exact=FALSE)
#             out <- data.table(tf_gene = rn_bmat[i], target_gene = rn_emat[j], pvalue = test$p.value, corr_coef = test$estimate)
#             return(out)
#         }, error = function(err){
#             return(data.table(tf_gene = NA, target_gene = NA, pvalue = NA, corr_coef = NA))
#         })
#         return(res)
#     }