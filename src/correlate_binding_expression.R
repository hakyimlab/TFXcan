# Author: Temi
# Date: Tuesday February 24 2026
# Description: correlates a TF-binding-score array (gene x sample x tf_tissue) against a gene-expression matrix,
#   per tf_tissue, via a per-gene linear regression (expression ~ binding)
# Usage: Rscript correlate_binding_expression.R --binding_matrix <rds> --gene_expression <tsv> --output_file_basename <prefix>

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--binding_matrix", help='rds file with the TF-binding-score array (gene x sample x tf_tissue)'),
    make_option("--gene_expression", help='tsv(.gz) gene expression matrix; first column "gene", rest are samples'),
    make_option("--output_file_basename", help='prefix for the output file path'),
    # NOTE: --test is effectively dead -- the actual per-gene test below always runs lm(), regardless of this flag
    make_option("--test", default = "correlation", help='What test? correlation or regression?')
)

opt <- parse_args(OptionParser(option_list=option_list))


# opt <- list()
# opt$binding_matrix <- "/scratch/midway3/temi/enpact_scores/expression_vs_binding.enpact_scores.geuvadis.chippeak_models.data.rds.gz"

# opt$gene_expression <- "/beagle3/haky/users/temi/projects/TFXcan/experiments/gene_modules/files/geuvadis.lcl.gene_expression.tsv.gz"

# opt$output_file_basename <- "/scratch/midway3/temi/enpact_scores/expression_vs_binding.enpact_scores.geuvadis.chippeak_models"

print(opt)

library(glue)
library(data.table)
library(tidyverse)
library(foreach)

# TF binding matrix
mat_TF <- readRDS(opt$binding_matrix)
mat_gexpr <- data.table::fread(opt$gene_expression) %>%
    tibble::column_to_rownames('gene') %>% as.matrix()

print(dim(mat_TF))
print(dim(mat_gexpr))

query <- dimnames(mat_TF)[[3]]
print(length(query))

# correlation or regression

#  lm(MY[i, ] ~ MX[i, ]) %>% summary() %>% coef() %>% as.data.table(keep.rownames = 'term')

#saveRDS(mat_TF[,,1:2], file = glue::glue("{opt$output_file_basename}.data.rds.gz"), compress = TRUE)

# register cores since I want to parallelize
num_clusters <- 64 #- 5 # 12 - 5
doParallel::registerDoParallel(num_clusters)

cat('INFO - Registering', num_clusters, 'clusters for a parallel run\n')

results <- purrr::map(query, function(qq){

    ma <- mat_TF[,,qq]
    mb <- mat_gexpr[,]
    common_genes <- intersect(rownames(ma), rownames(mb))
    common_names <- intersect(colnames(ma), colnames(mb))

    MX <- ma[common_genes, common_names]
    MY <- mb[common_genes, common_names]
    
    out <- foreach::foreach(i = seq_len(nrow(MX)), .inorder = T, .combine = 'rbind') %dopar% {
        res <- tryCatch({
            #test <- cor.test(MX[i, ], MY[i, ], method = 'spearman', exact=FALSE)
            # per-gene regression: expression (MY, gene_expression matrix) ~ binding score (MX, this tf_tissue's slice)
            test <- lm(MY[i, ] ~ MX[i, ]) %>% summary() %>% coef() %>% as.data.table(keep.rownames = 'term') %>% dplyr::mutate(tf_tissue = qq, gene = rownames(MX)[i])
            #return(c(test$p.value, test$estimate))
            return(test)
        }, error = function(err){
            return(c(NA, NA, NA, NA, NA, NA, NA))
        })

        return(res)
    }

    #rownames(out) <- rownames(MX)
    # out <- as.data.table(out) %>%
    #     #stats::setNames(nm = c('gene', 'spearman_pvalue', 'spearman_corr')) %>%
    #     dplyr::mutate(tf = qq)

    return(out)

}, .progress = TRUE)

# test <- with(lt, lm(expr ~ binding)) %>% summary() %>% coef() %>% as.data.table(keep.rownames = 'term')

doParallel::stopImplicitCluster()

results <- Filter(Negate(is.null), results) %>%
    dplyr::bind_rows()

print(head(results))

data.table::fwrite(results, file = glue::glue("{opt$output_file_basename}.regression.tsv.gz"), sep = '\t', col.names = T, row.names = F, compress = 'gzip', quote = F)








# data.table::fwrite(out, file = glue::glue("{opt$output_file_basename}.correlation.tsv.gz"), sep = '\t', col.names = T, row.names = F, compress = 'gzip', quote = F)

# calculate_spearman_correlation_of_two_matrices <- function(ma, mb, matrix_names = c('matrixA', 'matrixB')){
#     # match the names
#     common_genes <- intersect(rownames(ma), rownames(mb))
#     common_names <- intersect(colnames(ma), colnames(mb))

#     MX <- ma[common_genes, common_names]
#     MY <- mb[common_genes, common_names]

#     if(all((nrow(MX) == nrow(MY)) & (ncol(MX) == ncol(MY)))){
#         splitX <- base::asplit(MX, 1)
#         splitY <- base::asplit(MY, 1)

#         test_results <- mapply(function(x1, x2){
#             tryCatch({
#                 test <- cor.test(x1, x2, method = 'spearman', exact=FALSE)
#                 return(c(test$p.value, test$estimate))
#             }, error = function(err){
#                 return(c(NA, NA))
#             })
#         }, splitX, splitY)   
#         dt <- test_results %>% t() %>% 
#             as.data.frame() %>% 
#             stats::setNames(nm = c('spearman_pvalue', 'spearman_r'))
#         return(dt)
#     } else {
#         stop("ERROR - Dimensions or names of input matrices are not the same")
#     }
# }


# calculate_regression_of_two_matrices <- function(ma, mb, matrix_names = c('matrixA', 'matrixB')){
#     # match the names
#     common_genes <- intersect(rownames(ma), rownames(mb))
#     common_names <- intersect(colnames(ma), colnames(mb))

#     MX <- ma[common_genes, common_names]
#     MY <- mb[common_genes, common_names]

#     if(all((nrow(MX) == nrow(MY)) & (ncol(MX) == ncol(MY)))){
#         splitX <- base::asplit(MX, 1)
#         splitY <- base::asplit(MY, 1)
        
#         test_results <- mapply(function(x1, x2){
#             tryCatch({
#                 lt <- data.table(expr = x1, binding = x2)
#                 test <- with(lt, lm(expr ~ binding)) %>% summary() %>% coef() %>% as.data.table(keep.rownames = 'term')
#                 return(test)
#             }, error = function(err){
#                 return(NULL)
#             })
#         }, splitX, splitY, SIMPLIFY = FALSE)

#         dt_results <- Filter(Negate(is.null), test_results) 
        
#         return(dt_results)
#     } else {
#         stop("ERROR - Dimensions or names of input matrices are not the same")
#     }
# }



# out <- purrr::map(query, function(qq){
#     tryCatch({
#         dt <- calculate_regression_of_two_matrices(mat_TF[,,qq], mat_gexpr) %>% dplyr::bind_rows(.id = "gene")
#         #dt <- dt %>% tibble::rownames_to_column('gene')
#         dt$tf_tissue <- qq
#         return(dt)
#     }, error = function(err){
#         return(NULL)
#     })
# }, .progress = TRUE)

# out <- Filter(Negate(is.null), out) %>% dplyr::bind_rows()
# print(out[1:5, ])



# out <- purrr::map(query, function(qq){
#     tryCatch({
#         dt <- calculate_spearman_correlation_of_two_matrices(mat_TF[,,qq], mat_gexpr) %>% dplyr::bind_rows(.id = "gene")
#         #dt <- dt %>% tibble::rownames_to_column('gene')
#         dt$tf_tissue <- qq
#         return(dt)
#     }, error = function(err){
#         return(NULL)
#     })
# }, .progress = TRUE)

# out <- Filter(Negate(is.null), out) %>% dplyr::bind_rows()
# print(out[1:5, ])