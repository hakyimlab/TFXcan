suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--binding_matrix", help='a list of loci'),
    make_option("--gene_expression", help='rds file of the splits'),
    make_option("--output_file", help='rds file of the splits')
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

library(glue)
library(data.table)
library(tidyverse)

# TF binding matrix
mat_TF <- readRDS(opt$binding_matrix)
mat_gexpr <- data.table::fread(opt$gene_expression) %>%
    tibble::column_to_rownames('gene') %>% as.matrix()

print(dim(mat_TF))
print(dim(mat_gexpr))

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
#             test <- cor.test(x1, x2, method = 'spearman', exact=FALSE)
#             return(c(test$p.value, test$estimate))
#         }, splitX, splitY)
    
#         dt <- test_results %>% t() %>% 
#             as.data.frame() %>% 
#             stats::setNames(nm = c('spearman_pvalue', 'spearman_r'))
#         return(dt)
#     } else {
#         stop("ERROR - Dimensions or names of input matrices are not the same")
#     }
# }


calculate_regression_of_two_matrices <- function(ma, mb, matrix_names = c('matrixA', 'matrixB')){
    # match the names
    common_genes <- intersect(rownames(ma), rownames(mb))
    common_names <- intersect(colnames(ma), colnames(mb))

    MX <- ma[common_genes, common_names]
    MY <- mb[common_genes, common_names]

    if(all((nrow(MX) == nrow(MY)) & (ncol(MX) == ncol(MY)))){
        splitX <- base::asplit(MX, 1)
        splitY <- base::asplit(MY, 1)
        
        test_results <- mapply(function(x1, x2){
            tryCatch({
                lt <- data.table(expr = x1, binding = x2)
                test <- with(lt, lm(expr ~ binding)) %>% summary() %>% coef() %>% as.data.table(keep.rownames = 'term')
                return(test)
            }, error = function(err){
                return(NULL)
            })
        }, splitX, splitY, SIMPLIFY = FALSE)

        dt_results <- Filter(Negate(is.null), test_results) 
        
        return(dt_results)
    } else {
        stop("ERROR - Dimensions or names of input matrices are not the same")
    }
}

query <- dimnames(mat_TF)[[3]]
print(length(query))

out <- purrr::map(query, function(qq){
    tryCatch({
        dt <- calculate_regression_of_two_matrices(mat_TF[,,qq], mat_gexpr) %>% dplyr::bind_rows(.id = "gene")
        #dt <- dt %>% tibble::rownames_to_column('gene')
        dt$tf_tissue <- qq
        return(dt)
    }, error = function(err){
        return(NULL)
    })
}, .progress = TRUE)

out <- Filter(Negate(is.null), out) %>% dplyr::bind_rows()
print(out[1:5, ])

data.table::fwrite(out, file = opt$output_file, sep = '\t', col.names = T, row.names = F, compress = 'gzip', quote = F)