# Author: Temi
# Description: Fits a random forest (ranger) predicting one gene's expression from predicted TF
#   binding at its promoter/TSS locus, across individuals. Run once per gene (see
#   promoter_grns.tree.sbatch for the parallel loop over a gene list). Reports out-of-bag R^2 as a
#   proxy for how well TF binding at that locus explains expression — used to build promoter GRNs.
# Usage: Rscript promoter_grns.tree.R --gene <locus id, e.g. chrX_71068333_71068335> --binding_mat <predicted TF binding tsv.gz> --expression_mat <gene expression tsv.gz> --output_file_basename <path prefix>

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--gene", type="character", help='gene/TSS locus id to fit (e.g. chrX_71068333_71068335); must be a row name in --expression_mat'),
    make_option("--binding_mat", type="character", help='predicted TF binding matrix (individuals x TF-tissue features, gzipped tsv) used as predictors'),
    make_option("--expression_mat", type="character", help='gene expression matrix (genes x individuals, gzipped tsv) used as the outcome'),
    make_option("--output_file_basename", type="character", help='output path prefix; result written to {output_file_basename}.{gene}.promoterGRN.GATA1.tsv')
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

library(dplyr)
library(data.table)
library(tidyverse)
library(glue)
library(ranger)

tfxcan_directory <- '/beagle3/haky/users/temi/projects/TFXcan'

# tss_metadata <- data.table::fread('/beagle3/haky/users/temi/projects/dgtex/files/canonical_tss.matched_genes.tsv')

# # Remove duplicate loci
# dups <- tss_metadata %>%
#     dplyr::group_by(locus) %>%
#     dplyr::summarize(cnt = n(), .groups = 'drop') %>%
#     dplyr::filter(cnt > 1)

# dt.tss_metadata <- tss_metadata %>%
#     dplyr::select(ensembl_gene_id, locus, external_gene_name) %>%
#     dplyr::filter(!locus %in% dups$locus) %>%
#     dplyr::distinct()

# cat("Protein-coding genes with unique TSS locus:", nrow(dt.tss_metadata), "\n")

# opt <- list()
# opt$binding_mat="/beagle3/haky/users/temi/projects/TFXcan/experiments/gene_modules/files/predictedTFbinding.GATA1.tsv.gz"
# opt$expression_mat="/beagle3/haky/users/temi/projects/TFXcan/experiments/gene_modules/files/geuvadis.lcl.gene_expression.tsv.gz"
# opt$gene="chrX_71068333_71068335"

# predictors
# file.path(tfxcan_directory, "experiments/gene_modules/files", "predictedTFbinding.GATA1.tsv.gz")
Xpredictors <- data.table::fread(opt$binding_mat) %>% dplyr::select(-tf_tissue)
mat_binding <- tibble::column_to_rownames(as.data.frame(Xpredictors), colnames(Xpredictors)[1])

#Youtcome <- data.table::fread(file.path(gene_modules_dir, 'gene_expression.matrix.tsv'))
# file.path(tfxcan_directory, "experiments/gene_modules/files", "geuvadis.lcl.gene_expression.tsv.gz")
Youtcome <- data.table::fread(opt$expression_mat)
mat_expr <- tibble::column_to_rownames(as.data.frame(Youtcome), colnames(Youtcome)[1])

# Align individuals
shared_inds <- intersect(rownames(mat_binding), colnames(mat_expr))
mat_binding <- mat_binding[shared_inds, ]
mat_expr <- mat_expr[, shared_inds]

cat("Predicted binding matrix dimensions:", nrow(mat_binding), "individuals x", ncol(mat_binding), "TSS locus\n")
cat("Expression matrix dimensions:", nrow(mat_expr), "genes x", ncol(mat_expr), "individuals\n")

# Genes present in both expression matrix and TSS metadata
# genes_to_test <- rownames(mat_expr) #intersect(rownames(mat_expr), dt.tss_metadata$locus)
# cat("Eligible genes:", length(genes_to_test), "\n")

# Sample for speed — set to NULL to run all (slow: ~30s/gene)
# set.seed(42)
# genes_sampled <- sample(genes_to_test, 10)

# fits one gene: regress its expression across individuals on predicted TF binding features
# (random forest via ranger), skip genes with too few complete individuals or no variance
fit_tree <- function(gene, individuals, Xmat, Ymat){

    y_g  <- unlist(Ymat[gene, individuals])
    keep <- !is.na(y_g) & !is.na(rowSums(Xmat))
    if (sum(keep) < 50 || var(y_g[keep]) == 0) return(NULL)

    # print(dim(Xmat))
    # print(length(y_g))

    fit <- ranger(
        x           = Xmat[keep, ],
        y           = y_g[keep],
        importance  = "none",   # skip for speed; enable for GRN edge weights
        num.trees   = 200,
        mtry        = floor(ncol(Xmat) / 10),
        num.threads = 8
    )

    out <- data.frame(
        locus = gene,
        oob_r2          = fit$r.squared,
        n               = sum(keep)
    )
    return(out)
}


results <- tryCatch({
    out <- fit_tree(gene = opt$gene, individuals = shared_inds, Xmat = mat_binding, Ymat = mat_expr)
    out
}, error = function(err){
    out <- data.table(locus = opt$gene, oob_r2 = NA, n = NA)
    out
})

results_rf <- as.data.table(results)

# cat("Genes fitted:", nrow(results_rf), "\n")
# summary(results_rf$oob_r2)

#outfile <- file.path(tfxcan_directory, "experiments/gene_modules/files", "promoterGRN.GATA1.tsv.gz")
outfile <- glue::glue("{opt$output_file_basename}.{opt$gene}.promoterGRN.GATA1.tsv")
data.table::fwrite(results_rf, outfile, row.names = F, col.names = T, quote = F, sep = '\t')

# num of cores
# ncores <- parallelly::availableCores() #detectCores()
# cat("Found ", ncores, " cores...\n")
# num_clusters <- ncores - 4 #- 5 # 12 - 5
# doParallel::registerDoParallel(num_clusters)

# results <- foreach::foreach(i = seq_len(length(genes_to_test)), .inorder = T, .combine = 'rbind', .packages = c("data.table")) %dopar% {
#     res <- tryCatch({
#         out <- fit_tree(gene = genes_to_test[i], individuals = shared_inds, Xmat = mat_binding, Ymat = mat_expr)
#         return(out)
#     }, error = function(err){
#         return(data.table(locus = NA, oob_r2 = NA, n = NA))
#     })
#     return(res)
# } 

# doParallel::stopImplicitCluster()
# print(glue::glue("INFO - Finished running tree fitting"))


# set.seed(2026)
# results_rf <- purrr::map_dfr(genes_to_test, function(g) {
#     y_g  <- unlist(mat_expr[g, shared_inds])
#     keep <- !is.na(y_g) & !is.na(rowSums(mat_binding))
#     if (sum(keep) < 50 || var(y_g[keep]) == 0) return(NULL)

#     fit <- ranger::ranger(
#         x           = mat_binding[keep, ],
#         y           = y_g[keep],
#         importance  = "none",   # skip for speed; enable for GRN edge weights
#         num.trees   = 200,
#         mtry        = floor(ncol(mat_binding) / 10),
#         num.threads = 8
#     )

#     data.frame(
#         locus = g,
#         oob_r2          = fit$r.squared,
#         n               = sum(keep)
#     )
# }, .progress = TRUE)

