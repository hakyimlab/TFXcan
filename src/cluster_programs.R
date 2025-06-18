

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--flash_results", help='[input] An rds (list) of multiple iterations of flashier'),
    make_option("--output_basename", help='[output] A dataframe file of the cluster assignments'),
    make_option("--normalization", default='L2', help='[input] The normalization method to use. Default is f (L2 norm). Other options are Linf, L1')
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

library(glue)
library(data.table)
library(flashier)
library(magrittr)
library(factoextra)

# opt <- list(
#     flash_results = '/beagle3/haky/users/temi/projects/Enpact/data/tenerife/prca_flashier_repeats.ebnm_point_exponential-ebnm_point_exponential.1000Iters.rds.gz',
#     output_file = '/beagle3/haky/users/temi/projects/Enpact/data/tenerife/prca_program_clusters.1000.txt.gz',
#     mat_file = '/beagle3/haky/users/temi/projects/Enpact/data/tenerife/prca_program_correlations.1000.rds.gz',
#     loci_assignment = '/beagle3/haky/users/temi/projects/Enpact/data/tenerife/prca_loci_assignments.1000.txt.gz',
#     normalization = 'L2'
# )

# create the output directory
if(!dir.exists(dirname(opt$output_basename))){
    dir.create(dirname(opt$output_basename), recursive = TRUE)
}

# read in the flash results
rFlash_results <- readRDS(opt$flash_results)

print(glue('INFO - Input data has {length(rFlash_results)} flash subsamples'))

norm_type <- switch(opt$normalization,
                    'L2' = 'f',
                    'L1' = 'o',
                    'Linf' = 'i')
# do the ldf decomposition
print(glue('INFO - Performing LDF decomposition with {opt$normalization} norm...'))
rFlash_results <- purrr::map(rFlash_results, function(x){
    xd <- flashier::ldf(x, norm_type) # 'f' for L2 norm, 
    xd$LD <- xd$L %*% diag(xd$D)
    colnames(xd$L) <- colnames(xd$F) <- colnames(xd$LD) <- paste0("k", 1:ncol(xd$F))
    return(xd)
}, .progress = TRUE)

names(rFlash_results) <- paste0('iter', 1:length(rFlash_results))

# collect all the matrices into one big matrix
rownames_order <- rownames(rFlash_results[[1]]$F)
ymat <- lapply(seq_along(rFlash_results), function(i){
    b <- rFlash_results[[i]]$F
    b <- b[rownames_order, ]
    colnames(b) <- glue::glue('iter{i}_k{1:ncol(b)}')
    b
}) %>% do.call('cbind', .) %>% t()


print(glue('INFO - ymat has {nrow(ymat)} programs and {ncol(ymat)} TF/tissue pairs across {length(rFlash_results)} iterations.'))

# calculate correlations
print(glue('INFO - Calculating correlations...'))

n <- nrow(ymat)
# scale and center the data
ymat_minus_mean <- (ymat - rowMeans(ymat)) 
# calculate std dev
ymat_sdev <- sqrt(rowSums(ymat_minus_mean^2)/(n-1))
ymat_scaled <- ymat_minus_mean / ymat_sdev

# cross product i.e t(X) %*% X
yclust <- tcrossprod(ymat_scaled) / (n-1)

print(glue('INFO - Correlation matrix has dimensions: {nrow(yclust)} x {ncol(yclust)}'))
print(glue('INFO - Calculating the best number of clusters...'))
# sample the data
set.seed(2025)
# srows <- sample(1:nrow(yclust), 1000)
# yclust_sample <- yclust[srows, srows]
# sample 100
# srows <- sample(1:length(rFlash_results), 200)
# ldf.100.ylist <- rFlash_results[srows]

# concatenate the Fs
# collect all the matrices into one big matrix
ymat.1000 <- lapply(seq_along(rFlash_results), function(i){
    b <- rFlash_results[[i]]$F
    b <- b[rownames_order, ]
    colnames(b) <- glue::glue('iter{i}_k{1:ncol(b)}')
    b
}) %>% do.call('cbind', .) %>% t()

# find the best cluster using clara
fg_clusters <- factoextra::fviz_nbclust(ymat.1000, cluster::clara, method = "silhouette", k.max = 20, samples = 10, print.summary = TRUE) 

fgd <- fg_clusters$data %>% 
    dplyr::mutate(clusters = as.integer(clusters))

data.table::fwrite(fgd, file = glue::glue('{opt$output_basename}.programs_clara_silhouette.{length(rFlash_results)}samples.{opt$normalization}_norm.txt.gz'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE, compress = 'gzip')

best_k <- base::subset(fgd, fgd$y == max(fgd$y))$clusters

print(glue('INFO - Best number of clusters (k) is {best_k}'))

# cluster the data
print(glue('INFO - Clustering...'))
yclust_pam <- cluster::clara(yclust, k = best_k, metric = "euclidean", stand = FALSE, cluster.only=TRUE, samples = 10)
dt.yclust_pam <- data.table::as.data.table(yclust_pam, keep.rownames = 'subprogram') 
data.table::setnames(dt.yclust_pam, c('subprogram', 'cluster'))

# save the correlations
cluster_matrix_file <- glue::glue('{opt$output_basename}.programs_matrix.{length(rFlash_results)}samples.{opt$normalization}_norm.rds.gz')
saveRDS(ymat, file = cluster_matrix_file, compress = TRUE)

corrmat_file <- glue::glue('{opt$output_basename}.programs_corrmatrix.{length(rFlash_results)}samples.{opt$normalization}_norm.rds.gz')
saveRDS(yclust, file = corrmat_file, compress = TRUE)

# save the cluster assignments
clusters_file <- glue::glue('{opt$output_basename}.program_clusters.{length(rFlash_results)}samples.{opt$normalization}_norm.txt.gz')
data.table::fwrite(dt.yclust_pam, file = clusters_file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE, compress = 'gzip')


# drt <- data.table::fread('/beagle3/haky/users/temi/projects/Enpact/misc/prca_program_clusters.1000.txt.gz')
# yclust_pam <- drt$cluster
# yclust_pam <- setNames(yclust_pam, drt$subprogram)

# rFlash_results <- readRDS('/beagle3/haky/users/temi/projects/Enpact/data/tenerife/prca_flashier_repeats.ebnm_point_exponential-ebnm_point_exponential.1000Iters.rds.gz')
# rFlash_results <- purrr::map(rFlash_results, function(x){
#     xd <- flashier::ldf(x, 'i')
#     colnames(xd$L) <- colnames(xd$F) <- paste0("k", 1:ncol(xd$F))
#     return(xd)
# }, .progress = TRUE)


print(glue('INFO - Assigning loci to clusters'))
loci_assignments <- purrr::map(c(1:best_k), function(i){
    clx <- yclust_pam[which(yclust_pam == i)] %>%
        as.data.frame() %>%
        tibble::rownames_to_column('subprogram') %>%
        stats::setNames(c('subprogram', 'cluster')) %>%
        tidyr::separate(subprogram, c('subsample', 'program'), sep = '_', remove = F)
    
    clusterLs <- purrr::map(1:nrow(clx), function(j){
        subsample <- clx$subsample[j]
        program <- clx$program[j]
        mt <- rFlash_results[[subsample]][['LD']][, program, drop = FALSE]
        mt <- as.data.frame(mt)
        mt <- tibble::rownames_to_column(mt, 'locus')
        return(mt)
    })

    names(clusterLs) <- clx$subsample

    # full_join
    dt <- purrr::reduce(clusterLs, dplyr::full_join, by = 'locus') %>%
        tibble::column_to_rownames('locus') 
    
    dt_loci_clx <- dt %>%
        rowMeans(na.rm = TRUE) %>%
        as.data.frame() %>%
        tibble::rownames_to_column('locus') %>%
        stats::setNames(c('locus', glue('cluster{i}')))

    return(dt_loci_clx)
})

dt_loci_assignments <- purrr::reduce(loci_assignments, dplyr::full_join, by = 'locus')

loci_file <- glue::glue("{opt$output_basename}.loci_assignments.{length(rFlash_results)}samples.{opt$normalization}_norm.txt.gz")

data.table::fwrite(dt_loci_assignments, file = loci_file, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE, compress = 'gzip')