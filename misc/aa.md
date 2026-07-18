

```{r}
results_directory <- '/beagle3/haky/users/temi/projects/TFXcan-snakemake/data/prostate_cancer_risk_2025-03-13/factorize'
```



```{r}
outfile <- file.path(results_directory, 'prostate_cancer_risk.prostate_cancer_risk_2025-03-13.program_clusters.txt.gz')
mat_file <- file.path(results_directory, 'prostate_cancer_risk.prostate_cancer_risk_2025-03-13.programs_matrix.rds.gz')
corr_file <- file.path(results_directory, 'prostate_cancer_risk.prostate_cancer_risk_2025-03-13.programs_corrmatrix.rds.gz')
loci_assignment <- file.path(results_directory, 'prostate_cancer_risk.prostate_cancer_risk_2025-03-13.loci_assignments.txt.gz')
```

```{r}
clust1000 <- data.table::fread(outfile)
cluster_assignments <- clust1000$cluster
cluster_assignments <- setNames(cluster_assignments, clust1000$subprogram)
```

```{r}
yclust <- readRDS(mat_file)
dim(yclust); yclust[1:5, 1:5]
```

```{r}
ycorr_mat <- readRDS(corr_file)
dim(ycorr_mat); ycorr_mat[1:5, 1:5]
```


```{r}

nc_programs <- length(table(cluster_assignments))

pclusts <- purrr::map(1:nc_programs, function(i){
    X <- select_ks(cluster_assignments, i, yclust, samp = FALSE)
    X <- rowMeans(X, na.rm = TRUE) |> as.matrix()
    return(X)
}, .progress = TRUE) %>% setNames(paste0('Program_', 1:nc_programs))

prca_tfxcan_programs <- lapply(seq_along(pclusts), function(i){
   pclusts[[i]] %>%
        data.table::as.data.table(keep.rownames = TRUE) %>%
        setNames(c('tf_tissue', 'contribution')) %>%
        dplyr::arrange(desc(contribution)) %>%
        data.table::setnames(c('tf_tissue', glue("Program_{i}")))
}) %>% purrr::reduce(dplyr::inner_join, by = 'tf_tissue')

dt_loci_assignments <- data.table::fread(loci_assignment)
tfxcan_locus_effects <- dt_loci_assignments %>%
    tibble::column_to_rownames('locus') %>%
    setNames(paste0('k', 1:ncol(.))) %>%
    tibble::rownames_to_column('locus')

colnames(tfxcan_locus_effects) <- gsub('k', 'Program_', colnames(tfxcan_locus_effects))
```

```{r}
lapply(pclusts, function(x){
   x %>%
        data.table::as.data.table(keep.rownames = TRUE) %>%
        setNames(c('tf_tissue', 'contribution')) %>%
        dplyr::arrange(desc(contribution)) %>%
        head(n=30)
})
```

```{r}
ep <- prca_tfxcan_programs %>% tibble::column_to_rownames('tf_tissue') %>% as.matrix() %>% apply(2, DescTools::Entropy)
```

```{r}
ep
```


```{r}
data.table::fwrite(prca_tfxcan_programs, file = file.path(files_directory, 'prca.tf_tissue.programs.tsv'), row.names = FALSE, col.names = TRUE, sep = '\t', quote = F)
data.table::fwrite(tfxcan_locus_effects, file = file.path(files_directory, 'prca.loci.assignments.tsv'), row.names = FALSE, col.names = TRUE, sep = '\t', quote = F)
```
