




```{r}

pp_above <- purrr::map(top10, function(tt){
    dt.ort_corr[[tt]] %>%
        dplyr::filter(connectivity > 0) %>%
        ggplot(aes(x = abs(zscore), y = connectivity)) +
        geom_point() + theme_simple +
        labs(subtitle = tt)
})

```

```{r}
pp_above <- patchwork::wrap_plots(pp_above, nrow = 3, ncol = 3) + patchwork::plot_layout(axis_titles = "collect")


ggsave(file.path(figures_directory, 'scatterplot.tf_tissue_GRN.positive.correlation.abs.png'), pp_above, width = 6500, height = 6500, dpi = 800, units = 'px')
```


```{r}
ngenes_top10 <- purrr::map(top10, function(tt){
    dt.ort_corr[[tt]] %>% nrow()
})
```

```{r}
pp_corrcoef <- ggplot() +
    geom_histogram(aes(ort_corr$corr), col = 'black') +
    theme_simple +
    labs(y = 'Number of TFs', x = 'Correlation Coefficient')

ggsave(file.path(figures_directory, 'histogram.tf_tissue_GRN.correlation.png'), pp_corrcoef, width = 3000, height = 3000, dpi = 800, units = 'px')
```

```{r}
dt_all <- dplyr::bind_rows(dt.ort_corr, .id = 'tf_tissue')
```

```{r}
mat_LCL_grn['FOXA1', ] %>% as.numeric() %>% hist()
```

```{r}
ppgrn <- ggplot() +
    geom_histogram(aes(as.numeric(mat_LCL_grn['FOXA1', ])), col = 'black') +
    theme_simple +
    labs(x = 'Connectivity', y = 'Number of Genes', subtitle = 'FOXA1 Gene Regulatory Network')

ggsave(file.path(figures_directory, 'histogram.FOXA1_GRN.png'), ppgrn, width = 3000, height = 3000, dpi = 800, units = 'px')
```

```{r}
out_gene <- dt_all %>% 
    dplyr::group_by(gene) %>% 
    dplyr::group_map(.f = function(.x, .y){
        dt <- .x %>% 
            tidyr::separate_wider_delim(col = tf_tissue, names = c('tf', 'tissue'), delim = '_', cols_remove = F) %>%
            dplyr::group_by(tf) %>%
            dplyr::summarize(sumZ = sum(zscore), meanC = mean(connectivity))
        ct <- with(dt, cor.test(sumZ, meanC))
        return(data.table(tf_tissue = .y, pvalue = ct$p.value, corr = ct$estimate))
    }) 
```

```{r}
dt.out_gene <- dplyr::bind_rows(out_gene) %>%
    dplyr::mutate(fdrsignif = qvalue::qvalue(pvalue, fdr.level = 0.05)$significant)
```

```{r}
dt.out_gene$pvalue |> hist()
```

```{r}
dt.out_gene %>% dplyr::filter(fdrsignif == TRUE)
```

```{r}
dt_blood_grn %>% dplyr::filter(V1 == 'CTCF')
```

```{r}
ctcf_blood_corr <- corr.expr_binding %>% dplyr::filter(tf_tissue == 'CTCF_Blood') %>% dplyr::inner_join(tmed, by = c('gene' = 'locus'))
ctcf_blood_corr[1:5, 1:5]
```

```{r}
cm_genes <- intersect(colnames(dt_blood_grn), ctcf_blood_corr$ensembl_gene_id)
length(cm_genes)
```

```{r}
tf_expr_aa <- ctcf_blood_corr %>%
    dplyr::filter(ensembl_gene_id %in% cm_genes) %>% dplyr::select(gene = ensembl_gene_id, spearman_r)

network_aa <- dt_blood_grn %>% dplyr::select(all_of(c('V1', cm_genes))) %>% dplyr::filter(V1 == 'CTCF') %>% data.table::transpose(., keep.names = "V1", make.names = "V1") %>% setnames(c('gene', 'connectivity'))
```

```{r}
network_aa |> dim(); dim(tf_expr_aa)
```

```{r}
mgd <- dplyr::inner_join(tf_expr_aa, network_aa, by = 'gene')
```

```{r}
with(mgd, cor.test(DescTools::FisherZ(spearman_r), connectivity))
```

```{r}

```

```{r}
t(network_aa) |> head()
```




```{r}
knitr::knit_exit()
```

```{r}
idvalid <- c()
predicted_binding <- purrr::map(dt_ids$individual[1:3], function(idd){
    ff <- file.path(predictions_directory, glue("{idd}.enpact_scores.tsv.gz"))
    if(file.exists(ff)){
        idvalid <<- append(idvalid, idd)
        mat_tf <- data.table::fread(ff) %>% tibble::column_to_rownames('locus') %>% as.matrix()
        # select from gexpr
        mat_gexpr.id <- mat_gexpr[, idd, drop = F]
        # find common genes
        common_tss <- intersect(rownames(mat_gexpr.id), rownames(mat_tf))

        # arrange
        #mat.GEXPR <- mat_gexpr.id[common_tss, , drop = F]
        mat.TFB <- mat_tf[common_tss, ]
        #print(dim(mat.GEXPR)); print(dim(mat.TFB))
        return(mat.TFB)
    } else {
        return(NULL)
    }  
}, .progress = T) %>% Filter(Negate(is.null), .)
```

```{r}
length(predicted_binding); predicted_binding[[1]][1:5, 1:5]
```

```{r}
library(abind)
```

```{r}
#combine into an array
myarray <- abind::abind(predicted_binding, along=3)
dimnames(myarray)[[3]] <- idvalid

reshapedarray <- aperm(myarray, c(1, 3, 2), resize=TRUE)
reshapedarray[1:5, , 1:4]; dim(reshapedarray)
# saveRDS(reshapedarray, file = "/scratch/midway3/temi/enpact_scores/enpact_scores.geuvadis.rds.gz", compress = "gzip")

# # process and prepare for predictDB; create a predictDB input EnpactScores file
# tarray <- as.array(reshapedarray)
# tarray <- aperm(tarray, c(3, 1, 2), resize=TRUE)
# dta <- dim(tarray)
# farray <- as.array(tarray)
# dim(farray) <- c(dta[1] * dta[2], dta[3])
```




## Correlation with gene expression
```{r}
dim(mat_gexpr)
```

```{r}

```
```{r}
calculate_spearman_correlation_of_two_matrices <- function(ma, mb, matrix_names = c('matrixA', 'matrixB')){
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
                return(NA)
            })
        }, splitX, splitY, SIMPLIFY = FALSE)

        dt_results <- Filter(Negate(is.null), test_results) 
        
        return(dt_results)
    } else {
        stop("ERROR - Dimensions or names of input matrices are not the same")
    }
}
```

```{r}
out <- calculate_spearman_correlation_of_two_matrices(reshapedarray[1:5, 1:3, 'CTCF_Brain'], mat_gexpr[1:5, 1:3]) %>% 
    dplyr::bind_rows(.id = "locus")
out
```

```{r}
out %>% dplyr::bind_rows(.id = "locus")
```


```{r}
firstOne <- base::asplit(reshapedarray[1:5, 1:3, 'CTCF_Brain'], 1)
secondOne <- base::asplit(mat_gexpr[1:5, 1:3], 1)
```

```{r}
firstOne
```



```{r}
sapply(out, coef)
```

```{r}
out |> dim()
```

```{r}
head(out, n = 20)
```






















```{r}
DTH %>% dplyr::filter(tf_tissue == 'AR_Prostate') %>% dplyr::pull(spearman_pvalue) %>% hist()
```


```{r}
# correlation
    out <- apply(mat.TFB, 2, function(ea){
        ctest <- cor.test(ea, mat.GEXPR[, 1])
        data.table(pvalue = ctest$p.value, corr = ctest$estimate)
    }) %>% dplyr::bind_rows() %>% dplyr::mutate(locus = colnames(mat.TFB))
```

```{r}
predicted_binding[[3]] %>% summary()
```

```{r}
tempout <- purrr::map(valid_tss, function(tss_locus){
    tryCatch({
        exc <- rhdf5::h5read("/beagle3/haky/data/Geuvadis_TSS_enformer_prediction/NA20783.h5", glue('{tss_locus}_predictions'))
        epigenome <- t(matrix(rowMeans(exc)))
        return(epigenome)
    }, error = function(err){
        return(matrix(NA, ncol = 5313, nrow = 1))
    })
    
}, .progress = TRUE) %>% do.call('rbind', .)

rownames(tempout) <- valid_tss
```

```{r}
tempout[1:3, 1:5]; dim(tempout)
```

```{r}
wgts <- data.table::fread(file.path('/beagle3/haky/users/temi/projects/TFXcan', 'data/enpact/weights/ENPACT_734_2025-04-24.compiled_weights.with_intercept.lambda.1se.txt.gz')) %>% tibble::column_to_rownames('feature') %>% as.matrix()

Xpred <- as.matrix(tempout %*% wgts[-1,,drop = FALSE]) + wgts[1,]
```

```{r}
Xpred[1:5, 1:5]
```


```{r}
data_dir <- "/beagle3/haky/users/temi/data/geuvadis_promoters"
```

```{r}
geuvadis_epigenome <- data.table::fread(file.path(data_dir, 'NA20828_aggByMean_AR_Prostate.csv'))
```


```{r}
dim(geuvadis_epigenome); geuvadis_epigenome[1:5, 1:5]
```








```{r}
dt_prca_grn <- data.table::fread(file.path(grn_directory, 'panda_tcga_prad.csv'))
dt_prostate_grn <- data.table::fread(file.path(grn_directory, 'Prostate.csv'))
dim(dt_prca_grn); dim(dt_prostate_grn)
```

```{r}
# match for common TFs
common_tfs <- intersect(dt_prca_grn$tf, dt_prostate_grn$V1)
common_genes <- intersect(colnames(dt_prca_grn), colnames(dt_prostate_grn))
length(common_tfs); length(common_genes)
```

```{r}
tf_grn.prca.mat <- dt_prca_grn %>% 
    tibble::column_to_rownames('tf') %>% 
    as.matrix()

tf_grn.prostate.mat <- dt_prostate_grn %>% 
    tibble::column_to_rownames('V1') %>% 
    as.matrix()

tf_grn.prca.mat <- tf_grn.prca.mat[common_tfs, common_genes]
tf_grn.prostate.mat <- tf_grn.prostate.mat[common_tfs, common_genes]
```

```{r}
dim(tf_grn.prca.mat); dim(tf_grn.prostate.mat)
```

997 TFs to 21316 genes. The higher the edge weights, the stronger the regulation. 


```{r}
mtplot <- ComplexHeatmap::Heatmap(tf_grn.prca.mat[1:20, 1:20], col = brewer.pal(n = 8, name = "Reds"), border = 'black', cluster_rows = F, cluster_columns = F, column_names_rot = 90, heatmap_legend_param = list(title = "Edge\nWeight"), row_names_gp = gpar(fontsize = 7, fontfamily = 'Helvetica'), column_names_gp = gpar(fontsize = 7, fontfamily = 'Helvetica'), rect_gp = gpar(col = "black", lwd = 1)) %>% ComplexHeatmap::draw() %>%
    grid.grabExpr() %>%
    cowplot::ggdraw()

ggsave(filename=file.path(figures_directory, "prostate_cancer.grn.few_example.heatmap.png"), mtplot, width = 4000, height = 3000, units = "px", dpi = 800, limitsize = FALSE)
```

```{r}
mtplot <- ComplexHeatmap::Heatmap(tf_grn.prostate.mat[1:20, 1:20], col = brewer.pal(n = 8, name = "Reds"), border = 'black', cluster_rows = F, cluster_columns = F, column_names_rot = 90, heatmap_legend_param = list(title = "Edge\nWeight"), row_names_gp = gpar(fontsize = 7, fontfamily = 'Helvetica'), column_names_gp = gpar(fontsize = 7, fontfamily = 'Helvetica'), rect_gp = gpar(col = "black", lwd = 1)) %>% ComplexHeatmap::draw() %>%
    grid.grabExpr() %>%
    cowplot::ggdraw()

ggsave(filename=file.path(figures_directory, "prostate_tissue.grn.few_example.heatmap.png"), mtplot, width = 4000, height = 3000, units = "px", dpi = 800, limitsize = FALSE)
```

### How many genes on average does a TF regulate?

```{r}
tfxcan_grn.mat <- dt_grn %>% 
    dplyr::filter(tf %in% tfxcan_tfs) %>% 
    tibble::column_to_rownames('tf') %>% 
    as.matrix()

dim(tfxcan_grn.mat); tfxcan_grn.mat[1:5, 1:5]
```


Five point summary statistic
```{r}
tf_grn.mat <- dt_grn %>% 
    tibble::column_to_rownames('tf') %>% 
    as.matrix()
tf_grn.mat[1:5, 1:5]; dim(tf_grn.mat )
```

```{r}
pp1.AR <- ggplot() +
    geom_density(aes(x = tf_grn.prca.mat['AR', ], fill = 'Prostate Cancer'), alpha = 0.3, col = 'black') +
    geom_density(aes(x = tf_grn.prostate.mat['AR', ], fill = 'Normal Prostate'), alpha = 0.3, col = 'black') +
    scale_fill_manual(name ='Tissue Source', values = c('Prostate Cancer' = 'red', 'Normal Prostate' = 'green')) +
    labs(subtitle = 'Distribution of AR edge weights in prostate tissue', x = 'Strength of Regulation (Edge Weights)', y = 'Number of Genes') +
    theme_minimal() +
    theme(legend.position="bottom")

```

```{r}
pp1.AR
```

### How can I have a global edge weight cut off?

```{r}
median_edge_weights.prca <- apply(tf_grn.prca.mat, 1, median)
median_edge_weights.prostate <- apply(tf_grn.prostate.mat, 1, median)
summary(median_edge_weights.prca); summary(median_edge_weights.prostate)
```

```{r}
xx.prca <- as.data.table(median_edge_weights.prca, keep.rownames = 'tf')
xx.prostate <- as.data.table(median_edge_weights.prostate, keep.rownames = 'tf')

xx.merged <- dplyr::inner_join(xx.prca, xx.prostate, by = 'tf') %>% tidyr::pivot_longer(cols = !c(tf), names_to = 'dataset', values_to = 'edge_weight') 
```

```{r}
pp1.ALL_TF.boxplot <- ggplot(xx.merged, aes(y = edge_weight, fill = dataset)) +
    geom_boxplot(alpha = 0.3, col = 'black') +
    theme_minimal() +
    scale_fill_manual(name ='Tissue Source', values = c('median_edge_weights.prca' = 'red', 'median_edge_weights.prostate' = 'green'), labels = c('Cancerous Prostate', 'Normal Prostate')) +
    labs(subtitle = 'Boxplot of median edge weights in prostate tissues', y = 'Regulations Strength (Edge Weights)', y = '') +
    theme(legend.position="bottom")

ggsave(file.path(figures_directory, "boxplot.edge_weights.medians.pdf"), pp1.ALL_TF.boxplot, width = 4000, height = 4000, dpi = 800, units = 'px', limitsize = F)
```



## maximum
```{r}
max_edge_weights.prca <- apply(tf_grn.prca.mat, 1, max)
max_edge_weights.prostate <- apply(tf_grn.prostate.mat, 1, max)
summary(max_edge_weights.prca); summary(max_edge_weights.prostate)
```

```{r}
xx.prca <- as.data.table(max_edge_weights.prca, keep.rownames = 'tf')
xx.prostate <- as.data.table(max_edge_weights.prostate, keep.rownames = 'tf')

xx.merged <- dplyr::inner_join(xx.prca, xx.prostate, by = 'tf') %>% tidyr::pivot_longer(cols = !c(tf), names_to = 'dataset', values_to = 'edge_weight') 
```

```{r}
pp1.ALL_TF.boxplot <- ggplot(xx.merged, aes(y = edge_weight, fill = dataset)) +
    geom_boxplot(alpha = 0.3, col = 'black') +
    theme_minimal() +
    scale_fill_manual(name ='Tissue Source', values = c('max_edge_weights.prca' = 'red', 'max_edge_weights.prostate' = 'green'), labels = c('Cancerous Prostate', 'Normal Prostate')) +
    labs(subtitle = 'Boxplot of maximum edge weights in prostate tissues', y = 'Regulations Strength (Edge Weights)', y = '') +
    theme(legend.position="bottom")

ggsave(file.path(figures_directory, "boxplot.edge_weights.maximums.pdf"), pp1.ALL_TF.boxplot, width = 4000, height = 4000, dpi = 800, units = 'px', limitsize = F)
```


## minimum
```{r}
min_edge_weights.prca <- apply(tf_grn.prca.mat, 1, min)
min_edge_weights.prostate <- apply(tf_grn.prostate.mat, 1, min)
summary(min_edge_weights.prca); summary(min_edge_weights.prostate)
```

```{r}
xx.prca <- as.data.table(min_edge_weights.prca, keep.rownames = 'tf')
xx.prostate <- as.data.table(min_edge_weights.prostate, keep.rownames = 'tf')

xx.merged <- dplyr::inner_join(xx.prca, xx.prostate, by = 'tf') %>% tidyr::pivot_longer(cols = !c(tf), names_to = 'dataset', values_to = 'edge_weight') 
```

```{r}
pp1.ALL_TF.boxplot <- ggplot(xx.merged, aes(y = edge_weight, fill = dataset)) +
    geom_boxplot(alpha = 0.3, col = 'black') +
    theme_minimal() +
    scale_fill_manual(name ='Tissue Source', values = c('min_edge_weights.prca' = 'red', 'min_edge_weights.prostate' = 'green'), labels = c('Cancerous Prostate', 'Normal Prostate')) +
    labs(subtitle = 'Boxplot of minimum edge weights in prostate tissues', y = 'Regulations Strength (Edge Weights)', y = '') +
    theme(legend.position="bottom")

ggsave(file.path(figures_directory, "boxplot.edge_weights.minimum.pdf"), pp1.ALL_TF.boxplot, width = 4000, height = 4000, dpi = 800, units = 'px', limitsize = F)
```



```{r}
pp1.ALL_TF <- ggplot() +
    geom_boxplot(aes(x = median_edge_weights.prca, fill = 'Prostate Cancer'), alpha = 0.3, col = 'black') +
    geom_boxplot(aes(x = median_edge_weights.prostate, fill = 'Normal Prostate'), alpha = 0.3, col = 'black') +
    scale_fill_manual(name ='Tissue Source', values = c('Prostate Cancer' = 'red', 'Normal Prostate' = 'green')) +
    labs(subtitle = 'Distribution of median edge weights in prostate tissue', x = 'Median Strength of Regulation (Edge Weights)', y = 'Number of Genes') +
    theme_minimal() +
    theme(legend.position="bottom")

```

```{r}
pp1.ALL_TF
```

### Determinining a global edge weight cutoff

```{r}
# how many genes are greater than the average median weight
num_genes_regulated.prca <- apply(tf_grn.prca.mat, 1, function(err){
    sum(err > median(median_edge_weights.prca))
})

num_genes_regulated.prostate <- apply(tf_grn.prostate.mat, 1, function(err){
    sum(err > median(median_edge_weights.prostate))
})

```

```{r}
xx.prca <- as.data.table(num_genes_regulated.prca, keep.rownames = 'tf')
xx.prostate <- as.data.table(num_genes_regulated.prostate, keep.rownames = 'tf')

xx.merged <- dplyr::inner_join(xx.prca, xx.prostate, by = 'tf')
```

```{r}
pp1.num_regulated <- ggplot(xx.merged, aes(x = num_genes_regulated.prostate, y = num_genes_regulated.prca)) +
    geom_point() +
    theme_minimal() +
    geom_abline(intercept = 0, slope = 1, col = 'grey') +
    labs(x = 'Number of Genes Regulate (Normal Prostate)', y = 'Number of Genes Regulated (Prostate Cancer)') +
    theme(legend.position="bottom")
```

```{r}
summary(xx.merged$num_genes_regulated.prca);
summary(xx.merged$num_genes_regulated.prostate)
```

```{r}
ppn <- patchwork::wrap_plots(pp1.AR, pp1.ALL_TF, pp1.num_regulated, nrow = 1)
ggsave(file.path(figures_directory, "number_genes_regulated.medians.png"), ppn, width = 10000, height = 4500, dpi = 800, units = 'px', limitsize = F)
```


```{r}
ppn <- patchwork::wrap_plots(pp1.AR, pp1.ALL_TF, nrow = 1)
ggsave(file.path(figures_directory, "edge_weight_distribution.png"), ppn, width = 9000, height = 4500, dpi = 800, units = 'px', limitsize = F)
```

```{r}

```

### Connecting to TFXcan programs

```{r}
consensus_tfxcan.F <- readxl::read_excel(file.path(tfxcan_directory, 'results/files/enpact-paper.supplementary_tables.xlsx'), sheet = 'Table S12', skip = 1)

dt.consensusF <- consensus_tfxcan.F %>%
    tidyr::pivot_longer(
        cols = !tf_tissue, 
        names_to = "program", 
        values_to = "contribution"
  )
```

```{r}
tfxcan_tfs <- dt.consensusF %>% dplyr::select(tf_tissue) %>% distinct() %>% tidyr::separate_wider_delim(col = tf_tissue, names = c('tf', 'tissue'), delim = '_', cols_remove = F) %>% pull(tf) %>% unique()
```

### Using one program as an example

```{r}
prgs <- unique(dt.consensusF$program)
```

```{r}
prgs_goa <- purrr::map(prgs, function(prg){
    print(prg)
    if(prg %in% c('CTCF', 'REST')){
        return(NULL)
    }
    top10.drivers <- dt.consensusF %>% 
        dplyr::filter(program == prg) %>% 
        dplyr::arrange(desc(contribution)) %>% head(n = 10) %>% 
        tidyr::separate_wider_delim(col = tf_tissue, names = c('tf', 'tissue'), delim = '_', cols_remove = F) %>% 
        pull(tf) %>% 
        unique()

    weights.prca.drivers <- tf_grn.prca.mat[rownames(tf_grn.prca.mat) %in% top10.drivers, ] 
    driver_regulated_genes <- apply(weights.prca.drivers, 1, function(orr){
        ev <- orr[orr >= median(median_edge_weights.prca)]
        return(ev)
    }, simplify = FALSE)
    genes_drivers <- purrr::reduce(lapply(driver_regulated_genes, names), base::intersect)

    if(length(genes_drivers) == 0){
        return(NULL)
    }

    return(genes_drivers)
}, .progress = TRUE)

names(prgs_goa) <- prgs
prgs_goa <- Filter(Negate(is.null), prgs_goa)
```

```{r}
ora_res <- purrr::map(prgs_goa, function(ggs){
    res <- gprofiler2::gost(query = ggs,
            organism = "hsapiens", ordered_query = FALSE,
            multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
            measure_underrepresentation = F, evcodes = FALSE, 
            user_threshold = 0.05, correction_method = "fdr", custom_bg = colnames(tf_grn.prca.mat), numeric_ns = "", sources = c("GO:MF", "REAC", 'GO:BP', 'GO:CC', 'KEGG'), as_short_link = FALSE, highlight = TRUE)
    return(res$result)
})

ora_res <- Filter(Negate(is.null), ora_res)

```

```{r}
length(ora_res)
```

```{r}
lapply(ora_res, dim)
```

```{r}
prgs_goa <- Filter(Negate(is.null), prgs_goa)
```

```{r}
prgs_goa$`Hormone Receptors` |> View()
```

```{r}
prgs_goa$`Colon` |> View()
```

```{r}
prgs_goa[['Hormone Receptors']]
```

```{r}
enrichment_plots <- purrr::map(names(ora_res), function(pname){
    ep <- ora_res[[pname]] %>%
        dplyr::filter(significant == TRUE) %>%
        dplyr::arrange(p_value) %>%
        dplyr::select(term_name, p_value, intersection_size) %>%
        dplyr::mutate(neglogP = -log10(p_value))

    pp <- ggplot(ep, aes(x = neglogP, y = reorder(term_name, neglogP))) +
        geom_point(aes(size = intersection_size)) +
        geom_segment(aes(xend=0, yend=term_name)) +
        theme_minimal() +
        guides(size = guide_legend(title = 'Gene Count')) +
        labs(x = expression(-log[10]*"P"), y = 'Term', title = unique(pname))

    return(pp)
}, .progress = TRUE)

names(enrichment_plots) <- names(ora_res)
```

```{r}
ggsave(glue("{figures_directory}/enrichment_analysis_on_HormoneReceptors_gene_modules.png"), enrichment_plots[['Hormone Receptors']], width = 5000, height = 10200, dpi = 800, units = 'px')

ggsave(glue("{figures_directory}/enrichment_analysis_on_Liver_gene_modules.png"), enrichment_plots[['Liver']], width = 10000, height = 15000, dpi = 800, units = 'px')

ggsave(glue("{figures_directory}/enrichment_analysis_on_BoneMarrow_gene_modules.png"), enrichment_plots[['Bone Marrow']], width = 5000, height = 10500, dpi = 800, units = 'px')
```


```{r}
ggsave(glue("{figures_directory}/enrichment_analysis_on_Program2_gene_modules.png"), enrichment_plots[['Program 2']], width = 10000, height = 10500, dpi = 800, units = 'px')

ggsave(glue("{figures_directory}/enrichment_analysis_on_Blood_gene_modules.pdf"), enrichment_plots[['Blood']], width = 10000, height = 10500, dpi = 800, units = 'px')
```

```{r}
plot(1:5)
```


```{r}
top15.hrp <- dt.consensusF %>% 
    dplyr::filter(program == 'Program 2') %>% 
    dplyr::arrange(desc(contribution)) %>% head(n = 10) %>% 
    tidyr::separate_wider_delim(col = tf_tissue, names = c('tf', 'tissue'), delim = '_', cols_remove = F) %>% 
    pull(tf) %>% 
    unique()
```

```{r}
weights.prca.hrp <- tf_grn.prca.mat[rownames(tf_grn.prca.mat) %in% top15.hrp, ] 
sum_weights.prca.hrp <- weights.prca.hrp %>% colSums() 

weights.prostate.hrp <- tf_grn.prostate.mat[rownames(tf_grn.prostate.mat) %in% top15.hrp, ] 
sum_weights.prostate.hrp <- weights.prostate.hrp %>% colSums()
```

```{r}
hrp_regulated_genes <- apply(weights.prca.hrp, 1, function(orr){
    ev <- orr[orr >= median(median_edge_weights.prca)]
    return(ev)
}, simplify = FALSE)
```

```{r}
#lapply(hrp_regulated_genes, length)
genes.hrp <- purrr::reduce(lapply(hrp_regulated_genes, names), intersect)
```

```{r}
length(hrp_regulated_genes); hrp_regulated_genes[1:5]
```

```{r}
res <- gprofiler2::gost(query = genes.hrp,
                organism = "hsapiens", ordered_query = FALSE,
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = F, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "fdr", custom_bg = colnames(weights.prca.hrp), 
                numeric_ns = "", sources = c("GO:MF", "REAC", 'GO:BP', 'GO:CC', 'KEGG'), as_short_link = FALSE, highlight = TRUE)
res$result[1:5, 1:5]; dim(res$result)
```

```{r}
hrp_regulated_genes[hrp_regulated_genes == TRUE] |> length()
```

```{r}
# select those greater than the median of medians
sum_weights.prca.hrp[sum_weights.prca.hrp > median(median_edge_weights.prca)] |> length()
sum_weights.prostate.hrp[sum_weights.prostate.hrp > median(median_edge_weights.prostate)] |> length()
```

```{r}

```

```{r}
ppn <- as.data.table(num_genes_regulated, keep.rownames = 'TF') %>%
    ggplot(aes(x = reorder(TF, -num_genes_regulated), y = num_genes_regulated)) +
    theme_minimal() +
    geom_col(col = 'black', width = 0.8, position = 'dodge') +
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    labs(x = 'Transcription Factor', y = 'Number of Genes Regulated')

ggsave(file.path(figures_directory, "number_of_genes_regulated.GRN.pdf"), ppn, width = 80000, height = 6000, dpi = 800, units = 'px', limitsize = F)
```

```{r}
dt_grn.hrp <- dt_grn %>% dplyr::filter(tf %in% tfxcan_tfs)
```

```{r}
dt_grn.hrp[1:5, 1:5]; dim(dt_grn.hrp)
```

```{r}
top15.hrp %>% noquote()
```

```{r}
gene_mat.top15.hrp <- dt_grn.hrp %>% dplyr::filter(tf %in% top15.hrp) %>% tibble::column_to_rownames('tf') %>% as.matrix() %>% t()
dim(gene_mat.top15.hrp)
```

So, across all TFs, these have the strongest weights.
```{r}
# ranking by sum of effect
rank_by_sum <- rowSums(gene_mat.top15.hrp) %>% sort(decreasing = TRUE)

# most high
gene_names.rank_by_sum <- names(rank_by_sum)
```

