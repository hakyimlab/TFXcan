
# Matrix Factorization on enpact weights

```{r}
enpact_weights <- data.table::fread("/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2024-07-26/statistics/ENPACT_734_2024-07-26.compiled_weights.lambda.1se.txt.gz")
```

```{r}
# remove some weird tissues
weird_tissues <- c('HCT116', "HeLacontaminant", "Headandneck", "LNCaPcells")
remove_these <- rownames(enpact_weights)[grepl(paste0("_", weird_tissues, collapse="|"), rownames(enpact_weights))]
enpact_weights <- enpact_weights[!rownames(enpact_weights) %in% remove_these, ]
# enpact_weights <- tibble::rownames_to_column(enpact_weights, 'feature')

data.table::fwrite(enpact_weights, file=glue('{project_dir}/files/ENPACT_734_2024-07-26.compiled_weights.lambda.1se.txt.gz'), col.names=TRUE, row.names=FALSE, quote=F, compress='gzip', sep='\t')
enpact_weights[1:5, 1:5]
```

```{r}
datafile <- glue('{project_dir}/files/ENPACT_734_2024-07-26.compiled_weights.lambda.1se.txt.gz')
column_for_rownames <- 'feature'
output_basename <- glue('{project_dir}/files/runFlashier')
priorL <- "ebnm_point_exponential"
priorF <- "ebnm_point_exponential"
transpose <- "TRUE"
greedy_Kmax <- 500L
```

```{r}
# different combinations of priorL and priorF
possible_priors <- c("ebnm_point_exponential", "ebnm_point_normal", "ebnm_point_laplace", "ebnm_normal", "ebnm_normal_scale_mixture")
combination_priors <- expand.grid(possible_priors, possible_priors)
```

```{r}
for(i in 1:nrow(combination_priors)){
    priorL <- combination_priors[i, 1]
    priorF <- combination_priors[i, 2]
    cmd <- glue('sbatch {project_dir}/src/runFlashier.sbatch {datafile} {column_for_rownames} {output_basename} {priorL} {priorF} {transpose} {greedy_Kmax}')
    print(cmd);
    system(cmd)
}
```

```{r}
# if needed, you should edit the sbatch file: {project_dir}/src/predictCWASscores.sbatch
cmd <- glue('sbatch {project_dir}/src/runFlashier.sbatch {datafile} {column_for_rownames} {output_basename} {priorL} {priorF} {transpose} {greedy_Kmax}')
cmd
```

```{r}
system(cmd) ; system('squeue -u temi')
```

```{r}
ff <- readRDS(file=glue('{project_dir}/files/testrunFlashier.ebnm_point_exponential-ebnm_point_exponential.rds.gz'))
```

```{r}
names_out <- c() #glue("{priorL}-{priorF}")
out <- apply(combination_priors, 1, function(each_priors){
    priorL <- each_priors[1]
    priorF <- each_priors[2]
    names_out <<- c(names_out, glue("{priorL}-{priorF}"))
    ff <- glue('{project_dir}/files/runFlashier.{priorL}-{priorF}.rds.gz')
    if(file.exists(ff)){
        flash_object <- readRDS(file=ff)
        return(flash_object)
    } else {
        return(NULL)
    }
    
})
names(out) <- names_out
out <- Filter(Negate(is.null), out)
```

```{r}
out$`ebnm_point_exponential-ebnm_point_exponential`$pve
```

```{r}
lapply(out, function(each_out){
    sum(each_out$pve)
})
```

```{r}
lapply(out, function(each_out){
    each_out$n_factors
})
```

```{r}
lapply(out, function(each_out){
    each_out$factor_weights
})
```

```{r}
lapply(out, function(each_out){
    each_out$factor_weights
})
```

```{r}
xt <- readRDS('/project2/haky/temi/projects/enpact-predict-snakemake/output/1KG_2024-10-04/1KG_2024-10-04.enpact_scores.array.rds.gz')
```
