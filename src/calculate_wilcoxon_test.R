# Author: Temi
# Date: Friday August 15 2025
# Description: run a Wilcoxon rank-sum test comparing Enpact's predicted binding scores (TFPred_score) between
#   bound (binding_class == 1) and unbound (binding_class == 0) loci, for every TF-tissue model listed in
#   models_list, on both the train and test splits.
# Usage: Rscript calculate_wilcoxon_test.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--models_list", help='A transcription factor e.g. AR, or a comma-separated list of TFs e.g. AR,AR,FOXA1 '),
    make_option("--path_pattern", help='pattern to match path e.g. /beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2025-04-24/evaluation/{1}_{2}_2025-04-24.logistic.{3}_eval.txt.gz'),
    make_option("--output_file", help = "The name of output file. Should be appended with .gz since this will be a compressed file.")
)

opt <- parse_args(OptionParser(option_list=option_list))

library(data.table)
library(PRROC)
library(stringr)
library(foreach)
library(glmnet)
library(glue)
library(ggplot2)


# opt <- list()
# opt$transcription_factors <- 'FOXA1,MYC,YY1,CTCF,GATA4,TAL1,SP1,MYC,MYC,VDR,CTCF,MYC,NEUROG2,CTCF,GATA2,PPARG,FOSL2,CTCF,CTCF,CTCF,CTCF,NANOG,FOXA2,ERG,REST,FOXA1,POU5F1,GATA1,FOXM1,GATA3,GATA3,CTCF,FOXA1,CTCF,CTCF,SPI1,RUNX1,CTCF,HOXB13,AR,CTCF,GATA2,E2F1,ESR1,NR3C1,ESR1,RELA,NR3C1'
# opt$tissues <- 'MammaryGland,Breast,Blood,EmbryonicKidney,Embryo,BoneMarrow,Colon,Blood,Colon,Blood,MammaryGland,BoneMarrow,FetalLung,Liver,BoneMarrow,Adipose,Lung,Brain,Skin,UmbilicalVein,Embryo,Embryo,Skin,Prostate,Blood,Breast,Embryo,BoneMarrow,Breast,MammaryGland,Breast,Lung,Prostate,Colon,BoneMarrow,Blood,Blood,Breast,Prostate,Prostate,Blood,Prostate,Prostate,MammaryGland,Lung,Breast,Blood,Breast'
# opt$path_pattern <- '/beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2025-04-24/evaluation/{1}_{2}_2025-04-24.logistic.{3}_eval.txt.gz'

dt_models <- data.table::fread(opt$models_list)
tflist <- dt_models$transcription_factor
tissuelist <- dt_models$tissue

pp <- c("\\{1\\}", "\\{2\\}", "\\{3\\}")

qMatrix <- rbind(
    cbind(tflist, tissuelist, 'train'),
    cbind(tflist, tissuelist, 'test')
) |> as.data.frame() |> setNames(c('transcription_factor', 'tissue', 'type')) %>%
    dplyr::arrange(transcription_factor, tissue, type)


calculate_wilcoxon_pvalues <- function(dt, transcription_factor, tissue, type){

    tryCatch({
        dt <- dt[complete.cases(dt), ]

        # wilcoxon
        wiltest <- wilcox.test(TFPred_score ~ binding_class, data=dt)
        wiltest_pvalue <- wiltest$p.value
        
        return(cbind(transcription_factor, tissue, type, wiltest_pvalue))
    }, error = function(e){
        print(e)
        return(cbind(transcription_factor, tissue, type, NA))
    })
   
}

num_clusters <- 4 #- 5 # 12 - 5
doParallel::registerDoParallel(num_clusters)

evaluations <- foreach::foreach(i=seq_len(nrow(qMatrix)), .combine='rbind', .inorder=T) %dopar% {
    # do something
    qi <- qMatrix[i, ] |> unlist()
    tfile <- Sys.glob(stringr::str_replace_all(opt$path_pattern, purrr::set_names(qi, pp)))
    #mpath <- Sys.glob(stringr::str_replace_all(opt$model_path, purrr::set_names(qi[1:2], pp[1:2])))

    if(file.exists(tfile)){
        dt <- data.table::fread(tfile)
        dt <- dt[complete.cases(dt), ]
        dt <- dt[!duplicated(dt[[c("locus")]]), ]
        message(glue::glue("INFO - {qi[1]} in {qi[2]} ({qi[3]}) has {nrow(dt)} rows and {ncol(dt)} columns."))
        if(any(nrow(dt) == 0)){
            evaluations <- data.frame(assay = qi[1], context = qi[2], type = qi[3], wiltest_pvalue = NA)
        } else if (any(nrow(dt) > 0)) {
            evaluations <- calculate_wilcoxon_pvalues(dt, qi[1], qi[2], qi[3])
            evaluations <- evaluations %>% as.data.frame() %>%
                stats::setNames(c('assay', 'context', 'type', 'wiltest_pvalue')) 
        }
    }
    return(evaluations)
}

data.table::fwrite(evaluations, file=opt$output_file, sep='\t', row.names=F, col.names=T, quote=F, compress = 'gzip')

doParallel::stopImplicitCluster()