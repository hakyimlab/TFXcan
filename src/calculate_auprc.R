# Author: Temi
# Date: Friday August 15 2025
# Description: compute AUPRC (area under the precision-recall curve) of Enpact model predictions vs. ground
#   truth binding labels, for every TF-tissue model listed in models_list, on both the train and test splits.
#   Optionally saves a precision-recall plot per TF-tissue pair.
# Usage: Rscript calculate_auprc.R [options]

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--models_list", help='A transcription factor e.g. AR, or a comma-separated list of TFs e.g. AR,AR,FOXA1 '),
    make_option("--path_pattern", help='pattern to match path e.g. /beagle3/haky/users/temi/projects/TFPred-snakemake/data/ENPACT_734_2025-04-24/evaluation/{1}_{2}_2025-04-24.logistic.{3}_eval.txt.gz'),
    make_option("--output_file", help = "The name of output file. Should be appended with .gz since this will be a compressed file."),
    make_option("--plots_directory", help = "A folder for saving AUPRC plots. If supplied, will save a plot for each TF-tissue pair.", default = NULL)
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



theme_simple <- theme_bw() + theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.title = element_text(family = 'Helvetica', colour = 'black', size = 9.33),
    text=element_text(size=9.33, colour = 'black', family = 'Helvetica'),
    axis.text = element_text(hjust = 1), 
    axis.text.x = element_text( hjust=0.5), 
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
)


plot_pr_curve.ggplot <- function(pr_object, tf_tissue){
    dt_pr_curve <- pr_object$curve %>% as.data.table() %>% setnames(c('recall', 'precision', 'unk'))

    tlabel <- glue("AUPRC: {round(pr_object$auc.integral, 3)}")
    ttitle <- glue("Precision-Recall curve for\nEnpact's {tf_tissue} Model")

    pp <- dt_pr_curve %>% ggplot(aes(x = recall, y = precision)) +
        ggplot2::geom_line() +
        ylim(0, 1) +
        xlim(0, 1) +
        labs(x = 'Recall', y = 'Precision', title = ttitle) +
        theme_simple +
        geom_text(label = tlabel, col = 'darkgreen', x = 0.2, y = 0.2, size = 2.0)

    return(pp)

    #
}


calculate_auprc <- function(dt, transcription_factor, tissue, type, plots_dir = NULL){
    tryCatch({
        cl0 <- base::subset(dt, binding_class == 0)$TFPred_score # background
        cl1 <- base::subset(dt, binding_class == 1)$TFPred_score # foreground
        cl0 <- cl0[!is.na(cl0)]
        cl1 <- cl1[!is.na(cl1)]
        pr <- PRROC::pr.curve(scores.class0 = cl1, scores.class1=cl0, curve = TRUE)
        if(!is.null(plots_dir)){
            prplot <- plot_pr_curve.ggplot(pr, glue("{transcription_factor}-{tissue}"))
            ggsave(glue("{plots_dir}/AUPRC.{type}.{transcription_factor}-{tissue}.png"), prplot, width = 2000, height = 2000, dpi = 800, units = 'px')
        }
        return(cbind(transcription_factor, tissue, type, pr$auc.integral))
    }, error = function(e){
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
            evaluations <- data.frame(assay = qi[1], context = qi[2], type = qi[3], prauc = NA)
        } else if (any(nrow(dt) > 0)) {
            evaluations <- calculate_auprc(dt, qi[1], qi[2], qi[3], plots_dir = opt$plots_directory)
            evaluations <- evaluations %>% as.data.frame() %>%
                stats::setNames(c('assay', 'context', 'type', 'prauc')) 
        }
    }
    return(evaluations)
}

data.table::fwrite(evaluations, file=opt$output_file, sep='\t', row.names=F, col.names=T, quote=F, compress = 'gzip')

doParallel::stopImplicitCluster()