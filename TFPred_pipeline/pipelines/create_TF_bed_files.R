

#Rscript -e 'install.packages(c("yaml", "data.table", "glue", "tidyverse"), repos="https://cloud.r-project.org")'

# TF <- 'FOXA1'
# tissue_type <- 'Breast'

args <- commandArgs(trailingOnly = TRUE)

TF <- args[1]
tissue_type <- args[2]
output_dir <- args[3]

print(TF)
print(tissue_type)
print(output_dir)

library(yaml)
library(data.table)
library(glue)
library(tidyverse)

setwd('/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/pipelines')
directives <- yaml::yaml.load_file("../config_files/pipeline.yaml")

tf_db <- data.table::fread(directives$TF_table)
# tf_db %>% dplyr::filter(Factor == TF) %>% dplyr::group_by(Cell_line) %>% summarise(n_files=n()) %>% arrange(desc(n_files))

#base::subset(x=tf_db, subset = Factor == TF)

# get the ids 
TF_data <- base::subset(x=tf_db, subset = (Factor == TF) & (Tissue_type == tissue_type))
# 
if( ! (nrow(TF_data) > 1 & ncol(TF_data) > 1)){
    stop('')
} else {
    TF_files <- list.files(glue('{directives$data_dir}/human_factor'), pattern=paste0('^', TF_data$DCid, collapse='_*|'), full.names=T)
}

sapply(TF_files, function(each_file){
    R.utils::createLink(link=each_file, target=glue("{output_dir}/{basename(each_file)}"), overwrite=TRUE)
})

# find homer motifs
# tf_lower <- tolower(TF)
# potential_motif_files <- list.files(glue('{directives$homer$dir}/data/knownTFs/motifs'), glue('^{tf_lower}'), full.names=T)[1]


#glue("{directives$project_dir}/scripts/utilities/scan_for_motifs.sh")


