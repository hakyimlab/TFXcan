
library(data.table)
library(dplyr)
library(tidyr)
library(glue)

imlab_dir <- '/lus/grand/projects/covid-ct/imlab'
project_dir <- glue('{imlab_dir}/users/temi/projects/TFXcan/TFpred_pipeline')

enformer_db <- data.table::fread(glue('{project_dir}/../modeling_pipeline/metadata/enformer_tracks_annotated-resaved.txt'))
enformer_db <- enformer_db[enformer_db$assay == 'TF ChIP-seq', ] # using only TF chip-seq data

cistrome_db <- data.table::fread(glue('{imlab_dir}/data/cistrome/raw/human_factor_full_QC.txt'))

homer_db <- sapply(strsplit(toupper(list.files(glue('/lus/grand/projects/covid-ct/imlab/users/temi/software/homer/motifs'))), '\\.'), getElement, 1) |> unique() |> data.frame()
colnames(homer_db) <- 'homer_TF'

# curate enformer
ue <- unique(enformer_db$target) |> sort()
enformer_present <- ue[ue %in% homer_db$homer_TF]
enformer_present <- cbind(enformer_present, 1)
enformer_absent <- ue[!ue %in% homer_db$homer_TF]
enformer_absent <- cbind(enformer_absent, 0)

enformer_TF_db <- as.data.frame(rbind(enformer_present, enformer_absent))
colnames(enformer_TF_db) <- c('TF', 'status')

# curate cistrome
uc <- unique(cistrome_db$Factor)
cistrome_present <- uc[uc %in% homer_db$homer_TF]
cistrome_present <- cbind(cistrome_present, 1)
cistrome_absent <- uc[!uc %in% homer_db$homer_TF]
cistrome_absent <- cbind(cistrome_absent, 0)

cistrome_TF_db <- as.data.frame(rbind(cistrome_present, cistrome_absent))
colnames(cistrome_TF_db) <- c('TF', 'status')

# curate homer_db
homer_TF_db <- cbind(homer_db, 1) |> as.data.frame()
colnames(homer_TF_db) <- c('TF', 'status')

# curate TF database
#queries <- list(enformer_TF_db, homer_TF_db, cistrome_TF_db)

#TF_db <- queries %>% purrr::reduce(full_join, by = c('TF')) 
#colnames(TF_db) <- c('TF', 'enformer_db', 'homer_db', 'cistrome_db')

a <- base::merge(homer_TF_db, enformer_TF_db, by='TF', all=T, suffixes=c('_homer_db', '_enformer_db'))
TF_db <- merge(a, cistrome_TF_db, by='TF', all=T)
colnames(TF_db) <- c('TF', 'homer_db', 'enformer_db', 'cistrome_db')

TF_db <- TF_db %>% 
    as.data.frame() %>%
    dplyr::mutate(across(c('homer_db', 'enformer_db', 'cistrome_db'), as.numeric))
TF_db[is.na(TF_db)] <- 0

cistrome_cnts <- cistrome_db[cistrome_db$Factor %in% TF_db$TF[TF_db$cistrome_db == 1], ] %>%
    dplyr::group_by(Factor) %>%
    dplyr::summarize(cistrome_files = n()) %>%
    dplyr::arrange(desc(cistrome_files))

TF_db <- merge(TF_db, cistrome_cnts, by.x='TF', by.y='Factor', all=T)
TF_db[is.na(TF_db)] <- 0

data.table::fwrite(TF_db, glue('{project_dir}/metadata/transcription_factor_avaialability.txt'), quote=F, sep='\t')
# ===== 
TF_db <- data.table::fread(glue('{project_dir}/metadata/transcription_factor_avaialability.txt'), fill=T)
base::subset(TF_db, (homer_db == 1 & cistrome_db == 1 & enformer_db == 0)) %>% dplyr::arrange(desc(cistrome_files)) 

base::subset(TF_db, (homer_db == 1 & cistrome_db == 1 & enformer_db == 1))


