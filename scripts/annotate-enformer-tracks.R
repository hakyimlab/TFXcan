

require(glue)
require(rjson)
require(readr)
require(stringr)
require(data.table)
require(GenomicRanges)
require(tidyverse)

# BiocManager::install("ENCODExplorerData")
# BiocManager::install("AnnotationHub")
# BiocManager::install("GEOquery")

library(GEOquery)
library(AnnotationHub)
library(ENCODExplorerData)


# these are ENFORMER'S input data
enformer_features_metadata <- data.table::fread('https://raw.githubusercontent.com/calico/basenji/master/manuscripts/cross2020/targets_human.txt')

# And the  accession numbers for the files are in the `identifier` column
sum(startsWith(enformer_features_metadata$identifier, 'ENCFF'))

# I need to filter for those features used in impact
# impact_features_metadata[impact_features_metadata$mark == 'DNase', ]
# list_impact_features <- split(impact_features_metadata, f = impact_features_metadata$mark)
# list_impact_features$DNase

enformer_features_metadata <- enformer_features_metadata %>%
    tidyr::separate(col=description, into = c('mark', 'detail'), sep = ':', extra = 'merge')

enformer_features_metadata$mark |> table()


# You can use AnnotationHub to retrieve encode_df_lite.
# http://bioconductor.org/packages/release/data/annotation/manuals/ENCODExplorerData/man/ENCODExplorerData.pdf
hub <- AnnotationHub()

query_encode <- query(hub, 'ENCODE', 'metadata')
query_encode$title[grepl('ENCODE', query_encode$title)]

# subsetting for ENCODE data only ==========
# these are the ENCODE builds avaialable. 
# encode_match <- !(str_match(hub, pattern = regex('.*(ENCODE).*'))[, 1] |> is.na())# need to find all ENCODE ANNOTATIONS
# hub[encode_match, ]

encode_meta_df <- subset(hub, title=='ENCODE File Metadata (Light, 2019-10-13 build)')
encode_meta_df <- encode_meta_df[[1]] |> as.data.frame()

# filter for those in the enformer 
# match with `ENCFF...` for others, match with `GSM` and `C...`

info_cols <- c('file_accession', 'target', 'assay')

# first lets filter out for the DNase data in both Enformer and Impact
enformer_dnase <- base::subset(enformer_features_metadata, mark=='DNASE')
#impact_dnase <- list_impact_features$DNase
enformer_dnase_encode <- encode_meta_df %>% dplyr::filter(file_accession %in% enformer_dnase$identifier) %>% select(all_of(info_cols))
enformer_dnase_encode$mark <- 'DNASE'
names(enformer_dnase_encode) <- c('identifier', 'target', 'assay', 'mark')

# CHIP DATA ENCODE ===
enformer_chip <- base::subset(enformer_features_metadata, mark=='CHIP')
enformer_chip_encode <- encode_meta_df %>% dplyr::filter(file_accession %in% enformer_chip$identifier) %>% select(all_of(info_cols))
enformer_chip_encode$mark <- 'CHIP'
names(enformer_chip_encode) <- c('identifier', 'target', 'assay', 'mark')

# CHIP DATA GEO ===
enformer_chip_geo <- base::subset(enformer_features_metadata, mark=='CHIP' & startsWith(identifier, 'GSM'))

target_type <- lapply(enformer_chip_geo$identifier, function(each_id){
    temp <- getGEO(each_id) |> Meta()
    temp$title
})

enformer_chip_geo <- data.frame(identifier = enformer_chip_geo$identifier, target=unlist(target_type), mark='CHIP')

a <- enformer_chip_geo[startsWith(enformer_chip_geo$target, 'batch'), ] %>% tidyr::separate(target, into=c(NA, NA, NA, 'target'), sep='_')
b <- enformer_chip_geo[startsWith(enformer_chip_geo$target, 'ChIP_'), ] %>% tidyr::separate(target, into=c(NA, 'target'), sep='_') %>% tidyr::separate(target, into=c('target', NA), sep=' ') %>% tidyr::separate(target, into=c('target', NA), sep='-')
#c <- enformer_chip_geo[startsWith(enformer_chip_geo$target, 'ChIP-seq'), ] %>% tidyr::separate(target, into=c(NA, 'target'), sep=' ')
c <- enformer_chip_geo[endsWith(enformer_chip_geo$target, 'DMI'), ] %>% tidyr::separate(target, into=c('target', NA), sep='_') %>% tidyr::separate(target, into=c(NA, 'target'), sep=' ')
d <- enformer_chip_geo[endsWith(enformer_chip_geo$target, 'uL'), ] %>% tidyr::separate(target, into=c(NA, 'target', NA), sep='_')
e <- enformer_chip_geo[!(endsWith(enformer_chip_geo$target, 'DMI') | endsWith(enformer_chip_geo$target, 'uL') | startsWith(enformer_chip_geo$target, 'ChIP_') | startsWith(enformer_chip_geo$target, 'batch')) | startsWith(enformer_chip_geo$target, 'anti-'), ] %>% tidyr::separate(col=target, into=c('target', NA, NA), sep='_')
f <- rbind(a, b, c, d, e) |> as.data.frame()

enformer_chip_geo <- dplyr::left_join(x=enformer_chip_geo, y=f, by=c('identifier' = 'identifier'))
enformer_chip_geo <- enformer_chip_geo[, -c(2, 5)]
names(enformer_chip_geo) <- c('identifier', 'mark', 'target')

# either this `enformer_chip_encode$target[enformer_chip_encode$assay == 'Histone ChIP-seq'] |> unique()`
epigenomic_marks <- c('H3K4me3', 'H3K27ac', 'H3K27me3', 'H3K4me1', 'H3K4me2', 'H3K9ac', 'H4K20me1', 'H3K36me3', 'H3K9me1', 'H3K79me2', 'H3K9me3', 'H2AFZ', 'H4K8ac', 'H3K9me2', 'H4K5ac', 'H3K56ac', 'H2AK9ac', 'H3K4ac', 'H2BK20ac', 'H3K14ac', 'H3K23me2', 'H3K18ac', 'H2BK15ac', 'H3F3A', 'H2BK5ac', 'H4K91ac', 'H2BK120ac', 'H3K79me1', 'H2BK12ac', 'H3K23ac', 'H2AK5ac', 'H3T11ph', 'H4K12ac')
enformer_chip_geo$assay <- ifelse(tolower(enformer_chip_geo$target) %in% tolower(epigenomic_marks), 'Histone ChIP-seq', 'TF ChIP-seq')


# ATAC SEQ = Not encode; need geo ===
query_encode <- query(hub, 'ncbi', 'metadata')
query_encode$title[grepl('NCBI', query_encode$title)]


enformer_atac_geo <- base::subset(enformer_features_metadata, mark=='ATAC' & startsWith(identifier, 'GSM'))
target_type <- lapply(enformer_atac_geo$identifier, function(each_id){
    temp <- getGEO(each_id) |> Meta()
    temp$title
})

enformer_atac_geo <- data.frame(identifier = enformer_atac_geo$identifier, target=unlist(target_type), mark='ATAC', assay='ATAC-seq')

# collate everything

tracks_info_list <- list(enformer_atac_geo, enformer_chip_encode, enformer_chip_geo, enformer_dnase_encode)

tracks_info_dt <- base::Reduce(function(a, b){
    rbind(a, b[, match(colnames(a), colnames(b))])
}, tracks_info_list)

info_cols <- c('index', 'identifier', 'sum_stat', 'mark' = 'mark.x', 'detail', 'target', 'assay')
tracks_info_dt_full <- dplyr::full_join(
    enformer_features_metadata[, c('index', 'identifier', 'sum_stat', 'mark', 'detail')],
    tracks_info_dt, 
    by=c('identifier' = 'identifier')
) %>% dplyr::select(all_of(info_cols))

tracks_info_dt_full$assay[tracks_info_dt_full$mark=='CAGE'] <- 'CAGE experiment'
tracks_info_dt_full$feature_names <- paste0('f', '_', tracks_info_dt_full$index)

write.table(tracks_info_dt_full, '/projects/covid-ct/imlab/users/temi/projects/TFXcan/info-files/enformer_tracks_annotated-resaved.txt', quote=F, row.names=F, sep='\t')