

# This part of the script analyses the intersected kawakami data

setwd('/projects/covid-ct/imlab/users/temi/projects/TFXcan/scripts/')

rm(list=ls())

library(glue)
library(GenomicRanges)
library(reticulate)
library(R.utils)
library(data.table)
library(tidyverse)

# kawakami intersections

TF <- 'FOXA1'
# where is the TF chip-seq data, as well as other directories
work_dir <- '/projects/covid-ct/imlab/users/temi/projects/TFXcan'
kawakami_dir <- glue('{work_dir}/data/kawakami-human') # there is a 'detailed_info/' folder that holds the files

tf_intersected <- data.table::fread(glue('{kawakami_dir}/{TF}_intersected_files.bed'))
dim(tf_intersected)

# get the names of the bedfiles
bedfiles_names <- colnames(tf_intersected)[!colnames(tf_intersected) %in% c("chrom","start","end","num","list",'bedcounts')]

# count the number of files that have the peaks in common
chip_counts <- sapply(tf_intersected$list, function(f){
    getElement(strsplit(f, split=','), 1) |> length()
})
# how many of them have 1 overlap, 2 overlaps, 3, 4, 5, 6, up to the length of bedfiles_names
regionCounts <- sapply(1:length(bedfiles_names), function(f){
    which(chip_counts == f) |> length()
})


# KAWAKAMI before intersection with motif


# plot the distribution of common peaks across different ovelaps of the files
bp <- barplot(regionCounts, names.arg=1:length(bedfiles_names), ylim=c(0, max(regionCounts + 100000)), ylab='Number of overlapping ChIP peaks', xlab='Number of files (individuals)')
#bp <- axis(1, at=1:length(bedfiles_names), labels=F, tick=T)
text(x=bp, y=regionCounts, label=regionCounts, cex=0.8, pos=4, srt=45) 
mtext('Number of ChIP-seq peaks overlapped across number of individuals')

# from the plot, there is one peak common to all 60 (and 0 common to 61 or 62 files)

# ===
regionPeaks <- sapply(1:length(bedfiles_names), function(f){
    which(chip_counts == f)
})

# which peaks are common to all vs to any number of files/individuals
overlap_n <- length(bedfiles_names) # i.e. 62

overlap_62 <- tf_intersected[regionPeaks[overlap_n] |> unlist(), ]
dim(overlap_62)

# I need to know more details about these experiments === or maybe not
gsm_files <- sapply(str_split(string=bedfiles_names, pattern='_'), getElement, 2)

gsm_files_collapsed <- paste(gsm_files, collapse=' ')


kawakami_files_dir <- glue('{kawakami_dir}/kawakami-files-details')
if(!dir.exists(kawakami_files_dir)){
    dir.create(kawakami_files_dir)
} else {
    print('Directory exists.')
}

ffq_tool <- '/home/temi/miniconda3/envs/compbio-tools/bin/ffq'
#system(glue('{ffq_tool} {gsm_files_collapsed} -o {kawakami_dir}/kawakami-files-details --split'))

lapply(gsm_files, function(each_gsm){
    system(glue('{ffq_tool} {each_gsm} -o {kawakami_files_dir}/{each_gsm}_details.json'))
})

# overlap_60 <- tf_intersected[regionPeaks[60] |> unlist(), ]
# dim(overlap_60)

# overlap_29$chrom |> table() %>% sort(decreasing=T)

# overlap_dir <- glue('{work_dir}/processed-data/FOXA1/overlapped_29')
# write.table(overlap_29[, 1:3], file=glue('{overlap_dir}/input_chip.bed'), sep='\t', col.names=F, row.names=F)


# # read the data
# overlap_29 <- data.table::fread(glue('{overlap_dir}/input_chip.bed'))


# common to 42
# which peaks are common to all vs to any number of files/individuals
overlap_42 <- 42
overlap_42 <- tf_intersected[regionPeaks[overlap_42] |> unlist(), ]
dim(overlap_42)

tf_chip_intersected <- tf_intersected[, c("chrom","start","end","num")]

# need to intesect these with motifs for FOXA1
tf_chip_intersected <- lapply(split(tf_chip_intersected, f=tf_chip_intersected$num), function(each_tf){
    each_tf$id <- paste('peak', each_tf$num, 1:nrow(each_tf), sep='_')
    return(each_tf)
}) 

tf_chip_intersected <- do.call(rbind, tf_chip_intersected)

kawakami_impact_dir <- glue('{kawakami_dir}/kawakami-impact')
if(!dir.exists(kawakami_impact_dir)){
    dir.create(kawakami_impact_dir)
} else {
    print('Directory exists.')
}
#overlap_dir <- glue('{kawakami_dir}/processed-data/FOXA1/overlapped_29')
write.table(tf_chip_intersected, file=glue('{kawakami_impact_dir}/{TF}_intersected_chip.bed'), sep='\t', col.names=F, row.names=F)


## KAWAKAMI data after intersecting with the Motif

# pick a transcription factor
TF <- 'FOXA1' #
cell_line <- 'LuCaP' # a prostrate cancer cell line
foxa1_motif <- 'foxa1.lncap.motif'

# where is the TF chip-seq data, as well as other directories
work_dir <- '/projects/covid-ct/imlab/users/temi/projects/TFXcan'
kawakami_dir <- glue('{work_dir}/data/kawakami-human')
kawakami_impact_dir <- glue('{kawakami_dir}/kawakami-impact')

tp <- read.table(glue("{kawakami_impact_dir}/{TF}_{cell_line}_train_test_positive_bed.txt"))
peak_info <- read.table(glue("{kawakami_impact_dir}/{TF}_{cell_line}_true_positive_peak_info.txt"))
peak_match <- match(peak_info$peak_id, tp$V4)
peak_info_details <- peak_info[peak_info$peak_id %in% tp$V4, ]
peak_distribution <- sapply(strsplit(peak_info_details$kawakami_id, split='_'), getElement, 2)
peak_distribution <- table(peak_distribution) |> as.data.frame() 
colnames(peak_distribution) <- c('count', 'freq')
peak_distribution <- peak_distribution %>% dplyr::mutate(count=as.numeric(count)) %>% dplyr::arrange(count)

bp <- barplot(peak_distribution$freq, names.arg=peak_distribution$count, 
    ylim=c(0, max(peak_distribution$freq + 10000)), 
    ylab='No. of overlapping true ChIP peaks', xlab='Number of files (individuals)')

#bp <- axis(1, at=1:length(bedfiles_names), labels=F, tick=T)
#text(x=bp, y=peak_distribution$freq, label=peak_distribution$freq, cex=0.8, pos=4, srt=90) 
mtext('Number of ChIP-seq peaks overlapped across number of Kawakami files')


layout(matrix(c(1, 2), nrow=2))
# plot the distribution of common peaks across different ovelaps of the files
bp <- barplot(regionCounts, names.arg=1:length(bedfiles_names), ylim=c(0, max(regionCounts + 100000)), ylab='Number of overlapping ChIP peaks', xlab='Number of files (individuals)')
#bp <- axis(1, at=1:length(bedfiles_names), labels=F, tick=T)
text(x=bp, y=regionCounts, label=regionCounts, cex=0.8, pos=4, srt=45) 
mtext('Number of ChIP-seq peaks overlapped across number of individuals')

bp <- barplot(peak_distribution$freq, names.arg=peak_distribution$count, 
    ylim=c(0, max(peak_distribution$freq + 10000)), 
    ylab='No. of overlapping true ChIP peaks', xlab='Number of files (individuals)')
#bp <- axis(1, at=1:length(bedfiles_names), labels=F, tick=T)
text(x=bp, y=peak_distribution$freq, label=peak_distribution$freq, cex=0.8, pos=4, srt=90) 
mtext('Number of overlapping true positives after intersecting with known motif sites')
