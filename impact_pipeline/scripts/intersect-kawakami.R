
setwd('/projects/covid-ct/imlab/users/temi/projects/TFXcan/scripts/')

rm(list=ls())

library(glue)
library(GenomicRanges)
library(reticulate)
library(R.utils)
library(data.table)
library(tidyverse)

# pick a transcription factor
TF <- 'FOXA1' #
# cell_line <- 'LuCaP' # a prostrate cancer cell line
# foxa1_motif <- 'foxa1.lncap.motif'

# where is the TF chip-seq data, as well as other directories
work_dir <- '/projects/covid-ct/imlab/users/temi/projects/TFXcan'
kawakami_dir <- glue('{work_dir}/data/kawakami-human') # there is a 'detailed_info/' folder that holds the files

# all the files seem to be sorted so I don't need to worry about sorting them

# # sorted directory
# sort_dir <- glue('{kawakami_dir}/sorted_bed_files')
# if(!dir.exists(sort_dir)){
#     dir.create(sort_dir, recursive=T)
# } else {
#     print('Dir exists.')
# }

peak_files_dir <- glue('{kawakami_dir}/detailed_info')
homer_dir <- '~/miniconda3/envs/r-env/share/homer'
project_dir <- '../'
#output_dir <- glue('../processed-data/{TF}/{cell_line}/freedman')

peak_files_full_path <- list.files(glue('{peak_files_dir}'), full.names=T)
FOXA1_bed_files <- peak_files_full_path[startsWith(list.files(glue('{peak_files_dir}')), TF)] # only the TF files
FOXA1_files <- lapply(FOXA1_bed_files, data.table::fread) # read in the files into like a list

#FOXA1_bed_files <- peak_files_full_path[startsWith(list.files(glue('{peak_files_dir}')), TF)] # only the TF files
file_names <- sapply(strsplit(x = FOXA1_bed_files, split='/', fixed=T), function(ea){
    ea[12]
})

# create a folder to hold the 'sorted' bed files
# sorted directory
sort_dir <- glue('{kawakami_dir}/sorted_bed_files')
if(!dir.exists(sort_dir)){
    dir.create(sort_dir, recursive=T)
} else {
    print('Dir exists.')
}

# save the files without the headers
sapply(seq_along(FOXA1_bed_files), function(i){
    write.table(x=FOXA1_files[[i]], file=glue('{sort_dir}/{file_names[i]}'), col.names=F, quote=F, row.names=F, sep='\t')
})

sorted_dir <- glue('{kawakami_dir}/sorted_bed_files')
peak_files_full_path <- list.files(glue('{sorted_dir}'), full.names=T)
FOXA1_bed_files <- peak_files_full_path[startsWith(list.files(glue('{sorted_dir}')), TF)] # only the TF files

# dump these into a file
sink(file=glue('{kawakami_dir}/chip-file-{TF}.txt'))
cat(FOXA1_bed_files)
sink()

# # need to ensure the files are sorted
# for(i in 1:length(FOXA1_bed_files)){
#     cmd <- glue('sort -k1,1V -k2,2n {FOXA1_bed_files[i]} > {sort_dir}/{file_names[i]}')
#     system(cmd)
# }

# find intersections
bedtools_cmd <- '/home/temi/miniconda3/envs/compbio-tools/bin/bedtools'
cmd <- glue('{bedtools_cmd} multiinter -header -i $(cat {kawakami_dir}/chip-file-{TF}.txt) > {kawakami_dir}/{TF}_intersected_files.bed')
system(cmd)


# start from here
bedfiles_names <- sapply(strsplit(file_names, split='_'), function(each){
    paste0(each[1:2], collapse='_')
})

FOXA1_intersected <- data.table::fread(glue('{kawakami_dir}/FOXA1_intersected_files.bed'), col.names=c('chrom', 'start', 'end', 'num', 'list', bedfiles_names))

chip_counts <- sapply(FOXA1_intersected$list, function(f){
    getElement(strsplit(f, split=','), 1) |> length()
})

FOXA1_intersected$bedcounts <- chip_counts
write.table(FOXA1_intersected, file=glue('{kawakami_dir}/FOXA1_intersected_chip.bed'), sep='\t', col.names=T, row.names=F)



# 

FOXA1_intersected <- data.table::fread(glue('{kawakami_dir}/FOXA1_intersected_files.bed'))
dim(FOXA1_intersected)

bedfiles_names <- colnames(FOXA1_intersected)[!colnames(FOXA1_intersected) %in% c("chrom","start","end","num","list",'bedcounts')]

regionCounts <- sapply(1:length(bedfiles_names), function(f){
    which(chip_counts == f) |> length()
})

bp <- barplot(regionCounts, names.arg=1:length(bedfiles_names), ylim=c(0, max(regionCounts + 100000)), ylab='Number of overlapping ChIP peaks', xlab='Number of files (individuals)')
#bp <- axis(1, at=1:length(bedfiles_names), labels=F, tick=T)
text(x=bp, y=regionCounts, label=regionCounts, cex=0.8, pos=4, srt=45) 
mtext('Number of ChIP-seq peaks overlapped across number of individuals')

regionPeaks <- sapply(1:length(bedfiles_names), function(f){
    which(chip_counts == f)
})

overlap_29 <- FOXA1_intersected[regionPeaks[29] |> unlist(), ]

overlap_29$chrom |> table() %>% sort(decreasing=T)

overlap_dir <- glue('{work_dir}/processed-data/FOXA1/overlapped_29')
write.table(overlap_29[, 1:3], file=glue('{overlap_dir}/input_chip.bed'), sep='\t', col.names=F, row.names=F)


# read the data
overlap_29 <- data.table::fread(glue('{overlap_dir}/input_chip.bed'))

