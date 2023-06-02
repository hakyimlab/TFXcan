

rm(list=ls())

setwd('/projects/covid-ct/imlab/users/temi/projects/TFXcan/scripts/')
source('/projects/covid-ct/imlab/users/temi/projects/TFXcan/scripts/custom-impact.R')

library(glue)
library(GenomicRanges)
library(reticulate)
library(R.utils)
library(data.table)
library(tidyverse)
library(parallel)

# pick a transcription factor
TF <- 'FOXA1' #
cell_line <- 'LuCaP' # a prostrate cancer cell line
foxa1_motif <- 'foxa1.lncap.motif'

# where is the TF chip-seq data, as well as other directories
work_dir <- '/projects/covid-ct/imlab/users/temi/projects/TFXcan'
freedman_dir <- glue('{work_dir}/data/freedman')
freedman_impact_dir <- glue('{work_dir}/data/freedman')
freedman_impact_split_dir <- glue('{freedman_dir}/freedman-impact-split')
chip_data_dir <- glue('{freedman_impact_split_dir}/chip_data_files')
common_files_dir <- glue('{work_dir}/data/common-files/{TF}')
homer_dir <- glue('~/miniconda3/share/homer')
project_dir <- glue('{work_dir}/data/common-files')
#freedman_impact_split_dir <- glue('../processed-data/{TF}/{cell_line}/freedman')

# ensure that the directories are all created

lapply(c(freedman_dir, freedman_impact_split_dir, chip_data_dir), function(each_dir){

    if(dir.exists(paths=each_dir)){
        cat(paste(each_dir, 'exists.', sep=' '))
    } else {
        dir.create(path=each_dir, recursive=T)
        cat(paste(each_dir, 'created.', sep=' '))
    }
})

# ensure all needed files are in common files
# make sure that the TF's Chip data and motif data exist
# check for the motif in Homer
if(!file.exists(glue('{homer_dir}/data/knownTFs/motifs/{tolower(foxa1_motif)}'))){
    print(glue('No motif information for {TF} available'))
} else {
    print(glue('Motif information for {TF} available. Moving on.'))
}

# different cell line but using it as a proxy
if(!file.exists(glue('{common_files_dir}/{TF}_{cell_line}.motif'))){
    # copy the motif file from homer
    cmd <- glue('cp {homer_dir}/data/knownTFs/motifs/{tolower(foxa1_motif)} {common_files_dir}/{TF}_{cell_line}.motif')
    system(cmd)
    print(glue('{common_files_dir}/{TF}_{cell_line}.motif copied!'))
} else {
    print(glue('{common_files_dir}/{TF}_{cell_line}.motif exists'))
}

# remake the motif information file
pwm <- read.table(glue('{common_files_dir}/{TF}_{cell_line}.motif'), sep = "\t", header = F, stringsAsFactors = FALSE, fill=T)
pwm[1, 4] <- 'PH'
pwm[1, 2] <- glue('{TF}_{cell_line}') # NOTE THIS PART CAREFULLY
write.table(pwm[, 1:4], glue('{common_files_dir}/{TF}_{cell_line}.motif'), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

#Input Motifs
MotifThreshvals <- pwm[which(pwm[,4]=="PH"), 3]
MotifThreshnames <- pwm[which(pwm[,4]=="PH"), 2]
MotifThreshinfo <- cbind(MotifThreshnames, MotifThreshvals)

write.table(MotifThreshinfo, glue('{common_files_dir}/{TF}_{cell_line}_motif_threshold_info.txt'), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

# top 15000 motifs for the TF are already available; there is no point re-running this 
if(!file.exists(glue('{common_files_dir}/{TF}_{cell_line}_scanMotifsGenomewide_sort_15000.txt'))){
    cmd <- glue('perl {homer_dir}/bin/scanMotifGenomeWide.pl {common_files_dir}/{TF}_{cell_line}.motif {homer_dir}/data/genomes/hg19 > {common_files_dir}/{TF}_{cell_line}_scanMotifsGenomewide.txt')
    system(cmd)
    cmd <- glue("sort -t $'\t' -k6,6rn {common_files_dir}/{TF}_{cell_line}_scanMotifsGenomewide.txt > {common_files_dir}/{TF}_{cell_line}_scanMotifsGenomewide_sort.txt")
    system(cmd)
    cmd <- glue('head -n 15000 {common_files_dir}/{TF}_{cell_line}_scanMotifsGenomewide_sort.txt > {common_files_dir}/{TF}_{cell_line}_scanMotifsGenomewide_sort_15000.txt')
    system(cmd)
} else {
    print(glue('{common_files_dir}/{TF}_{cell_line}_scanMotifsGenomewide_sort_15000.txt exists.'))
}


# read in the intersected file
chip_peak_headers <- c('chrom', 'start', 'end', 'num')
chip_peak_file <- data.table::fread(glue('{freedman_impact_dir}/{TF}_intersected_chip.bed'), select=chip_peak_headers)
chip_peak_file[, 'strand'] <- '*'

# check if the intersection files already have an id column if not, create one
if(!'id' %in% colnames(chip_peak_file)){
    # need to intesect these with motifs for FOXA1
    tf_chip_intersected <- lapply(split(chip_peak_file, f=chip_peak_file$num), function(each_tf){
        each_tf$id <- paste('peak', each_tf$num, 1:nrow(each_tf), sep='_')
        return(each_tf)
    }) 

    chip_peak_file <- do.call(rbind, tf_chip_intersected)
}

peak_info <- cbind(chip_peak_file$id, paste('peak', 1:length(chip_peak_file$id), sep='')) |> as.data.frame()
colnames(peak_info) <- c('data_id', 'peak_id')
write.table(peak_info, glue("{freedman_impact_split_dir}/{TF}_{cell_line}_true_positive_peak_info.txt"), sep = "\t", quote = F, row.names=F)

# split by num
chip_peak_file_split <- split(chip_peak_file, f=chip_peak_file$num)

# using mclapply
# numCores <- detectCores() - 10
#run_output <- mclapply(X=seq_along(chip_peak_file_split), FUN=annotate_with_impact, output_dir=freedman_impact_split_dir, common_dir=common_files_dir, mc.cores=numCores)
run_output <- lapply(X=seq_along(chip_peak_file_split), FUN=annotate_with_impact, output_dir=freedman_impact_split_dir, common_dir=common_files_dir)


# =======================
# read in the peak files
peak_info <- data.table::fread(glue("{freedman_impact_split_dir}/{TF}_{cell_line}_true_positive_peak_info.txt"))

allfiles <- list.files(glue("{freedman_impact_split_dir}"))
positive_files <- allfiles[endsWith(x=allfiles, suffix='FOXA1_LuCaP_train_test_positive_bed.txt')]
#peak_info <- allfiles[endsWith(x=allfiles, suffix='FOXA1_LuCaP_true_positive_peak_info.txt')]

# create enformer beds =====
tp_list <- lapply(positive_files, function(filename){
    read.table(glue("{freedman_impact_split_dir}/{filename}"))
})

names(tp_list) <- sapply(strsplit(positive_files, split='_'), function(each){paste(each[1], each[2], sep='_')})

# peak_info_list <- lapply(peak_info, function(filename){
#     read.table(glue("{kawakami_impact_split_dir}/{filename}"))
# })
# names(peak_info_list) <- sapply(strsplit(peak_info, split='_'), function(each){paste(each[1], each[2], sep='_')})
#tn <- read.table(glue("{kawakami_impact_split_dir}/{TF}_{cell_line}_train_test_negative_bed.txt"))

# compare with the peak information file

#peak_match <- match(peak_info$V2, tp$V4)

peak_distribution <- lapply(tp_list, function(each_tp){

    peak_info_details <- peak_info[peak_info$peak_id %in% each_tp$V4, ]
    peak_distribution <- sapply(strsplit(peak_info_details$data_id, split='_'), getElement, 2)
    peak_distribution <- table(peak_distribution) |> as.data.frame() 
    colnames(peak_distribution) <- c('count', 'freq')
    return(peak_distribution)

})

names(peak_distribution) <- names(tp_list)
pkd <- do.call(rbind, peak_distribution)

pkd$count <- sapply(strsplit(names(peak_distribution), split='_'), function(each){
    getElement(each, 2)
})

pkd <- pkd %>% dplyr::mutate(count=as.numeric(count)) %>% dplyr::arrange(count)

chip_distribution <- sapply(chip_peak_file_split, nrow) |> as.data.frame()
chip_distribution$count <- row.names(chip_distribution)
colnames(chip_distribution) <- c('freq_chip', 'count')

pkd <- merge(pkd, chip_distribution, by='count')

row.names(pkd) <- pkd$count
pkd$count <- NULL

colnames(pkd) <- c('chip+motif', 'chip')

bp <- barplot(t(pkd), ylab='Number of peaks', xlab='Number of individuals intersected', beside=T, legend = colnames(pkd), ylim=c(0, max(pkd)), font.axis=2)

#bp <- axis(1, at=1:length(bedfiles_names), labels=F, tick=T)
#text(x=bp, y=pkd$freq, label=pkd$freq, cex=0.8, pos=4, srt=90) 
mtext('Freedman (FOXA1): Number of ChIP-seq peaks overlapped with true, known motifs')

# ====
true_positives <- lapply(seq_along(tp_list), function(i){

    each_tp <- tp_list[[i]][, 1:3]
    each_tp$binding_count <- getElement(strsplit(names(tp_list)[i], split='_'), 1)[2]

    return(each_tp)

})

true_positives <- do.call(rbind, true_positives)
colnames(true_positives) <- c('chr', 'start', 'end', 'binding_counts')

true_positives <- true_positives %>% dplyr::arrange(chr, start)


# false positives
false_positives <- chip_peak_file[, 1:3]

fp <- false_positives[!(chrom %in% true_positives$chr & start %in% true_positives$start & end %in% true_positives$end), ]



# cl <- parallel::makeCluster(25)
# registerDoParallel(cl)
# parallel::clusterExport(cl, c('chip_regions_split', 'motif_instances'))

# st_time <- proc.time()
# out_peaks <- lapply(chip_regions_split, function(each_chip){

#     parallel::clusterExport(cl, c('each_chip'), envir=environment())

#     each_chip_newpeaks <- foreach(i = icount(n_row), .combine='rbind') %dopar% {
#         print(paste0("row ", i, " of ", nrow(each_chip)))
#         w <- which(motif_instances$PositionID == each_chip[i,4]) #match peakID, gauranteed to have a match 
#         #only choose 1 motif, take the one with the top score
#         w1 <- w[match(max(motif_instances$MotifScore[w]), motif_instances$MotifScore[w])]
#         if (motif_instances$Strand[w1] == "+"){
#             centerpeak <- motif_instances$Offset[w1] + floor(0.5*(nchar(motif_instances$Sequence[w1]))) #this is relative to start 
#         }else{
#             centerpeak <- motif_instances$Offset[w1] - floor(0.5*(nchar(motif_instances$Sequence[w1])))
#         }
#         newstart <- each_chip$starts[i]+centerpeak #previous version of impact -windowsize/2 #establish left bound of peaks
#         newend <- newstart + 1 #previous version of impact +windowsize #establish right bound of peaks 
#         out_append <- c(as.character(each_chip$seqnames[i]), newstart, newend, as.character(each_chip$names[i]), as.character(each_chip$scores[i]), as.character(each_chip$strands[i]))
#         #chip_regions_newpeaks[i, ] <- out_append
#     }

#     each_chip_newpeaks

# })

# stopCluster(cl)
# en_time <- proc.time()


# chip_regions_newpeaks <- foreach(i = icount(n_row), .combine='rbind') %dopar% {
#     print(paste0("row ", i, " of ", nrow(chip_regions)))
#     w <- which(motif_instances$PositionID == chip_regions[i,4]) #match peakID, gauranteed to have a match 
#     #only choose 1 motif, take the one with the top score
#     w1 <- w[match(max(motif_instances$MotifScore[w]), motif_instances$MotifScore[w])]
#     if (motif_instances$Strand[w1] == "+"){
#         centerpeak <- motif_instances$Offset[w1] + floor(0.5*(nchar(motif_instances$Sequence[w1]))) #this is relative to start 
#     }else{
#         centerpeak <- motif_instances$Offset[w1] - floor(0.5*(nchar(motif_instances$Sequence[w1])))
#     }
# 	newstart <- chip_regions$starts[i]+centerpeak #previous version of impact -windowsize/2 #establish left bound of peaks
#     newend <- newstart + 1 #previous version of impact +windowsize #establish right bound of peaks 
#     out_append <- c(as.character(chip_regions$seqnames[i]), newstart, newend, as.character(chip_regions$names[i]), as.character(chip_regions$scores[i]), as.character(chip_regions$strands[i]))
#     #chip_regions_newpeaks[i, ] <- out_append
# }

# stopCluster(cl)

# en_time <- proc.time()

# en_time - st_time

# 2+2