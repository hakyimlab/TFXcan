
# This script defines true positive and true negative motifs for the Kawakami data


rm(list=ls())

setwd('/projects/covid-ct/imlab/users/temi/projects/TFXcan/scripts/')
source('/projects/covid-ct/imlab/users/temi/projects/TFXcan/scripts/custom-impact.R')

library(glue)
library(GenomicRanges)
library(reticulate)
library(R.utils)
library(data.table)
library(tidyverse)

# pick a transcription factor
TF <- 'FOXA1' #
cell_line <- 'LuCaP' # a prostrate cancer cell line
foxa1_motif <- 'foxa1.lncap.motif'

# where is the TF chip-seq data, as well as other directories
work_dir <- '/projects/covid-ct/imlab/users/temi/projects/TFXcan'
#kawakami_dir <- glue('{work_dir}/data/kawakami-human')
kawakami_impact_dir <- glue('{work_dir}/data/kawakami-human')
kawakami_impact_split_dir <- glue('{kawakami_impact_dir}/kawakami-impact-split')
chip_data_dir <- glue('{kawakami_impact_split_dir}/chip_data_files')
common_files_dir <- glue('{work_dir}/data/common-files/{TF}')
homer_dir <- glue('~/miniconda3/share/homer')
project_dir <- '../'
#kawakami_impact_split_dir <- glue('../processed-data/{TF}/{cell_line}/freedman')

# ensure that the kawakami_impact_split_dir is created
lapply(c(kawakami_impact_dir, kawakami_impact_split_dir, chip_data_dir, common_files_dir), function(each_dir){

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
chip_peak_file <- data.table::fread(glue('{kawakami_impact_dir}/{TF}_intersected_files.bed'), select=chip_peak_headers)
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
write.table(peak_info, glue("{kawakami_impact_split_dir}/{TF}_{cell_line}_true_positive_peak_info.txt"), sep = "\t", quote = F, row.names=F)

# split by num
chip_peak_file_split <- split(chip_peak_file, f=chip_peak_file$num)

run_output <- lapply(X=seq_along(chip_peak_file_split), FUN=annotate_with_impact, output_dir=kawakami_impact_split_dir, common_dir=common_files_dir)
















# out <- lapply(seq_along(chip_peak_file_split), function(n){


# })

# read in the true positives

allfiles <- list.files(glue("{kawakami_impact_split_dir}"))
positive_files <- allfiles[endsWith(x=allfiles, suffix='FOXA1_LuCaP_train_test_positive_bed.txt')]
#peak_info <- allfiles[endsWith(x=allfiles, suffix='FOXA1_LuCaP_true_positive_peak_info.txt')]

# create enformer beds =====
tp_list <- lapply(positive_files, function(filename){
    read.table(glue("{kawakami_impact_split_dir}/{filename}"))
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
    peak_distribution <- sapply(strsplit(peak_info_details$kawakami_id, split='_'), getElement, 2)
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

bp <- barplot(pkd$freq, names.arg=pkd$count, 
    ylim=c(0, max(pkd$freq + 30000)), 
    ylab='No. of overlapping true ChIP peaks', xlab='Number of files')

#bp <- axis(1, at=1:length(bedfiles_names), labels=F, tick=T)
text(x=bp, y=pkd$freq, label=pkd$freq, cex=0.8, pos=4, srt=90) 
mtext('Kawakami: Number of ChIP-seq peaks overlapped with known motifs')





































# out <- lapply(seq_along(chip_peak_file_split), function(n){

#     cat(glue('Starting {n} \n'))

#     chip_peak_file <- chip_peak_file_split[[n]]

#     peak_info <- cbind(chip_peak_file$id, paste('peak', 1:length(chip_peak_file$id), sep='')) |> as.data.frame()
#     colnames(peak_info) <- c('kawakami_id', 'peak_id')
#     write.table(peak_info, glue("{kawakami_impact_split_dir}/common_{n}_{TF}_{cell_line}_true_positive_peak_info.txt"), sep = "\t", quote = F)

#     motif <- paste0(TF, '_', cell_line)
#     TF_BED <- unique(cbind(chip_peak_file[, c(1:3)], 0, motif, "+"))
#     TF_BED[, 4] <- paste0("peak", seq(1,nrow(TF_BED), 1))

#     colnames(TF_BED)[1:6] <- c('chr','start','end','id','score','strand') #score may be number of reads
#     bed_Granges <- with(TF_BED, GRanges(chr, IRanges(start,end), strand, score, id = id))
#     bed_Granges <- reduce(bed_Granges)

#     tf_bed <- cbind(as.character(seqnames(bed_Granges)), start(bed_Granges), end(bed_Granges), paste0("peak", seq(1,length(bed_Granges),1)),0, "+")
#     write.table(tf_bed, glue("{kawakami_impact_split_dir}/common_{n}_{TF}_{cell_line}_train_df.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

#     #homer find instances of motifs
#     if(!file.exists(glue('{kawakami_impact_split_dir}/common_{n}_{TF}_{cell_line}_findinstances.txt'))){
#         cmd <- glue('perl {homer_dir}/bin/findMotifsGenome.pl {kawakami_impact_split_dir}/common_{n}_{TF}_{cell_line}_train_df.txt {homer_dir}/data/genomes/hg19 {project_dir}/processed-data/homer_output/ -size given -find {kawakami_impact_split_dir}/{TF}_{cell_line}.motif > {kawakami_impact_split_dir}/common_{n}_{TF}_{cell_line}_findinstances.txt')
#         system(cmd)
#     } else {
#         print(glue('{kawakami_impact_split_dir}/common_{n}_{TF}_{cell_line}_findinstances.txt exists.'))
#     }

#     # STEP 2 : prepare the positive and negative sets
#     motif_instances <- read.table(glue('{kawakami_impact_split_dir}/common_{n}_{TF}_{cell_line}_findinstances.txt'), sep = "\t", header = T, stringsAsFactors = FALSE)
#     chip_regions <- read.table(glue('{kawakami_impact_split_dir}/common_{n}_{TF}_{cell_line}_train_df.txt'), sep='\t', header=F, stringsAsFactors=F) #as.data.frame(tf_bed) # or read.table(glue('{kawakami_impact_split_dir}/{TF}_train_df.txt'), sep='\t', header=T, stringsAsFactors=F)
#     motif_library <- read.table(glue('{kawakami_impact_split_dir}/{TF}_{cell_line}_motif_threshold_info.txt'), sep = "\t", header = F, stringsAsFactors = FALSE)

#     colnames(chip_regions) <- c("seqnames","starts","ends","names","scores","strands")
#     #select ChIP peaks with a matching motif under them.
#     chip_regions_full <- chip_regions
#     m <- match(chip_regions[,4], motif_instances[,1]) # what chips are in motifs data?
#     chip_regions <- chip_regions_full[which(is.na(m) == F),]

#     chip_regions_newpeaks <- matrix(0,nrow(chip_regions),6) #dummy var number of rows 

#     if(!file.exists(glue("{kawakami_impact_split_dir}/common_{n}_{TF}_{cell_line}_train_test_positive_bed.txt"))){
#         for (i in 1:nrow(chip_regions)){
#             print(paste0("row ", i, " of ", nrow(chip_regions)))
#             w <- which(motif_instances$PositionID == chip_regions[i,4]) #match peakID, gauranteed to have a match 
#             #only choose 1 motif, take the one with the top score
#             w1 <- w[match(max(motif_instances$MotifScore[w]), motif_instances$MotifScore[w])]
#             if (motif_instances$Strand[w1] == "+"){
#                 centerpeak <- motif_instances$Offset[w1] + floor(0.5*(nchar(motif_instances$Sequence[w1]))) #this is relative to start 
#             }else{
#                 centerpeak <- motif_instances$Offset[w1] - floor(0.5*(nchar(motif_instances$Sequence[w1])))
#             }
#             newstart <- chip_regions$starts[i]+centerpeak #previous version of impact -windowsize/2 #establish left bound of peaks
#             newend <- newstart + 1 #previous version of impact +windowsize #establish right bound of peaks 
#             chip_regions_newpeaks[i,1] <- as.character(chip_regions$seqnames[i])
#             chip_regions_newpeaks[i,2] <- newstart
#             chip_regions_newpeaks[i,3] <- newend
#             chip_regions_newpeaks[i,4] <- as.character(chip_regions$names[i]) #preserve old information
#             chip_regions_newpeaks[i,5] <- as.character(chip_regions$scores[i])
#             chip_regions_newpeaks[i,6] <- as.character(chip_regions$strands[i])
#         }

#         write.table(chip_regions_newpeaks, glue("{kawakami_impact_split_dir}/common_{n}_{TF}_{cell_line}_train_df_newpeaks.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
    
#     } else {
#         print(glue('train df newpeaks file exists. Moving on...\n'))
#     }

    
#     chip_regions_newpeaks <- read.table(glue('{kawakami_impact_split_dir}/common_{n}_{TF}_{cell_line}_train_df_newpeaks.txt'), sep = "\t", header = T, stringsAsFactors = FALSE)

#     chip_regions_newpeaks_realdf <- data.frame(seqnames=chip_regions_newpeaks[,1],
#                                         starts=as.numeric(chip_regions_newpeaks[,2]),
#                                         ends=as.numeric(chip_regions_newpeaks[,3]),
#                                         names=chip_regions_newpeaks[,4],
#                                         scores=0,
#                                         strands="+")

#     set.seed(2022)
#     s <- sample(seq(1,nrow(chip_regions_newpeaks_realdf),1), nrow(chip_regions_newpeaks_realdf),replace = F)

#     chip_regions_newpeaks_realdf <- chip_regions_newpeaks_realdf[s,]
#     #put training and test together! 
#     #use partition function later to make it easier to sample many times! 
#     train_test_peaks <- unique(chip_regions_newpeaks_realdf)
#     colnames(train_test_peaks) <- c("seqnames","starts","ends","names","scores","strands")

#     #clean up 
#     if (length(which(train_test_peaks$start < 0)) > 0){
#         train_test_peaks$start[which(train_test_peaks$start < 0)] <- 0
#     }

#     t <- train_test_peaks[,3]-train_test_peaks[,2]
#     w <- which(t %% 2 != 0) #odd
#     if (length(w) > 0){
#         train_test_peaks[w,3] <- train_test_peaks[w,3] + 1
#     }

#     w1 <- which(train_test_peaks$chr == "chrM")
#     w2 <- which(nchar(train_test_peaks$chr) > 5)
#     if (length(c(w1,w2)) > 0){
#         train_test_peaks <- train_test_peaks[-c(w1,w2),]
#     }

#     write.table(train_test_peaks, glue("{kawakami_impact_split_dir}/common_{n}_{TF}_{cell_line}_train_test_positive_bed.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)


#     # negative sets

#     print('Preparing negative sets \n')
#     train_df_newpeaks <- read.table(glue("{kawakami_impact_split_dir}/common_{n}_{TF}_{cell_line}_train_df_newpeaks.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)

#     MotifThreshold <- read.table(glue("{kawakami_impact_split_dir}/{TF}_{cell_line}_motif_threshold_info.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)

#     motifs <- MotifThreshold[1,1]
#     instscan <- read.table(glue("{kawakami_impact_split_dir}/{TF}_{cell_line}_scanMotifsGenomewide_sort_15000.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
#     instscan <- instscan[order(instscan[,6], decreasing = T),]   
#     train_df <- read.table(glue("{kawakami_impact_split_dir}/common_{n}_{TF}_{cell_line}_train_df.txt"), sep = "\t", header = F, stringsAsFactors = FALSE) #old version: ChIPbed1.txt
#     colnames(train_df) <- c('chr','start','end','id','score','strand') #score may be number of reads
#     Granges <- with(train_df, GRanges(chr, IRanges(start,end), strand, score, id = id))
#     assign(paste0(motifs,"_Granges"),Granges)

#     allinstscan <- instscan
#     sapply_instscan <- sapply(1:nrow(allinstscan), function(x) strsplit(allinstscan[x,1], "-")[[1]][1])
#     allinstscan[,1] <- sapply_instscan
#     allinstscan[,6] <- allinstscan[,6]/MotifThreshold[1,2]

#     w <- which(allinstscan[,6] < 1) 
#     if(length(w)>0){allinstscan <- allinstscan[-w,]}

#     motifmatrix <- matrix(0,nrow(MotifThreshold),3)
#     #col1: motif name with cell type (1 per row)
#     #col2: degenerate motif name (used in scanMotifsGenome.pl), get with table(instances_scan[,1])
#     #col3: Granges files
#     motifmatrix[,1] <- motifs 
#     motifmatrix[,2] <- motifs
#     motifmatrix[,3] <- paste0(motifs,"_Granges")

#     #windowsize <- 380
#     print(paste0("working on ",motifmatrix[1,1]))
#     w_mcp <- which(train_df_newpeaks[,5] == motifmatrix[1,2]) #motifcentered peaks
#     motifcentered_specificTF <- train_df_newpeaks[w_mcp,]
#     colnames(motifcentered_specificTF) <- c('chr','start','end','id','score','strand') 
#     Granges_obj_mc <- with(motifcentered_specificTF, GRanges(chr, IRanges(start,end), strand, score, id = id))

#     w_inst <- which(allinstscan[,1] == motifmatrix[1,2]) #instances pertaining to specific motif
#     instances_specificTF <- allinstscan[w_inst,]
#     inst_bed <- instances_specificTF[,c(2,3,4)]
#     # if(dim(inst_bed)[1] == 0){
#     #     inst_bed <- rbind(inst_bed, c(0, 0, 0))
#     # }
#     inst_bed <- cbind(inst_bed, motifmatrix[1,1], instances_specificTF[,6], "+") 
#     #region as wide as a motif 
#     motif_width <- as.numeric(inst_bed[1,3]) - as.numeric(inst_bed[,2])
#     inst_bed[,2] <- as.numeric(inst_bed[,2]) + floor(motif_width/2)
#     inst_bed[,3] <- as.numeric(inst_bed[,2]) + 1

#     #inst_bed[,3] <- as.numeric(inst_bed[,3])#previous version of impact: +(windowsize-10)/2
#     #inst_bed[,2] <- as.numeric(inst_bed[,2])#previous version of impact: -(windowsize-10)/2
#     inst_bed <- unique(inst_bed)
#     colnames(inst_bed) <- c('chr','start','end','id','score','strand') 
#     inst_bed_Granges <- with(inst_bed, GRanges(chr, IRanges(start,end), strand, score, id = id)) #instances Granges obj

#     int_1 <- GenomicRanges::intersect(inst_bed_Granges, get(motifmatrix[1,3]))
#     int_2 <- GenomicRanges::intersect(inst_bed_Granges, Granges_obj_mc)

#     ms_1 <- match(start(int_1)-1, inst_bed$start)
#     me_1 <- match(end(int_1), inst_bed$end)
#     a <- c(ms_1,me_1)
#     a <- unique(a[is.na(a) == F])

#     ms_2 <- match(start(int_2)-1, inst_bed$start)
#     me_2 <- match(end(int_2), inst_bed$end)
#     b <- c(ms_2,me_2)
#     b <- unique(b[is.na(b) == F])
#     ab <- unique(c(a,b))

#     if (length(ab) > 0){
#         noChIP <- inst_bed[-ab,] #good negative peaks
#     }else{
#         noChIP <- inst_bed
#     }

#     write.table(noChIP, glue("{kawakami_impact_split_dir}/common_{n}_{TF}_{cell_line}_noChIP.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
#     colnames(noChIP) <- c("seqnames","starts","ends","names","scores","strands")

#     #shuffle 
#     s <- sample(seq(1,nrow(noChIP),1),nrow(noChIP),replace = F)
#     noChIP <- noChIP[s,]
#     noChIP <- unique(noChIP)

#     ##clean up 
#     if (length(which(noChIP$start < 0))){
#     noChIP$start[which(noChIP$start < 0)] <- 0
#     }

#     t <- noChIP[,3]-noChIP[,2]
#     w <- which(t %% 2 != 0) #odd
#     if (length(w) > 0){noChIP[w,3] <- noChIP[w,3] + 1}

#     w1 <- which(noChIP$chr == "chrM")
#     w2 <- which(nchar(noChIP$chr) > 5)
#     if (length(c(w1,w2)) > 0){noChIP <- noChIP[-c(w1,w2),]}

#     write.table(noChIP, glue("{kawakami_impact_split_dir}/common_{n}_{TF}_{cell_line}_train_test_negative_bed.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

#     cat(glue('Ended {n} \n'))

# })















