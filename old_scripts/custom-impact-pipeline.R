
# This pipeline works with Cistrome data
# Author: Temi
# Date: Circa July 23 2022


setwd('/projects/covid-ct/imlab/users/temi/projects/TFXcan/scripts/')

rm(list=ls())

library(glue)
library(GenomicRanges)
library(reticulate)
library(R.utils)
library(data.table)


# I very much prefer to extract the bed files this way because of the way the files have been archived i.e. zip:gzip...
# and using the python code I wrote makes it easy - also a good way to learn how to use reticulate ^^
reticulate::use_python("~/miniconda3/envs/r-env/bin/python3.10", required = T)
reticulate::py_config()
reticulate::source_python('./load-cistrome-data.py')

py$a <- as.integer(c(1, 8))
py$b <- 'TF'

#pbfile <- py$get_bed_file(py$a, py$b)


# pick a transcription factor
TF <- 'GATA3'

# where is the chip-seq data, as well as other directories
cistrome_dir <- '/projects/covid-ct/imlab/data/cistrome/compressed'
homer_dir <- '~/miniconda3/share/homer'
project_dir <- '../'
output_dir <- glue('../processed-data/{TF}')

list.files(cistrome_dir)

# make sure that the TF's Chip data and motif data exist
# check for the motif in Homer
if(!file.exists(glue('{homer_dir}/data/knownTFs/motifs/{tolower(TF)}.motif'))){
    print(glue('No motif information for {TF} available'))
} else {
    print(glue('Motif information for {TF} available. Moving on.'))
}

# check that the chip-seq data is avaialable
# hc_info <- data.table::fread(glue('{cistrome_dir}/human_ca_full_QC.txt'))
# hf_info <- data.table::fread(glue('{cistrome_dir}/human_factor_full_QC.txt'))
# hm_info <- data.table::fread(glue('{cistrome_dir}/human_hm_full_QC.txt'))
hf_info <- data.table::fread(glue('{cistrome_dir}/human_factor_full_QC.txt'))

mtdata_names <- c('Cell_line', 'Cell_type', 'Tissue_type')

# https://genome.ucsc.edu/FAQ/FAQformat.html#format12

if(!TF %in% hf_info$Factor){
    print(glue('No ChIP information for {TF} available'))
} else {
    print(glue('ChIP information for {TF} available. Loading the data...'))

    TF_INFO <- hf_info[hf_info$Factor == TF, c('DCid', 'Cell_line', 'Cell_type', 'Tissue_type')] # across different cell lines and tissues
    py$TF_DCID <- as.integer(TF_INFO$DCid)
    py$b <- 'TF'

    TF_CHIP <- py$get_bed_file(py$TF_DCID, py$b)
    # convert each object to an R object and retain only the first 3 columns
    TF_CHIP <- lapply(seq_along(TF_CHIP), function(i){
        each_chip <- TF_CHIP[[i]]
        each_chip_dcid <- as.numeric(names(TF_CHIP)[i])

        temp <- reticulate::py_to_r(each_chip)
        each_chip_extra_info <- TF_INFO[TF_INFO$DCid == each_chip_dcid, c('Cell_line', 'Cell_type', 'Tissue_type')]
        cbind(temp[, 1:3], each_chip_extra_info)
    })
    TF_CHIP <- do.call(rbind, TF_CHIP)
}


# TF_INFO <- hf_info[hf_info$Factor == TF, c('DCid', 'Cell_line', 'Cell_type', 'Tissue_type')] # across different cell lines and tissues
# py$TF_DCID <- as.integer(TF_INFO$DCid)
# py$b <- 'TF'

# TF_CHIP <- py$get_bed_file(py$TF_DCID, py$b)
# # convert each object to an R object and retain only the first 3 columns
# TF_CHIP <- lapply(seq_along(TF_CHIP), function(i){
#     each_chip <- TF_CHIP[[i]]
#     each_chip_dcid <- as.numeric(names(TF_CHIP)[i])

#     temp <- reticulate::py_to_r(each_chip)
#     each_chip_extra_info <- TF_INFO[TF_INFO$DCid == each_chip_dcid, c('Cell_line', 'Cell_type', 'Tissue_type')]
#     cbind(temp[, 1:3], each_chip_extra_info)
# })
# TF_CHIP <- do.call(rbind, TF_CHIP)


if(!dir.exists(output_dir)){
    dir.create(output_dir)
} else {
    print(glue('{output_dir} exists.'))
}

# TF_DCID <- hf_info$DCid[hf_info$Factor == TF] # across different cell lines and tissues

# py$STAT3_DCID <- as.integer(TF_DCID)
# py$b <- 'TF'

# STAT3_CHIP <- py$get_bed_file(py$STAT3_DCID, py$b)
# # convert each object to an R object and retain only the first 3 columns
# STAT3_CHIP <- lapply(STAT3_CHIP, function(each_chip){
#     temp <- reticulate::py_to_r(each_chip)
#     temp[, 1:3]
# })
# STAT3_CHIP <- do.call(rbind, STAT3_CHIP)

with(hf_info[hf_info$Factor == TF, ], setNames(as.character(Cell_type), DCid))


## TO-DO
# 1. I could merge the chip files that have the same cell lines, types and so on
# 2. 
# 3. 
# 4. 

length(TF_CHIP) # 42 different TF signals

motif <- paste0(TF, '_many_cell_lines')
TF_BED <- unique(cbind(TF_CHIP[, c(1:3)], 0, motif, "+", TF_CHIP[4:ncol(TF_CHIP)]))
TF_BED[, 4] <- paste0("peak", seq(1,nrow(TF_BED), 1))

colnames(TF_BED)[1:6] <- c('chr','start','end','id','score','strand') #score may be number of reads
bed_Granges <- with(TF_BED, GRanges(chr, IRanges(start,end), strand, score, id = id))
bed_Granges <- reduce(bed_Granges)
#elementMetadata(bed_Granges) <- TF_BED[, mtdata_names]

#ty <- cbind(bed_Granges, TF_BED[, mtdata_names])


tf_bed <- cbind(as.character(seqnames(bed_Granges)), start(bed_Granges), end(bed_Granges), paste0("peak",seq(1,length(bed_Granges),1)),0, "+")

write.table(tf_bed, glue("{output_dir}/{TF}_train_df.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

if(!file.exists(glue('{output_dir}/{TF}.motif'))){
    # copy the motif file from homer
    cmd <- glue('cp {homer_dir}/data/knownTFs/motifs/{tolower(TF)}.motif {output_dir}/{TF}.motif')
    system(cmd)
} else {
    print(glue('{output_dir}/{TF}.motif exists'))
}


# 
pwm <- read.table(glue('{output_dir}/{TF}.motif'), sep = "\t", header = F, stringsAsFactors = FALSE, fill=T)
pwm[1, 4] <- 'PH'
pwm[1, 2] <- glue('{TF}_all') # NOTE THIS PART CAREFULLY
write.table(pwm[, 1:4], glue('{output_dir}/{TF}.motif'), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

#Input Motifs
MotifThreshvals <- pwm[which(pwm[,4]=="PH"), 3]
MotifThreshnames <- pwm[which(pwm[,4]=="PH"), 2]
MotifThreshinfo <- cbind(MotifThreshnames, MotifThreshvals)

write.table(MotifThreshinfo, glue('{output_dir}/{TF}_motif_threshold_info.txt'), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

# create the genome wide motif files
if(!file.exists(glue('{output_dir}/{TF}_scanMotifsGenomewide_sort_15000.txt'))){

    cmd <- glue('perl ~/miniconda3/share/homer/bin/scanMotifGenomeWide.pl {output_dir}/{TF}.motif ~/miniconda3/share/homer/data/genomes/hg19 > {output_dir}/{TF}_scanMotifsGenomewide.txt')
    system(cmd)
    cmd <- glue("sort -t $'\t' -k6,6rn {output_dir}/{TF}_scanMotifsGenomewide.txt > {output_dir}/{TF}_scanMotifsGenomewide_sort.txt")
    system(cmd)
    cmd <- glue('head -n 15000 {output_dir}/{TF}_scanMotifsGenomewide_sort.txt > {output_dir}/{TF}_scanMotifsGenomewide_sort_15000.txt')
    system(cmd)
 
} else {
    print(glue('{output_dir}/{TF}_scanMotifsGenomewide_sort_15000.txt exists.'))
}

#homer find instances of motifs

if(!file.exists(glue('{output_dir}/{TF}_findinstances.txt'))){
    cmd <- glue('perl {homer_dir}/bin/findMotifsGenome.pl {output_dir}/{TF}_train_df.txt {homer_dir}/data/genomes/hg19 {project_dir}/processed-data/homer_output/ -size given -find {output_dir}/{TF}.motif > {output_dir}/{TF}_findinstances.txt')
    system(cmd)
} else {
    print(glue('{output_dir}/{TF}_findinstances.txt exists.'))
}

# STEP 2 : prepare the positive and negative sets
motif_instances <- read.table(glue('{output_dir}/{TF}_findinstances.txt'), sep = "\t", header = T, stringsAsFactors = FALSE)
chip_regions <- read.table(glue('{output_dir}/{TF}_train_df.txt'), sep='\t', header=F, stringsAsFactors=F) #as.data.frame(tf_bed) # or read.table(glue('{output_dir}/{TF}_train_df.txt'), sep='\t', header=T, stringsAsFactors=F)
motif_library <- read.table(glue('{output_dir}/{TF}_motif_threshold_info.txt'), sep = "\t", header = F, stringsAsFactors = FALSE)

colnames(chip_regions) <- c("seqnames","starts","ends","names","scores","strands")
#select ChIP peaks with a matching motif under them.
chip_regions_full <- chip_regions
m <- match(chip_regions[,4], motif_instances[,1]) 
chip_regions <- chip_regions_full[which(is.na(m) == F),]

chip_regions_newpeaks <- matrix(0,nrow(chip_regions),6) #dummy var number of rows 

for (i in 1:nrow(chip_regions)){
    print(paste0("row ", i, " of ", nrow(chip_regions)))
    w <- which(motif_instances$PositionID == chip_regions[i,4]) #match peakID, gauranteed to have a match 
    #only choose 1 motif, take the one with the top score
    w1 <- w[match(max(motif_instances$MotifScore[w]), motif_instances$MotifScore[w])]
    if (motif_instances$Strand[w1] == "+"){
        centerpeak <- motif_instances$Offset[w1] + floor(0.5*(nchar(motif_instances$Sequence[w1]))) #this is relative to start 
    }else{
        centerpeak <- motif_instances$Offset[w1] - floor(0.5*(nchar(motif_instances$Sequence[w1])))
    }
	newstart <- chip_regions$starts[i]+centerpeak #previous version of impact -windowsize/2 #establish left bound of peaks
    newend <- newstart + 1 #previous version of impact +windowsize #establish right bound of peaks 
    chip_regions_newpeaks[i,1] <- as.character(chip_regions$seqnames[i])
    chip_regions_newpeaks[i,2] <- newstart
    chip_regions_newpeaks[i,3] <- newend
    chip_regions_newpeaks[i,4] <- as.character(chip_regions$names[i]) #preserve old information
    chip_regions_newpeaks[i,5] <- as.character(chip_regions$scores[i])
    chip_regions_newpeaks[i,6] <- as.character(chip_regions$strands[i])
}


write.table(chip_regions_newpeaks, glue("{output_dir}/{TF}_train_df_newpeaks.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

chip_regions_newpeaks_realdf <- data.frame(seqnames=chip_regions_newpeaks[,1],
                                       starts=as.numeric(chip_regions_newpeaks[,2]),
                                       ends=as.numeric(chip_regions_newpeaks[,3]),
                                       names=chip_regions_newpeaks[,5],
                                       scores=0,
                                       strands="+")

set.seed(2022)
s <- sample(seq(1,nrow(chip_regions_newpeaks_realdf),1), nrow(chip_regions_newpeaks_realdf),replace = F)

chip_regions_newpeaks_realdf <- chip_regions_newpeaks_realdf[s,]
#put training and test together! 
#use partition function later to make it easier to sample many times! 
train_test_peaks <- unique(chip_regions_newpeaks_realdf)
colnames(train_test_peaks) <- c("seqnames","starts","ends","names","scores","strands")

#clean up 
if (length(which(train_test_peaks$start < 0)) > 0){
  train_test_peaks$start[which(train_test_peaks$start < 0)] <- 0
}

t <- train_test_peaks[,3]-train_test_peaks[,2]
w <- which(t %% 2 != 0) #odd
if (length(w) > 0){train_test_peaks[w,3] <- train_test_peaks[w,3] + 1}

w1 <- which(train_test_peaks$chr == "chrM")
w2 <- which(nchar(train_test_peaks$chr) > 5)
if (length(c(w1,w2)) > 0){train_test_peaks <- train_test_peaks[-c(w1,w2),]}

write.table(train_test_peaks, glue("{output_dir}/{TF}_train_test_positive_bed.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

# STEP 3 : prepare the negative sets

train_df_newpeaks <- read.table(glue("{output_dir}/{TF}_train_df_newpeaks.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)

MotifThreshold <- read.table(glue("{output_dir}/{TF}_motif_threshold_info.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)

motifs <- MotifThreshold[1,1]
instscan <- read.table(glue("{output_dir}/{TF}_scanMotifsGenomewide_sort_15000.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
instscan <- instscan[order(instscan[,6], decreasing = T),]   
train_df <- read.table(glue("{output_dir}/{TF}_train_df.txt"), sep = "\t", header = F, stringsAsFactors = FALSE) #old version: ChIPbed1.txt
colnames(train_df) <- c('chr','start','end','id','score','strand') #score may be number of reads
Granges <- with(train_df, GRanges(chr, IRanges(start,end), strand, score, id = id))
assign(paste0(motifs,"_Granges"),Granges)


allinstscan <- instscan
sapply_instscan <- sapply(1:nrow(allinstscan), function(x) strsplit(allinstscan[x,1], "-")[[1]][1])
allinstscan[,1] <- sapply_instscan
allinstscan[,6] <- allinstscan[,6]/MotifThreshold[1,2]

w <- which(allinstscan[,6] < 1) 
if(length(w)>0){allinstscan <- allinstscan[-w,]}

motifmatrix <- matrix(0,nrow(MotifThreshold),3)
#col1: motif name with cell type (1 per row)
#col2: degenerate motif name (used in scanMotifsGenome.pl), get with table(instances_scan[,1])
#col3: Granges files
motifmatrix[,1] <- motifs 
motifmatrix[,2] <- motifs
motifmatrix[,3] <- paste0(motifs,"_Granges")

#windowsize <- 380
print(paste0("working on ",motifmatrix[1,1]))
w_mcp <- which(train_df_newpeaks[,5] == motifmatrix[1,2]) #motifcentered peaks
motifcentered_specificTF <- train_df_newpeaks[w_mcp,]
colnames(motifcentered_specificTF) <- c('chr','start','end','id','score','strand') 
Granges_obj_mc <- with(motifcentered_specificTF, GRanges(chr, IRanges(start,end), strand, score, id = id))

w_inst <- which(allinstscan[,1] == motifmatrix[1,2]) #instances pertaining to specific motif
instances_specificTF <- allinstscan[w_inst,]
inst_bed <- instances_specificTF[,c(2,3,4)]
inst_bed <- cbind(inst_bed, motifmatrix[1,1], instances_specificTF[,6], "+") 
#region as wide as a motif 
motif_width <- as.numeric(inst_bed[1,3]) - as.numeric(inst_bed[,2])
inst_bed[,2] <- as.numeric(inst_bed[,2]) + floor(motif_width/2)
inst_bed[,3] <- as.numeric(inst_bed[,2]) + 1

#inst_bed[,3] <- as.numeric(inst_bed[,3])#previous version of impact: +(windowsize-10)/2
#inst_bed[,2] <- as.numeric(inst_bed[,2])#previous version of impact: -(windowsize-10)/2
inst_bed <- unique(inst_bed)
colnames(inst_bed) <- c('chr','start','end','id','score','strand') 
inst_bed_Granges <- with(inst_bed, GRanges(chr, IRanges(start,end), strand, score, id = id)) #instances Granges obj

int_1 <- intersect(inst_bed_Granges, get(motifmatrix[1,3]))
int_2 <- intersect(inst_bed_Granges, Granges_obj_mc)

ms_1 <- match(start(int_1)-1, inst_bed$start)
me_1 <- match(end(int_1), inst_bed$end)
a <- c(ms_1,me_1)
a <- unique(a[is.na(a) == F])

ms_2 <- match(start(int_2)-1, inst_bed$start)
me_2 <- match(end(int_2), inst_bed$end)
b <- c(ms_2,me_2)
b <- unique(b[is.na(b) == F])
ab <- unique(c(a,b))

if (length(ab) > 0){
    noChIP <- inst_bed[-ab,] #good negative peaks
}else{
    noChIP <- inst_bed
}

write.table(noChIP, glue("{output_dir}/{TF}_noChIP.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
colnames(noChIP) <- c("seqnames","starts","ends","names","scores","strands")

#shuffle 
s <- sample(seq(1,nrow(noChIP),1),nrow(noChIP),replace = F)
noChIP <- noChIP[s,]
noChIP <- unique(noChIP)

##clean up 
if (length(which(noChIP$start < 0))){
  noChIP$start[which(noChIP$start < 0)] <- 0
}

t <- noChIP[,3]-noChIP[,2]
w <- which(t %% 2 != 0) #odd
if (length(w) > 0){noChIP[w,3] <- noChIP[w,3] + 1}

w1 <- which(noChIP$chr == "chrM")
w2 <- which(nchar(noChIP$chr) > 5)
if (length(c(w1,w2)) > 0){noChIP <- noChIP[-c(w1,w2),]}

write.table(noChIP, glue("{output_dir}/{TF}_train_test_negative_bed.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

# create enformer beds =====
tp <- read.table(glue("{output_dir}/{TF}_train_test_positive_bed.txt"))
tn <- read.table(glue("{output_dir}/{TF}_train_test_negative_bed.txt"))

gsize <- read.table(glue('{output_dir}/../hg19.sizes'))


expand_region <- function(dt, genome_size, both=128, left=NULL, right=NULL, create_granges=F, save_as_bed=T){
    colnames(dt) <- c('chr','motif_center_start','motif_center_end','id','score','strand')
    colnames(genome_size) <- c('chr', 'size')
    
    if(both){
        left <- both/2
        right <- both/2
    }
    
    center <- (dt[, 'motif_center_start'] + dt[, 'motif_center_end'])/2
    dt$start <- center - left
    dt$end <- center + right
    
    # check if any is greater than the chromosome size
    #dt_split <- split(dt, f = dt[, 'chr'])
    dt_list <- lapply(genome_size$chr, function(each_chr){
        which_chr <- dt[dt$chr == each_chr, ]
        chr_size <- genome_size$size[genome_size$chr == each_chr]
        which_greater <- which(which_chr$end > chr_size)
        if(length(which_greater) > 0){
            which_chr[which_greater, ]$end <- chr_size
        }
        
        return(which_chr)
    })
    
    dt <- do.call(rbind, dt_list)
    
    return(dt)
}

# extend both the TP and TN by 
a <- expand_region(tp, gsize, both = 393216)
b <- expand_region(tn, gsize, both = 393216)

tp_dt <- cbind(a, motif_name=paste0('TP', 1:nrow(a)))
tn_dt <- cbind(b, motif_name=paste0('TN', 1:nrow(b)))

dt <- rbind(tp_dt, tn_dt)
#dt <- cbind(dt, pos=1:nrow(dt))

write.table(x=dt[!dt$chr %in% 'chrY', ], file=glue('{output_dir}/{TF}_motif_regions.txt'), quote = F, row.names = F, col.names = F)


# set.seed(2022)

# # just so I can work with a smaller dataset ===
# # also, I don't want chrY cos I don't have vcf files for it
# f1 <- tn_dt[sample(x=nrow(tn_dt), size = 5, replace = F), ]
# f1 <- cbind(f1, pos=row.names(f1))
# write.table(x=rbind(dt[1:5, ], f1), file=glue('{OUTPUT_DIR}/{TF}_motif_regions_TEMP.txt'), quote = F, row.names = F, col.names = F)

