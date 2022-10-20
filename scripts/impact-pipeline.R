

<<<<<<< HEAD
#rstudioapi::getActiveDocumentContext()$path |> dirname() |> setwd()

setwd('/projects/covid-ct/imlab/users/temi/projects/TFXcan/scripts/')
getwd()
=======
rstudioapi::getActiveDocumentContext()$path |> dirname() |> setwd()

setwd('~/imlab_folder/users/temi/projects/TFXcan/scripts/')

>>>>>>> 435f3a119b45722200d2bfa0887ac044811be129

rm(list=ls())

library(glue)
library(GenomicRanges)
<<<<<<< HEAD
library(reticulate)
=======
>>>>>>> 435f3a119b45722200d2bfa0887ac044811be129

# ARGUMENTS
make_motifs <- T
TF <- 'Gata3' 
cell_type <- 'Th2'

# directories
INPUT_CHIP <- '../data/input_chip_data'
OUTPUT_DIR <- '../processed-data'
REPRO_DIR <- '../../../projects-reproduce/IMPACT-repro/reproducing'
HOMER_DIR <- '../../../software/homer'

# check if files are available ==============
if(!dir.exists(INPUT_CHIP)){
    dir.create(path=INPUT_CHIP, recursive = T)
} else {
    print(glue('The directory ({INPUT_CHIP}) exists.'))
}

# download the files =======
if(!file.exists(glue('{INPUT_CHIP}/Input_{TF}_1.bed.gz'))){
    download.file(url=glue('https://github.com/immunogenomics/IMPACT/raw/master/Training/Input_{TF}_1.bed.gz'), 
                  destfile = glue('{INPUT_CHIP}/Input_{TF}_1.bed.gz'))
}
# if(!file.exists(glue('{INPUT_CHIP}/Input_{TF}_2.bed.gz'))){
#     download.file(url='https://github.com/immunogenomics/IMPACT/raw/master/Training/Input_{TF}_2.bed.gz', 
#                   destfile = glue('{INPUT_CHIP}/Input_{TF}_2.bed.gz'))
# }

# prepare the chip inputs = mostly gotten from the IMPACT code itself on github
#Input ChIP
chip_data <- do.call(rbind, chip_data)[, 1:3]

motif <- base::paste(TF, cell_type, sep = '_')
chip_bed <- unique(cbind(chip_data[,1:3], 0, motif, "+"))
chip_bed[,4] <- paste0("peak", seq(1, nrow(chip_bed), 1))

colnames(chip_bed) <- c('chr','start','end','id','score','strand') #score may be number of reads
bed_Granges <- with(chip_bed, GRanges(chr, IRanges(start,end), strand, score, id = id))
bed_Granges <- reduce(bed_Granges)
chip_bed <- cbind(as.character(seqnames(bed_Granges)), start(bed_Granges), end(bed_Granges), paste0("peak",seq(1,length(bed_Granges),1)),0, "+")
write.table(chip_bed, glue('{OUTPUT_DIR}/{TF}_train_df.txt'), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

# the motif
motif_pwm <- read.table(glue("{REPRO_DIR}/Motifs/Motif_{motif}.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
write.table(motif_pwm, glue('{OUTPUT_DIR}/{TF}_1_PWM.txt'), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

#Input Motifs
MotifThreshvals <- motif_pwm[which(motif_pwm[,4]=="PH"), 3]
MotifThreshnames <- motif_pwm[which(motif_pwm[,4]=="PH"), 2]
MotifThreshinfo <- cbind(MotifThreshnames, MotifThreshvals)

write.table(motif_pwm, glue("{OUTPUT_DIR}/{TF}_All_PWM.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
write.table(MotifThreshinfo, glue("{OUTPUT_DIR}/{TF}_MotifThresholdList.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

if(!file.exists(glue('{OUTPUT_DIR}/{TF}_ALL_PWM.motif'))){
    if(file.exists((glue('{OUTPUT_DIR}/{TF}_ALL_PWM.txt')))){
        system(glue('cp {OUTPUT_DIR}/{TF}_ALL_PWM.txt {OUTPUT_DIR}/{TF}_ALL_PWM.motif'))
    } else {
        stop('Position weighted matrix not in folder.')
    }
}

# find instances of the motifs by using HOMER
print('Finding instances...')
system(glue('perl {HOMER_DIR}/bin/findMotifsGenome.pl {OUTPUT_DIR}/{TF}_train_df.txt {HOMER_DIR}/data/genomes/hg19/ homer_output/ -size given -find {OUTPUT_DIR}/{TF}_ALL_PWM.motif > {OUTPUT_DIR}/{TF}_findinstances.txt'))

# PREPARE POSITVE SET ===
train_df <- read.table(glue("{OUTPUT_DIR}/{TF}_train_df.txt"), sep = "\t", header = F, stringsAsFactors = FALSE) # CHIP
motif_instances <- read.table(glue("{OUTPUT_DIR}/{TF}_findinstances.txt"), sep = "\t", header = T, stringsAsFactors = FALSE) # MOTIFS
motif_threshold <- read.table(glue("{OUTPUT_DIR}/{TF}_MotifThresholdList.txt"), sep = "\t", header = F, stringsAsFactors = FALSE) # threshold to determine a motif ???

colnames(train_df) <- c("seqnames","starts","ends","names","scores","strands")

#select ChIP peaks with a matching motif under them.
train_df_full <- train_df
m <- match(train_df[,4], motif_instances[,1]) 
train_df <- train_df_full[which(is.na(m) == F),] # get matrix of only peaks with motif instance, now recenter
#add window around peaks 
#windowsize <- 380
train_df_newpeaks <- matrix(0,nrow(train_df),6) # dummy var number of rows 

# Here, we define the true positives as well as the center of the motifs 
for (i in 1:nrow(train_df)){
    print(paste0("row ", i, " of ", nrow(train_df)))
    w <- which(motif_instances$PositionID == train_df[i,4]) #match peakID, guaranteed to have a match 
    #only choose 1 motif, take the one with the top score
    w1 <- w[match(max(motif_instances$MotifScore[w]), motif_instances$MotifScore[w])]
    if (motif_instances$Strand[w1] == "+"){
        centerpeak <- motif_instances$Offset[w1] + floor(0.5*(nchar(motif_instances$Sequence[w1]))) #this is relative to start 
    }else{
        centerpeak <- motif_instances$Offset[w1] - floor(0.5*(nchar(motif_instances$Sequence[w1])))
    }
    newstart <- train_df$starts[i] + centerpeak #previous version of impact -windowsize/2 #establish left bound of peaks
    newend <- newstart + 1 #previous version of impact +windowsize #establish right bound of peaks 
    train_df_newpeaks[i,1] <- as.character(train_df$seqnames[i])
    train_df_newpeaks[i,2] <- newstart
    train_df_newpeaks[i,3] <- newend
    train_df_newpeaks[i,4] <- as.character(train_df$names[i]) #preserve old information
    train_df_newpeaks[i,5] <- as.character(train_df$scores[i])
    train_df_newpeaks[i,6] <- as.character(train_df$strands[i])
}

write.table(train_df_newpeaks, glue("{OUTPUT_DIR}/{TF}_train_df_newpeaks.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

# cleaning up 
train_df_newpeaks_realdf <- data.frame(seqnames=train_df_newpeaks[,1],
                                        starts=as.numeric(train_df_newpeaks[,2]),
                                        ends=as.numeric(train_df_newpeaks[,3]),
                                        names=train_df_newpeaks[,5],
                                        scores=0,
                                        strands="+")
#shuffle regions

s <- sample(seq(1,nrow(train_df_newpeaks_realdf),1),nrow(train_df_newpeaks_realdf),replace = F)

train_df_newpeaks_realdf <- train_df_newpeaks_realdf[s,]
#put training and test together! 
#use partition function later to make it easier to sample many times! 
train_test_peaks <- unique(train_df_newpeaks_realdf)
colnames(train_test_peaks) <- c("seqnames","starts","ends","names","scores","strands")

#clean up 
if (length(which(train_test_peaks$starts < 0)) > 0){
    train_test_peaks$start[which(train_test_peaks$start < 0)] <- 0
}

t_ <- train_test_peaks[,3] - train_test_peaks[,2] # should all be 1 : this is the center of the peaks
w <- which(t_ %% 2 != 0) #odd
if (length(w) > 0){
    train_test_peaks[w,3] <- train_test_peaks[w,3] + 1
} # they should all now be 2 or even

# remove some chromosomes - M, e.t.c
w1 <- which(train_test_peaks$chr == "chrM")
w2 <- which(nchar(train_test_peaks$chr) > 5)
if (length(c(w1, w2)) > 0){
    train_test_peaks <- train_test_peaks[-c(w1, w2), ]
}

write.table(train_test_peaks, glue("{OUTPUT_DIR}/{TF}_train_test_positive_bed.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)


### ====== True NEGATIVEs ================
# these are either motifs alone or random sequences
TF2 <- 'Gata3'

train_df_newpeaks <- read.table(glue("{OUTPUT_DIR}/{TF2}_train_df_newpeaks.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)

motif_threshold <- read.table(glue("{OUTPUT_DIR}/{TF2}_MotifThresholdList.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
motifs <- motif_threshold[1,1]
instscan <- read.table(glue("{REPRO_DIR}/Motifs/scanMotifsgenomewide.{motifs}.15000.sort.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
instscan <- instscan[order(instscan[,6], decreasing = T), ]   
train_df <- read.table(glue("{OUTPUT_DIR}/{TF2}_train_df.txt"), sep = "\t", header = F, stringsAsFactors = FALSE) #old version: ChIPbed1.txt
colnames(train_df) <- c('chr','start','end','id','score','strand') #score may be number of reads
tf_Granges <- with(train_df, GRanges(chr, IRanges(start,end), strand, score, id = id))
assign(paste0(motifs,"_Granges"), tf_Granges)

allinstscan <- instscan
sapply_instscan <- sapply(1:nrow(allinstscan), function(x){
    strsplit(allinstscan[x,1], "-")[[1]][1]
})
allinstscan[,1] <- sapply_instscan
allinstscan[,6] <- allinstscan[,6]/motif_threshold[1,2] # divide or normalize

w <- which(allinstscan[,6] < 1) 
if(length(w)>0){allinstscan <- allinstscan[-w,]}

motifmatrix <- matrix(0,nrow(motif_threshold),3)
#col1: motif name with cell type (1 per row)
#col2: degenerate motif name (used in scanMotifsGenome.pl), get with table(instances_scan[,1])
#col3: Granges files
motifmatrix[,1] <- motifs 
motifmatrix[,2] <- motifs
motifmatrix[,3] <- paste0(motifs, "_Granges")

#windowsize <- 380 why?????
print(paste0("working on ",motifmatrix[1,1]))
w_mcp <- which(train_df_newpeaks[,5] == motifmatrix[1,2]) #motifcentered peaks !!! THIS LOOKS WRONG BUT LEMME GO ON ====
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
if (length(ab) > 0){noChIP <- inst_bed[-ab,] #good negative peaks
}else{
    noChIP <- inst_bed
}
write.table(noChIP, glue("{OUTPUT_DIR}/{TF}_noChIP.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
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

write.table(noChIP, glue("{OUTPUT_DIR}/{TF}_train_test_negative_bed.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

# create enformer beds
tp <- read.table(glue("{OUTPUT_DIR}/{TF}_train_test_positive_bed.txt"))
tn <- read.table(glue("{OUTPUT_DIR}/{TF}_train_test_negative_bed.txt"))

gsize <- read.table(glue('{INPUT_CHIP}/../hg19.sizes'))

expand_region <- function(dt, genome_size, both=128, left=NULL, right=NULL, create_granges=F, save_as_bed=T){
    colnames(dt) <- c('chr','start','end','id','score','strand')
    colnames(genome_size) <- c('chr', 'size')
    
    if(both){
        left <- both/2
        right <- both/2
    }
    
    
    center <- (dt[, 'start'] + dt[, 'end'])/2
    dt[, 'start'] <- center - left
    dt[, 'end'] <- center + right
    
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
tp_dt <- cbind(expand_region(tp, gsize, both = 393216), set='TP')
tn_dt <- cbind(expand_region(tn, gsize, both = 393216), set='TN')

dt <- rbind(tp_dt, tn_dt)
dt <- cbind(dt, pos=1:nrow(dt))

write.table(x=dt[!dt$chr %in% 'chrY', ], file=glue('{OUTPUT_DIR}/{TF}_motif_regions.txt'), quote = F, row.names = F, col.names = F)


set.seed(2022)

# just so I can work with a smaller dataset ===
# also, I don't want chrY cos I don't have vcf files for it
f1 <- tn_dt[sample(x=nrow(tn_dt), size = 5, replace = F), ]
f1 <- cbind(f1, pos=row.names(f1))
write.table(x=rbind(dt[1:5, ], f1), file=glue('{OUTPUT_DIR}/{TF}_motif_regions_TEMP.txt'), quote = F, row.names = F, col.names = F)

