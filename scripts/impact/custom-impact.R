

# This Generic pipeline works with split intersected Chip's data
# Author: Temi
# Date: Sep 13 2022
# 

annotate_with_impact <- function(n, output_dir='', common_dir=''){

    print(glue('Starting {n} ===== \n'))

    chip_peak_file <- chip_peak_file_split[[n]]

    motif <- paste0(TF, '_', cell_line)
    TF_BED <- unique(cbind(chip_peak_file[, c(1:3)], 0, motif, "+"))
    TF_BED[, 4] <- paste0("peak", seq(1,nrow(TF_BED), 1))

    colnames(TF_BED)[1:6] <- c('chr','start','end','id','score','strand') #score may be number of reads
    bed_Granges <- with(TF_BED, GenomicRanges::GRanges(chr, IRanges(start,end), strand, score, id = id))
    bed_Granges <- GenomicRanges::reduce(bed_Granges)

    tf_bed <- cbind(as.character(seqnames(bed_Granges)), start(bed_Granges), end(bed_Granges), paste0("peak", seq(1,length(bed_Granges),1)),0, "+")
    write.table(tf_bed, glue("{output_dir}/common_{n}_{TF}_{cell_line}_train_df.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

    cat('Dimension of chip peaks file for', n, 'is', dim(tf_bed))

    #homer find instances of motifs
    if(!file.exists(glue('{output_dir}/common_{n}_{TF}_{cell_line}_findinstances.txt'))){
        cmd <- glue('perl {homer_dir}/bin/findMotifsGenome.pl {output_dir}/common_{n}_{TF}_{cell_line}_train_df.txt {homer_dir}/data/genomes/hg19 {project_dir}/processed-data/homer_output/ -size given -find {common_dir}/{TF}_{cell_line}.motif > {output_dir}/common_{n}_{TF}_{cell_line}_findinstances.txt')
        system(cmd)
    } else {
        print(glue('{output_dir}/common_{n}_{TF}_{cell_line}_findinstances.txt exists.\n'))
    }

    # STEP 2 : prepare the positive and negative sets
    print('Preparing positive sets ===\n')
    motif_instances <- read.table(glue('{output_dir}/common_{n}_{TF}_{cell_line}_findinstances.txt'), sep = "\t", header = T, stringsAsFactors = FALSE)

    chip_regions <- read.table(glue('{output_dir}/common_{n}_{TF}_{cell_line}_train_df.txt'), sep='\t', header=F, stringsAsFactors=F) #as.data.frame(tf_bed) # or read.table(glue('{output_dir}/{TF}_train_df.txt'), sep='\t', header=T, stringsAsFactors=F)
    motif_library <- read.table(glue('{common_dir}/{TF}_{cell_line}_motif_threshold_info.txt'), sep = "\t", header = F, stringsAsFactors = FALSE)

    colnames(chip_regions) <- c("seqnames","starts","ends","names","scores","strands")
    #select ChIP peaks with a matching motif under them.
    chip_regions_full <- chip_regions
    m <- match(chip_regions[,4], motif_instances[,1]) # what chips are in motifs data?
    chip_regions <- chip_regions_full[which(is.na(m) == F),]

    if(nrow(chip_regions) == 0){
        train_test_peaks <- data.frame(matrix(ncol = 3, nrow = 0))
        #provide column names
        colnames(train_test_peaks) <- c('var1', 'var2', 'var3')

        write.table(train_test_peaks, glue("{output_dir}/common_{n}_{TF}_{cell_line}_train_test_positive_bed_EMPTY.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

        write.table(train_test_peaks, glue("{output_dir}/common_{n}_{TF}_{cell_line}_train_test_negative_bed_EMPTY.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)

    } else {

        chip_regions_newpeaks <- matrix(0,nrow(chip_regions),6) #dummy var number of rows 

        if(!file.exists(glue("{output_dir}/common_{n}_{TF}_{cell_line}_train_test_positive_bed.txt"))){
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

            write.table(chip_regions_newpeaks, glue("{output_dir}/common_{n}_{TF}_{cell_line}_train_df_newpeaks.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
        
        } else {
            print(glue('Train df newpeaks file exists. Moving on...\n'))
        }

        
        chip_regions_newpeaks <- read.table(glue('{output_dir}/common_{n}_{TF}_{cell_line}_train_df_newpeaks.txt'), sep = "\t", header = T, stringsAsFactors = FALSE)

        chip_regions_newpeaks_realdf <- data.frame(seqnames=chip_regions_newpeaks[,1],
                                            starts=as.numeric(chip_regions_newpeaks[,2]),
                                            ends=as.numeric(chip_regions_newpeaks[,3]),
                                            names=chip_regions_newpeaks[,4],
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
        if (length(w) > 0){
            train_test_peaks[w,3] <- train_test_peaks[w,3] + 1
        }

        w1 <- which(train_test_peaks$chr == "chrM")
        w2 <- which(nchar(train_test_peaks$chr) > 5)
        if (length(c(w1,w2)) > 0){
            train_test_peaks <- train_test_peaks[-c(w1,w2),]
        }

        write.table(train_test_peaks, glue("{output_dir}/common_{n}_{TF}_{cell_line}_train_test_positive_bed.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)


        # negative sets
        print('Preparing negative sets ===\n')
        train_df_newpeaks <- read.table(glue("{output_dir}/common_{n}_{TF}_{cell_line}_train_df_newpeaks.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)

        MotifThreshold <- read.table(glue("{common_dir}/{TF}_{cell_line}_motif_threshold_info.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)

        motifs <- MotifThreshold[1,1]
        instscan <- read.table(glue("{common_dir}/{TF}_{cell_line}_scanMotifsGenomewide_sort_15000.txt"), sep = "\t", header = F, stringsAsFactors = FALSE)
        instscan <- instscan[order(instscan[,6], decreasing = T),]   
        train_df <- read.table(glue("{output_dir}/common_{n}_{TF}_{cell_line}_train_df.txt"), sep = "\t", header = F, stringsAsFactors = FALSE) #old version: ChIPbed1.txt
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
        # if(dim(inst_bed)[1] == 0){
        #     inst_bed <- rbind(inst_bed, c(0, 0, 0))
        # }
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

        int_1 <- GenomicRanges::intersect(inst_bed_Granges, get(motifmatrix[1,3]))
        int_2 <- GenomicRanges::intersect(inst_bed_Granges, Granges_obj_mc)

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

        write.table(noChIP, glue("{output_dir}/common_{n}_{TF}_{cell_line}_noChIP.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
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

        write.table(noChIP, glue("{output_dir}/common_{n}_{TF}_{cell_line}_train_test_negative_bed.txt"), sep = "\t", quote = F, row.names = FALSE, col.names = FALSE)
    }

    cat(glue('\nEnded {n} \n'))

}



# === 


#elementMetadata(bed_Granges) <- TF_BED[, mtdata_names]

#ty <- cbind(bed_Granges, TF_BED[, mtdata_names])


# #chip_regions_newpeaks <- matrix(0,50,6)
# chip_regions_split <- split(chip_regions, chip_regions$seqnames)

# n_row <- nrow(chip_regions)

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



# STEP 3 : prepare the negative sets



# # create enformer beds =====
# tp <- read.table(glue("{output_dir}/{TF}_{cell_line}_train_test_positive_bed.txt"))
# tn <- read.table(glue("{output_dir}/{TF}_{cell_line}_train_test_negative_bed.txt"))


# # compare with the peak information file

# #peak_match <- match(peak_info$V2, tp$V4)
# peak_info_details <- peak_info[peak_info$V2 %in% tp$V4, ]
# peak_distribution <- sapply(strsplit(peak_info_details$V1, split='_'), getElement, 2)
# peak_distribution <- table(peak_distribution) |> as.data.frame() 
# colnames(peak_distribution) <- c('count', 'freq')
# peak_distribution <- peak_distribution %>% dplyr::mutate(count=as.numeric(count)) %>% dplyr::arrange(count)

# bp <- barplot(peak_distribution$freq, names.arg=peak_distribution$count, 
#     ylim=c(0, max(peak_distribution$freq + 10000)), 
#     ylab='No. of overlapping true ChIP peaks', xlab='Number of files (individuals)')

# #bp <- axis(1, at=1:length(bedfiles_names), labels=F, tick=T)
# #text(x=bp, y=peak_distribution$freq, label=peak_distribution$freq, cex=0.8, pos=4, srt=90) 
# mtext('Number of ChIP-seq peaks overlapped across number of Kawakami files')





# gsize <- read.table(glue('{output_dir}/../../../hg19.sizes'))


# expand_region <- function(dt, genome_size, both=128, left=NULL, right=NULL, create_granges=F, save_as_bed=T){
#     colnames(dt) <- c('chr','motif_center_start','motif_center_end','id','score','strand')
#     colnames(genome_size) <- c('chr', 'size')
    
#     if(both){
#         left <- both/2
#         right <- both/2
#     }
    
#     center <- (dt[, 'motif_center_start'] + dt[, 'motif_center_end'])/2
#     dt$start <- center - left
#     dt$end <- center + right
    
#     # check if any is greater than the chromosome size
#     #dt_split <- split(dt, f = dt[, 'chr'])
#     dt_list <- lapply(genome_size$chr, function(each_chr){
#         which_chr <- dt[dt$chr == each_chr, ]
#         chr_size <- genome_size$size[genome_size$chr == each_chr]
#         which_greater <- which(which_chr$end > chr_size)
#         if(length(which_greater) > 0){
#             which_chr[which_greater, ]$end <- chr_size
#         }
        
#         return(which_chr)
#     })
    
#     dt <- do.call(rbind, dt_list)
    
#     return(dt)
# }

# # extend both the TP and TN by 
# a <- expand_region(tp, gsize, both = 393216)
# b <- expand_region(tn, gsize, both = 393216)

# tp_dt <- cbind(a, motif_name=paste0('TP', 1:nrow(a)))
# tn_dt <- cbind(b, motif_name=paste0('TN', 1:nrow(b)))

# dt <- rbind(tp_dt, tn_dt)
# #dt <- cbind(dt, pos=1:nrow(dt))

# write.table(x=dt[!dt$chr %in% 'chrY', ], file=glue('{output_dir}/{TF}_{cell_line}_motif_regions.txt'), quote = F, row.names = F, col.names = F)


# # set.seed(2022)

# # # just so I can work with a smaller dataset ===
# # # also, I don't want chrY cos I don't have vcf files for it
# # f1 <- tn_dt[sample(x=nrow(tn_dt), size = 5, replace = F), ]
# # f1 <- cbind(f1, pos=row.names(f1))
# # write.table(x=rbind(dt[1:5, ], f1), file=glue('{OUTPUT_DIR}/{TF}_motif_regions_TEMP.txt'), quote = F, row.names = F, col.names = F)

# mr <- read.table('~/imlab_folder/users/temi/projects/TFXcan/processed-data/GATA3/T47D/GATA3_T47D_motif_regions.txt')

# # select 1500 TP and 1500 TN
# set.seed(2022)
# a <- mr[startsWith(mr$V9, 'TP'), ]
# a <- a[sample(1:nrow(a), 1500, replace=F), ]
# b <- mr[startsWith(mr$V9, 'TN'), ]
# b <- b[sample(1:nrow(b), 1500, replace=F), ]

# c <- rbind(a, b)
# c <- c[sample(1:nrow(c), nrow(c), replace=F), ]

# # pick a transcription factor
# TF <- 'GATA3'
# cell_line <- 'T47D'

# # where is the chip-seq data, as well as other directories
# cistrome_dir <- '/projects/covid-ct/imlab/data/cistrome/compressed'
# homer_dir <- '~/miniconda3/share/homer'
# project_dir <- '../'
# output_dir <- glue('../processed-data/{TF}/{cell_line}')

# write.table(x=c, file=glue('{output_dir}/{TF}_{cell_line}_motif_regions_3000_subset.txt'), quote = F, row.names = F, col.names = F)
