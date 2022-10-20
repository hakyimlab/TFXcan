


tp <- read.table(glue("{OUTPUT_DIR}/{TF}_train_test_positive_bed.txt"))
tn <- read.table(glue("{OUTPUT_DIR}/{TF}_train_test_negative_bed.txt"))

gsize <- read.table(glue('{INPUT_CHIP}/../hg19.sizes'))
dt <- tp
genome_size <- gsize

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

write.table(x=dt, file=glue('{OUTPUT_DIR}/{TF}_motif_regions.txt'), quote = F, row.names = F, col.names = F)


write.table(x=rbind(dt[1:5, ], tail(dt, 5)), file=glue('{OUTPUT_DIR}/{TF}_motif_regions_TEMP.txt'), quote = F, row.names = F, col.names = F)
