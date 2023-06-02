

rstudioapi::getActiveDocumentContext()$path |> dirname() |> setwd()

require(glue)
require(rjson)
require(readr)
require(stringr)
require(data.table)
require(GenomicRanges)

BiocManager::install("ENCODExplorerData")
BiocManager::install("AnnotationHub")
BiocManager::install("GEOquery")

library(GEOquery)

IMPACT_REPRO_DIR <- '../../../projects-reproduce/IMPACT-repro/reproducing'
FEATURES_RDATA <- '/Volumes/im-lab/nas40t2/temi/projects-reproduce/IMPACT-repro/reproducing/Rdata/Features.Rdata'
DATA_DIR <- '../data'
# this file also contains the betas for each mark

download.file('https://github.com/immunogenomics/IMPACT/raw/11f0a08cd37aa684b4867d61edde4a569a847df3/Features/Supplementary_Tables_Amariuta.xlsx', 
              destfile = glue('{DATA_DIR}/features_file.xlsx'))
impact_features_metadata <- readxl::read_xlsx(path = glue('{DATA_DIR}/features_file.xlsx'))

impact_features_metadata$mark |> table() |> sort() |> as.data.frame()

mt_data <- fromJSON(file = '/Volumes/im-lab/nas40t2/temi/projects/encode-to-impact/data/one-enformer-data.json')

mt_data <- fromJSON(file = 'https://www.encodeproject.org/search/?type=Biosample&biosample_ontology.classification=tissue&format=json')

lapply(mt_data$`@graph`, function(each_mtdata){
    cbind(each_mtdata$accession, each_mtdata$biosample_ontology$term_name)
})


# these are ENFORMER'S input data
enformer_features_metadata <- data.table::fread('https://raw.githubusercontent.com/calico/basenji/master/manuscripts/cross2020/targets_human.txt')

# And the  accession numbers for the files are in the `identifier` column
enformer_features_metadata$identifier


# I need to filter for those features used in impact
impact_features_metadata[impact_features_metadata$mark == 'DNase', ]
list_impact_features <- split(impact_features_metadata, f = impact_features_metadata$mark)
list_impact_features$DNase

enformer_features_metadata <- enformer_features_metadata %>%
    tidyr::separate(col=description, into = c('mark', 'detail'), sep = ':', extra = 'merge')

enformer_features_metadata$mark |> table()

# first lets filter out for the DNase data in both Enformer and Impact
enformer_dnase <- base::subset(enformer_features_metadata, mark=='DNASE')
impact_dnase <- list_impact_features$DNase


# impact granges features
load(file = FEATURES_RDATA)


# You can use AnnotationHub to retrieve encode_df_lite.
# http://bioconductor.org/packages/release/data/annotation/manuals/ENCODExplorerData/man/ENCODExplorerData.pdf
library(AnnotationHub)
hub <- AnnotationHub()

# subsetting for ENCODE data only ==========
# these are the ENCODE builds avaialable. 
encode_match <- !(str_match(title_hub, pattern = regex('.*(ENCODE).*'))[, 1] |> is.na())# need to find all ENCODE ANNOTATIONS
hub[encode_match, ]

encode_meta_df <- subset(hub, title=='ENCODE File Metadata (Light, 2019-10-13 build)')
encode_meta_df <- encode_meta_df[[1]]

# filter for those in the enformer 
# match with `ENCFF...` for others, match with `GSM` and `C...`
enformer_encode_mtdata <- encode_meta_df %>%
    dplyr::filter(file_accession %in% enformer_features_metadata[startsWith(enformer_features_metadata$identifier, 'ENCFF'), ]$identifier)

# subsetting for NCBI geo data ==========
# https://medicine.uiowa.edu/humangenetics/sites/medicine.uiowa.edu.humangenetics/files/wysiwyg_uploads/Workshop%20Follow-up_Day%202.pdf
(str_match(title_hub, pattern = regex('.*(NCBI).*'))[, 1])

hub$dataprovider |> unique()
hub[hub$dataprovider == "NCBI", ]$description

getGEO(GEO = 'GSM1208811') |> Table()

gsm_identifier <- enformer_features_metadata[startsWith(enformer_features_metadata$identifier, 'GSM'), ]$identifier
enformer_geo_mtdata <- lapply(1:length(gsm_identifier), function(i){
    tmp <- getGEO(gsm_identifier[i])
    tmp <- do.call(cbind, tmp |> Meta())
    tmp <- tmp[, colnames(tmp) %in% c('geo_accession', 'library_selection', 'library_strategy', 'source_name_ch1')]
    #cbind(identifier=gsm_identifier[i], tmp)
    tmp
})

names(enformer_geo_mtdata) <- enformer_features_metadata[startsWith(enformer_features_metadata$identifier, 'GSM'), ]$identifier # maybe not necessary really
enformer_geo_mtdata <- do.call(rbind, enformer_geo_mtdata)


# =========
chip <- matrix(0,2,3)
foxp3_chip_input <- data.table::fread(file='../data/Input_Foxp3_1.bed.gz', header = F)
chip <- rbind(chip, foxp3_chip_input[, 1:3])
chip <- chip[-c(1,2),]
bed <- unique(cbind(chip[,1:3], 0,'motif',"+"))
bed[,4] <- paste0("peak",seq(1,nrow(bed),1))

colnames(bed) <- c('chr','start','end','id','score','strand') #score may be number of reads
bed_Granges <- with(bed, GRanges(chr, IRanges(start,end), strand, score, id = id))
bed_Granges <- reduce(bed_Granges)
bed <- cbind(as.character(seqnames(bed_Granges)), start(bed_Granges), end(bed_Granges), paste0("peak",seq(1,length(bed_Granges),1)),0, "+")

as.numeric(bed[, 3]) - as.numeric(bed[, 2])


# ==== how are they overlapped

# this contains the chip seq data
train_df <- read.table("../../../projects-reproduce/IMPACT-repro/reproducing/output-files/train_df.txt", sep = "\t", header = F, stringsAsFactors = FALSE)

# this contains the motifs
instances <- read.table("../../../projects-reproduce/IMPACT-repro/reproducing/output-files/TFs.findinstances.txt", sep = "\t", header = T, stringsAsFactors = FALSE)
MotifLibrary <- read.table("../../../projects-reproduce/IMPACT-repro/reproducing/output-files/MotifThresholdList.txt", sep = "\t", header = F, stringsAsFactors = FALSE)

colnames(train_df) <- c("seqnames","starts","ends","names","scores","strands")
#select ChIP peaks with a matching motif under them; actually, this finds the first occurences
train_df_full <- train_df
m <- match(train_df[,4], instances[,1])

train_df <- train_df_full[which(is.na(m) == F),] #get matrix of only peaks with motif instance, now recenter
#add window around peaks 
#windowsize <- 380
train_df_newpeaks <- matrix(0,nrow(train_df),6)

i <- 1


for (i in 1:nrow(train_df)){
    print(paste0("row ", i, " of ", nrow(train_df)))
    w <- which(instances$PositionID == train_df[i,4]) #match peakID, gauranteed to have a match 
    #only choose 1 motif, take the one with the top score
    w1 <- w[match(max(instances$MotifScore[w]), instances$MotifScore[w])]
    if (instances$Strand[w1] == "+"){
        centerpeak <- instances$Offset[w1] + floor(0.5*(nchar(instances$Sequence[w1]))) #this is relative to start 
    }else{
        centerpeak <- instances$Offset[w1] - floor(0.5*(nchar(instances$Sequence[w1])))
    }
    newstart <- train_df$starts[i]+centerpeak #previous version of impact -windowsize/2 #establish left bound of peaks
    newend <- newstart + 1 #previous version of impact +windowsize #establish right bound of peaks 
    train_df_newpeaks[i,1] <- as.character(train_df$seqnames[i])
    train_df_newpeaks[i,2] <- newstart
    train_df_newpeaks[i,3] <- newend
    train_df_newpeaks[i,4] <- as.character(train_df$names[i]) #preserve old information
    train_df_newpeaks[i,5] <- as.character(train_df$scores[i])
    train_df_newpeaks[i,6] <- as.character(train_df$strands[i])
}

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

if (length(which(train_test_peaks$starts < 0)) > 0){
    train_test_peaks$start[which(train_test_peaks$start < 0)] <- 0
}

t <- train_test_peaks[,3]-train_test_peaks[,2]
w <- which(t %% 2 != 0) #odd
if (length(w) > 0){train_test_peaks[w,3] <- train_test_peaks[w,3] + 1}

# removes chromosome M and other weird chromosomes
w1 <- which(train_test_peaks$chr == "chrM")
w2 <- which(nchar(train_test_peaks$chr) > 5)
if (length(c(w1,w2)) > 0){
    train_test_peaks <- train_test_peaks[-c(w1,w2),]
}
