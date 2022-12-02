



test_granges <- files_granges[[1]]

kawakami_regions_granges <- kawakami_regions_granges |> unique()

# first overlaps between the individual's granges and with defined Kawakami regions to subset for those that are common to Kawakami
ind_kawakami_granges_overlap <- findOverlaps(query = test_granges, subject = kawakami_regions_granges, type = 'any')
overlaps_df <- data.frame(as.data.frame(test_granges)[queryHits(ind_kawakami_granges_overlap),], as.data.frame(kawakami_regions_granges)[subjectHits(ind_kawakami_granges_overlap),])

# convert these peaks to a granges object
freedman_df_granges <- with(overlaps_df, GenomicRanges::GRanges(seqnames, IRanges(start, end), strand))

# next overlap these subsets with the top 100000 motif regions
with_motif_granges_overlap <- findOverlaps(query = freedman_df_granges, subject = motif_cutoff_granges, type = 'any')

# the resulting subset are the truly bound regions
freedman_truly_bound_df <- data.frame(as.data.frame(freedman_df_granges)[queryHits(with_motif_granges_overlap),], motif_cutoff_granges[subjectHits(with_motif_granges_overlap),])
colnames(freedman_truly_bound_df)[c(1, 7, 8)] <- c('chr', 'start', 'end')
freedman_truly_bound_df <- freedman_truly_bound_df[, c(1, 7, 8)] 

# select distinct columns
freedman_truly_bound_df <- freedman_truly_bound_df %>% dplyr::distinct(chr, start, end)
#freedman_truly_bound_df |> head() ; freedman_truly_bound_df |> dim()

bound_with_kawakami_granges <- with(freedman_truly_bound_df, GenomicRanges::GRanges(chr, IRanges(start, end), strand='*'))
bound_with_kawakami_overlap <- findOverlaps(query = bound_with_kawakami_granges, subject = kawakami_regions_granges, type = 'any')

 # find the unbound regions (every other region shold be unbound really)
 # need to find overlap of the truly bound region with kawakami and then remove 
freedman_truly_unbound_df <- data.frame(kawakami_regions_granges[ - subjectHits(bound_with_kawakami_overlap),])
freedman_truly_unbound_df <- freedman_truly_unbound_df[c(1,2,3)]
colnames(freedman_truly_unbound_df) <- c('chr', 'start', 'end')
# select distinct columns
freedman_truly_unbound_df <- freedman_truly_unbound_df %>% dplyr::distinct(chr, start, end)

# select a random X3 of the truly bound regions
# set.seed(2022)
# freedman_truly_unbound_df <- freedman_truly_unbound_df[sample(1:nrow(freedman_truly_unbound_df), 3*nrow(freedman_truly_bound_df)), ]
#freedman_truly_unbound_df |> head() ; freedman_truly_unbound_df |> dim()

freedman_truly_unbound_df <- cbind(with(freedman_truly_unbound_df, paste(chr, start, end, sep='_')), 0)
freedman_truly_bound_df <- cbind(with(freedman_truly_bound_df, paste(chr, start, end, sep='_')), 1)

freedman_defined_regions <- rbind(freedman_truly_bound_df, freedman_truly_unbound_df)
#shuffle
freedman_defined_regions <- freedman_defined_regions[sample(1:nrow(freedman_defined_regions), nrow(freedman_defined_regions)), ]




 # first overlap with defined Kawakami regions to subset for those that are common to Kawakami
    granges_overlap <- findOverlaps(query = ind_granges, subject = kawakami_regions_granges, type = 'any')
    overlaps_df <- data.frame(as.data.frame(ind_granges)[queryHits(granges_overlap),], as.data.frame(kawakami_regions_granges)[subjectHits(granges_overlap),])

    # use the peaks
    freedman_df_granges <- with(overlaps_df, GenomicRanges::GRanges(seqnames, IRanges(start, end), strand))
    
    # next overlap these subsets with the top 10000 motif regions
    granges_overlap <- findOverlaps(query = freedman_df_granges, subject = motif_cutoff_granges, type = 'any')
    freedman_truly_bound_df <- data.frame(as.data.frame(freedman_df_granges)[queryHits(granges_overlap),], motif_cutoff_granges[subjectHits(granges_overlap),])
    colnames(freedman_truly_bound_df)[c(1, 7, 8)] <- c('chr', 'start', 'end')
    freedman_truly_bound_df <- freedman_truly_bound_df[, c(1, 7, 8)] 

    # select distinct columns
    freedman_truly_bound_df <- freedman_truly_bound_df %>% dplyr::distinct(chr, start, end)
    #freedman_truly_bound_df |> head() ; freedman_truly_bound_df |> dim()

    # find the unbound regions (every other region shold be unbound really)
    freedman_truly_unbound_df <- data.frame(motif_cutoff_granges[ - subjectHits(granges_overlap),])
    freedman_truly_unbound_df <- freedman_truly_unbound_df[c(1,2,3)]
    colnames(freedman_truly_unbound_df) <- c('chr', 'start', 'end')
    # select distinct columns
    freedman_truly_unbound_df <- freedman_truly_unbound_df %>% dplyr::distinct(chr, start, end)

    # select a random X3 of the truly bound regions
    set.seed(2022)
    freedman_truly_unbound_df <- freedman_truly_unbound_df[sample(1:nrow(freedman_truly_unbound_df), 3*nrow(freedman_truly_bound_df)), ]
    #freedman_truly_unbound_df |> head() ; freedman_truly_unbound_df |> dim()

    freedman_truly_unbound_df <- cbind(with(freedman_truly_unbound_df, paste(chr, start, end, sep='_')), 0)
    freedman_truly_bound_df <- cbind(with(freedman_truly_bound_df, paste(chr, start, end, sep='_')), 1)

    freedman_defined_regions <- rbind(freedman_truly_bound_df, freedman_truly_unbound_df)
    #shuffle
    freedman_defined_regions <- freedman_defined_regions[sample(1:nrow(freedman_defined_regions), nrow(freedman_defined_regions)), ]