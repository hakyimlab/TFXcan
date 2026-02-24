suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--ppi_db", help='STRING PPI tsv'),
    make_option("--query_db", help='Where to sample the nulls from? tsv'),
    make_option("--programs", help='TFXcan programs'),
    make_option("--outputfile", help='TFXcan programs'),
    make_option("--topN", help='', default = 15),
    make_option("--nReplicates", help='How many replicates?', default = 10000),
    make_option("--useseed", help='How many replicates?', default = 2025)
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

library(glue)
library(data.table)
library(tidyverse)


get_random_factors <- function(tf_database, nsize = 20){
    dt <- tf_database %>% dplyr::slice_sample(n = nsize)
    # create pairwise combinations
    pairwise_combinations <- t(combn(dt$Factor, 2)) %>% as.data.table() %>% setnames(c('partnerA', 'partnerB'))
    return(pairwise_combinations)
}

generate_n_interacting_factors <- function(tf_database, nsize, setseed, n_networks = 50){
    set.seed(setseed)
    purrr::map(1:n_networks, function(i){
        get_random_factors(tf_database, nsize)
    })
}


#tfxcan_directory <- '/beagle3/haky/users/temi/projects/TFXcan'
# These are the clusters
tf_programs <- readRDS(opt$programs)
ppi_database <- data.table::fread(opt$ppi_db)
query_database <- data.table::fread(opt$query_db)

# db_cistrome.mapped <- data.table::fread(file.path(tfxcan_directory, "experiments/validation/STRINGPPI.cistromeDB.TFXcan.stringIDsmapped.tsv"))

count_ppi <- function(x, y){
    mt <- ppi_database %>% 
        dplyr::filter(preferredName_A == x & preferredName_B == y | preferredName_A == y & preferredName_B == x)
    if(nrow(mt) > 0){
        return(list(1, mean(mt$score)))
    } else {return(list(0, 0))}
}


estimate_null_ppi_interactions <- function(program, nsize, useseed=opt$useseed, nreplicates = opt$nReplicates){
    number_of_interactions <- list()
    average_ppi_score <- list()
    dt_rand <- generate_n_interacting_factors(query_database, nsize, useseed, nreplicates)
    dt_interactions <- purrr::map(seq_along(dt_rand), function(i){
        test_dt <- dt_rand[[i]][, ]
        setDT(test_dt)
        out <- test_dt[, coefs := Map(count_ppi, partnerA, partnerB)]
        out <- out[, coefs[[1]], by=list(partnerA, partnerB)] %>% dplyr::rename(num_interactions = V1, avg_ppi_score = V2)
        number_of_interactions[[i]] <<- sum(out$num_interactions)
        average_ppi_score[[i]] <<- mean(out$avg_ppi_score)
        out
    }, .progress = F)

    null.number_of_interactions <- unlist(number_of_interactions)
    null.average_ppi_score <- unlist(average_ppi_score)

    data.table(program = program, null_n_interactions = null.number_of_interactions, null_avg_ppi_score = null.average_ppi_score)
}


consensus_programs <- lapply(seq_along(tf_programs), function(i){
    data.table::as.data.table(tf_programs[[i]], keep.rownames = "tf_tissue") %>%
        setnames(c('tf_tissue', glue::glue('consensus_contribution_{i}')))
}) %>% purrr::reduce(dplyr::inner_join, by = 'tf_tissue')


program_names <- names(tf_programs)

# opt <- list()

# opt$topN <- 10

top.n.TFs.programs <- purrr::map(1:(ncol(consensus_programs)-1), function(i){
    cpr <- glue("consensus_contribution_{i}")
    consensus_programs %>% dplyr::select(tf_tissue, all_of(cpr)) %>%
    dplyr::arrange(desc(.data[[cpr]])) %>% 
    dplyr::slice_head(n = opt$topN) %>% 
    tidyr::separate_wider_delim(col = tf_tissue, delim = '_', names = c('tf', 'tissue'), cols_remove = F) %>%
    dplyr::distinct(tf) %>% as.data.frame() #%>% dplyr::mutate(program = program_names[i])
}) #%>% dplyr::bind_rows() %>% dp lyr::relocate(program)
names(top.n.TFs.programs) <- program_names


full_null_dt <- purrr::map(names(top.n.TFs.programs), function(nxn){
    #outfile <- glue::glue("{tfxcan_directory}/experiments/validation/ppi/STRINGPPI.{nxn}.null_interactions.tsv")
    message(glue::glue("INFO - Simulating nulls for {nxn}..."))
    res <- tryCatch({
        dtd <- nrow(top.n.TFs.programs[[nxn]])
        ft <- estimate_null_ppi_interactions(nxn, dtd)
        return(ft)
    }, error = function(err){
        return(NULL)
    })

    if(!is.null(res)){
        jj <- res %>% data.table::fwrite(outfile, sep = '\t', quote = F, row.names = F, col.names = TRUE)
    }
    return(res)
    
}, .progress = F)

names(full_null_dt) <- names(top.n.TFs.programs)

saveRDS(full_null_dt, opt$outputfile)
