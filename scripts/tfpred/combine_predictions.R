
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--command", help='A list of files to combine'),
    make_option("--files", help='A list of files to combine'),
    make_option("--combined_file", help='The final output file.')
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

combine_predictions <- function(files, output_file){
    library(glue)
    library(data.table)
    library(tidyverse)

    files <- (files %>% strsplit(split=','))[[1]]
    print(glue('INFO - Found {length(files)} files to combine.'))

    files_dt <- purrr::map(files, function(each_file){
        if(file.exists(each_file)){
            dt <- data.table::fread(each_file) #%>% dplyr::rename(locus=V1)
            return(dt)
        }
    })

    files_dt <- purrr::reduce(files_dt, dplyr::left_join, by=c('locus'))
    #files_dt <- files_dt %>% tibble::column_to_rownames('locus') %>% as.matrix()
    data.table::fwrite(files_dt, file=output_file, quote=FALSE, row.names=FALSE, compress = 'gzip')
}


main <- function(command){
    run <- switch(as.character(command),
        'combine_predictions' = combine_predictions(files=opt$files, output_file=opt$combined_file),
        stop('The command does not exist.')
    )
    return(run)
}

main(opt$command)