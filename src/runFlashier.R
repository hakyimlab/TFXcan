

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--data", help='data to train with enet'),
    make_option("--column_for_rownames", help=''),
    make_option("--output_basename", help='.rds file to be created as the model'),
    make_option("--priorL", help='alpha value for enet', default="ebnm_point_normal"),
    make_option("--priorF", help='number of cores to use', default="ebnm_point_exponential"),
    make_option("--transpose", type="logical", action='store_true', help='Should the data be transposed?'),
    make_option("--greedy_Kmax", type="integer", default=100L, help='How many greedy iterations?')
)

opt <- parse_args(OptionParser(option_list=option_list))

library(glue)
library(data.table)
library(flashier)
library(magrittr)

# load the data
if(file.exists(opt$data)){
    print(glue('INFO - Reading data...'))
    dt <- data.table::fread(opt$data) %>%
        tibble::column_to_rownames(opt$column_for_rownames) %>%
        as.matrix()
    if(opt$transpose){
        dt <- t(dt)
    }
} else {
    stop(glue('ERROR - Data cannot be found.'))
}

choosePrior <- function(inputPrior){
    out <- switch(as.character(inputPrior),
        'ebnm_point_normal' = ebnm_point_normal,
        'ebnm_point_exponential' = ebnm_point_exponential,
        'ebnm_point_laplace' = ebnm_point_laplace,
        'ebnm_normal' = ebnm_normal,
        'ebnm_normal_scale_mixture' = ebnm_normal_scale_mixture,
        stop(glue('ERROR - Prior {inputPrior} not recognized.'))
    )
    return(out)
}

priors <- c(priorL = choosePrior(opt$priorL), priorF = choosePrior(opt$priorF))
flash_result <- flashier::flash(dt, ebnm_fn = priors, backfit = TRUE, greedy_Kmax = opt$greedy_Kmax, verbose = 0)
print(glue('INFO - Total explained variation using {opt$priorL}-{opt$priorF} is {sum(flash_result$pve)}'))
# write out result
rdsfile <- glue("{opt$output_basename}.{opt$priorL}-{opt$priorF}.rds.gz")
saveRDS(flash_result, file = rdsfile, compress = "gzip")
print(glue('INFO - Flashier result saved to {rdsfile}'))
