# Author: Temi
# Date: Wednesday June 18 2025
# Description: NOT part of the current active pipeline — an earlier single-run (non-batched)
#   flashier factorization script, superseded by repeat_flash.R (see the commented-out block at
#   the bottom of repeat_flash.sbatch for the old call pattern). Kept for reference only.
# Usage: Rscript runFlashier.R --data <matrix> --column_for_rownames <col> --output_basename <path> --priorL <ebnm_...> --priorF <ebnm_...> [--transpose] --greedy_Kmax <int>

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--data", help='data matrix to factorize (.txt/.csv read with fread, or already-formatted matrix); see column_for_rownames/transpose below'),
    make_option("--column_for_rownames", help='name of the column in --data to use as row names when converting it to a matrix'),
    make_option("--output_basename", help='output path prefix; result written to {output_basename}.{priorL}-{priorF}.rds.gz'),
    make_option("--priorL", help='flashier EBNM prior function name for the loadings (L) matrix, e.g. ebnm_point_normal (see choosePrior() below)', default="ebnm_point_normal"),
    make_option("--priorF", help='flashier EBNM prior function name for the factors (F) matrix, e.g. ebnm_point_exponential (see choosePrior() below)', default="ebnm_point_exponential"),
    make_option("--transpose", type="logical", action='store_true', help='Should the data matrix be transposed after loading?'),
    make_option("--greedy_Kmax", type="integer", default=100L, help='max number of factors flashier greedily adds before backfitting (flashier::flash greedy_Kmax arg)')
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
