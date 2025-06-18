 

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--data", help='data preferrably in .rds format of a matrix of GWAS loci by TF/tissue paris of ratios of z-scores'),
    make_option("--splits", help='rds file of the splits'),
    make_option("--batch", help='batch name'),
    make_option("--output_basename", help='.rds file to be created as the model'),
    make_option("--priorL", help='alpha value for enet', default="ebnm_point_exponential"),
    make_option("--priorF", help='number of cores to use', default="ebnm_point_exponential"),
    make_option("--transpose", type="logical", action='store_true', help='Should the data be transposed?'),
    make_option("--greedy_Kmax", type="integer", default=100L, help='How many greedy iterations?'),
    make_option("--subset", type="numeric", help='Subset of loci to use', default=0.8),
    make_option("--niterations", type="numeric", default = 500L, help='Number of iterations')
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

library(glue)
library(data.table)
library(flashier)
library(magrittr)

what_extension <- tools::file_ext(opt$data)
if(what_extension %in% c('rds', 'rds.gz')){
    print(glue('INFO - Reading data from {opt$data}'))
    dt <- readRDS(opt$data)
} 

what_extension <- tools::file_ext(opt$split)
if(what_extension %in% c('rds', 'rds.gz')){
    sp <- readRDS(opt$splits)[[1]]
    print(glue('INFO - {opt$batch} has {length(sp)} splits'))
} 

# else if(what_extension == 'txt'){
#     print(glue('INFO - Reading data...'))
#     dt <- data.table::fread(opt$data) %>%
#         tibble::column_to_rownames(opt$column_for_rownames) %>%
#         as.matrix()
#     if(opt$transpose){
#         dt <- t(dt)
#     }
# }

print(glue('INFO - Input data has {nrow(dt)} loci and {ncol(dt)} samples'))

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

# create a list of random loci

# random_loci <- purrr::map(1:opt$niterations, function(k){
#     jj <- sample(nrow(dt), round(nrow(dt) * 0.8))
#     return(jj)
# })

flash_result <- purrr::map(1:length(sp), function(k){

    picked_loci <- sp[[k]]
    xdt <- dt[picked_loci, ]

    print(glue('INFO - Running {opt$batch}, iteration {k}'))
    print(glue('INFO - Randomly selected {nrow(xdt)} loci'))
    out <- flashier::flash(xdt, ebnm_fn = priors, backfit = TRUE, greedy_Kmax = opt$greedy_Kmax, verbose = 0)
    print(glue('INFO - Ran {opt$batch}, iteration {k}'))
    print(glue('INFO - Randomly selected {nrow(xdt)} loci'))
    print(glue('INFO - PVE {opt$batch}, iteration {k} is {sum(out$pve)}'))

    return(out)
}, .progress = TRUE)

# print(glue('INFO - Total explained variation using {opt$priorL}-{opt$priorF} is {sum(flash_result$pve)}'))
# write out result
rdsfile <- glue("{opt$output_basename}.{opt$batch}.{opt$priorL}-{opt$priorF}.Iters.rds.gz")
# check if folder exists or create
if(!dir.exists(dirname(rdsfile))){
    dir.create(dirname(rdsfile), recursive = TRUE)
}

saveRDS(flash_result, file = rdsfile, compress = "gzip")
print(glue('INFO - Flashier result saved to {rdsfile}'))


# flash_result <- list()

# for(k in 1:opt$niterations){
#     # randomly subset the data
#     print(glue('INFO - Running iteration {k}'))
#     xdt <- dt[sample(nrow(dt), round(nrow(dt) * 0.8)), ]
#     print(glue('INFO - Randomly selected {nrow(xdt)} loci'))
#     flash_result[[k]] <- flashier::flash(xdt, ebnm_fn = priors, backfit = TRUE, greedy_Kmax = opt$greedy_Kmax, verbose = 0)
#     print(glue('INFO - Ran iteration {k}'))
#     print(glue('INFO - Randomly selected {nrow(xdt)} loci'))
#     print(glue('INFO - PVE iteration {k} is {sum(flash_result[[k]]$pve)}'))
# }
