suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option("--individuals", help='a list of loci'),
    make_option("--gene_expression", help='rds file of the splits'),
    make_option("--predictions_directory", help='rds file of the splits'),
    make_option("--output_file", help='.rds file to be created as the model')
)

opt <- parse_args(OptionParser(option_list=option_list))

print(opt)

library(glue)
library(data.table)
library(tidyverse)
library(abind)


dt_ids <- data.table::fread(opt$individuals, header = F)
mat_gexpr <- data.table::fread(opt$gene_expression) %>%
    tibble::column_to_rownames('gene') %>% as.matrix()

print(dim(mat_gexpr))
print(length(dt_ids$V1))

predicted_binding <- purrr::map(dt_ids$V1, function(idd){
    ff <- file.path(opt$predictions_directory, glue("{idd}.enpact_scores.tsv.gz"))
    if(file.exists(ff)){
        mat_tf <- data.table::fread(ff) %>% tibble::column_to_rownames('locus') %>% as.matrix()
        # select from gexpr
        mat_gexpr.id <- mat_gexpr[, idd, drop = F]
        # find common genes
        common_tss <- intersect(rownames(mat_gexpr.id), rownames(mat_tf))

        # arrange
        #mat.GEXPR <- mat_gexpr.id[common_tss, , drop = F]
        mat.TFB <- mat_tf[common_tss, ]
        #print(dim(mat.GEXPR)); print(dim(mat.TFB))
        return(mat.TFB)
    } else {
        return(NULL)
    }  
}, .progress = T) 

names(predicted_binding) <- dt_ids$V1
predicted_binding <- Filter(Negate(is.null), predicted_binding)

print(length(predicted_binding))

# combine into an array
myarray <- abind::abind(predicted_binding, along=3)
dimnames(myarray)[[3]] <- names(predicted_binding)
reshapedarray <- aperm(myarray, c(1, 3, 2), resize=TRUE)
lapply(dimnames(reshapedarray), length)
reshapedarray[1:3, 1:3, 1:3]
saveRDS(reshapedarray, file = opt$output_file, compress = "gzip")