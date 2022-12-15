

test_models <- function(model, X_test_set, y_test_set){
    features <- model$glmnet.fit$beta |> rownames()
    X_test_set <- X_test_set[, features]
    assess.glmnet(model, newx = X_test_set, newy = y_test_set) |> unlist()
}

# given a glmnet model, collect the coeffiecients into a dataframe
collate_coefficients <- function(glmnetCVmodel, return_best_beta=T){

    dimensions <- glmnetCVmodel$glmnet.fit$beta@Dim
    coef_mat <- as.data.frame(summary(glmnetCVmodel$glmnet.fit$beta))
    
    temp_mat <- matrix(data=NA, nrow=dimensions[1], ncol=dimensions[2])
    for(i in 1:nrow(coef_mat)){
        temp_mat[coef_mat[i, 'i'], coef_mat[i, 'j']] <- coef_mat[i, 'x']
    }

    temp_mat[is.na(temp_mat)] <- 0

    if(return_best_beta){
        min_error_index <- glmnetCVmodel$index['min', ]
        one_sd_index <- glmnetCVmodel$index['1se', ]

        df_beta <- cbind(temp_mat[, min_error_index], temp_mat[, one_sd_index]) |> as.data.frame()
        colnames(df_beta) <- c('min', '1se')
        return(df_beta)
    }

    return(temp_mat)
}

expand_region <- function(dt, genome_size, both=128, left=NULL, right=NULL, create_granges=F, save_as_bed=T){

    #gsize <- read.table(glue('{INPUT_CHIP}/../hg19.sizes'))

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
