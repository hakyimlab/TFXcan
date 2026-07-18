# Author: Temi
# Description: shared R helper library — plot themes and small transforms sourced by other scripts/notebooks.
#   Not a runnable script (no optparse, nothing executes on its own).
# Usage: source("modules.R")

# convert logit-scale values (e.g. model scores) to probabilities via the inverse logit
logitsToProbabilities <- function(x){
    if(any(is.na(x))){
        stop('ERROR - There are NAs in your supplied x. Remove them first.')
    }
    out <- exp(x)/(1 + exp(x))
    return(out)
}

# minimal ggplot theme: no gridlines, small Helvetica text, boxed panel border
theme_simple <- theme_bw() + theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    plot.title = element_text(family = 'Helvetica', colour = 'black', size = 9.33),
    text=element_text(size=9.33, colour = 'black', family = 'Helvetica'),
    axis.text = element_text(hjust = 1), 
    axis.text.x = element_text( hjust=0.5), 
    panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)
)

# theme add-on: draws a small legend inside the plot panel (top-right-ish) instead of beside it
theme_inside_legend <- theme(legend.position = "inside",
        legend.position.inside = c(0.8, 0.8), 
        legend.background = element_rect(fill = 'transparent'),
        legend.title = element_text(size = 3),          # title text size
        legend.text = element_text(size = 3),           # label text size
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.key.size = unit(0.3, "cm"),              # overall swatch box size (both dims)
        legend.key.height = unit(0.3, "cm"),            # override height only, if needed
        legend.key.width = unit(0.3, "cm"),             # override width only, if needed
        legend.key.spacing = unit(0.1, "cm"),           # gap between multiple keys (ggplot2 ≥3.5)
        legend.spacing = unit(0.1, "cm"),               # gap between legend and other elements
        legend.margin = margin(0, 0, 0, 0))             # padding around the whole legend box

# QQ-plot of observed vs. expected -log10(p) under a uniform null, thresholded to avoid -Inf at tiny p-values
qqunif.ggplot <- function(p, mlog10_p_thres=30){
    require(ggplot2)
    ## thresholded by default at 1e-30
    p=na.omit(p)
    nn = length(p)
    xx =  -log10((1:nn)/(nn+1))
    
    p_thres = 10^{-mlog10_p_thres}
    if(sum(p < p_thres)){
      warning(paste("thresholding p to ", p_thres))
      p = pmax(p, p_thres)
    }

    dt <- data.table(pval = xx, dist = -sort(log10(p)))

    pp <- ggplot2::ggplot(dt, aes(x = pval, y = dist)) +
        geom_point() +
        labs(x = expression(Expected~~-log[10](italic(p))),
          y = expression(Observed~~-log[10](italic(p)))) +
          geom_abline(aes(intercept = 0, slope = 1), col='gray', lty = 1)
        
    return(pp)
}

# plot a precision-recall curve (from PRROC::pr.curve) with the AUPRC printed on the panel
plot_pr_curve.ggplot <- function(pr_object, tf_tissue){
    dt_pr_curve <- pr_object$curve %>% as.data.table() %>% setnames(c('recall', 'precision', 'unk'))

    tlabel <- glue("AUPRC: {round(pr_object$auc.integral, 3)}")
    ttitle <- glue("Precision-Recall curve for\nEnpact's {tf_tissue} Model")

    pp <- dt_pr_curve %>% ggplot(aes(x = recall, y = precision)) +
        ggplot2::geom_line() +
        ylim(0, 1) +
        xlim(0, 1) +
        labs(x = 'Recall', y = 'Precision', title = ttitle) +
        theme_simple +
        geom_text(label = tlabel, col = 'darkgreen', x = 0.2, y = 0.2, size = 2.0)

    return(pp)

    #
}

# overlay two QQ-plots (e.g. observed vs. null p-values) on the same axes for comparison
qqcompare.ggplot <- function(p1, p2, names = c('observed', 'null'), mlog10_p_thres=30){
    require(ggplot2)
    ## thresholded by default at 1e-30
    p1 = na.omit(p1)
    p2 = na.omit(p2)
    p_thres = 10^{-mlog10_p_thres}
    if(!is.null(p_thres)){pvec1 = pmax(p1, p_thres); pvec2 = pmax(p2, p_thres)}
    nn = length(pvec1)
    xx =  -log10((1:nn)/(nn+1))

    pp <- ggplot2::ggplot() +
        geom_point(aes(x = xx, y = -sort(log10(pvec1)), col = names[1])) +
        geom_point(aes(x = xx, y = -sort(log10(pvec2)), col = names[2])) +
        labs(x = expression(Expected~~-log[10](italic(p))),
          y = expression(Observed~~-log[10](italic(p)))) +
          geom_abline(aes(intercept = 0, slope = 1), col='gray', lty = 1)
        
    return(pp)
}


# alias for qqunif.ggplot() above -- same QQ-plot, kept as a separate name since both are
# called from different notebooks; delegates instead of duplicating the body
gg_qq_empirical <- function(p, mlog10_p_thres=30){
    qqunif.ggplot(p, mlog10_p_thres)
}

# generic QQ-plot of any two numeric vectors x vs y (not p-value specific)
ggplot_qqplot <- function(x, y){
    qqdata <- qqplot(y = y, x = x, plot.it = F)
    ggplot2::ggplot() +
        geom_point(aes(x = qqdata$x, y = qqdata$y))
}

# for one cluster, pick one representative column per K (from names like "..._K") out of mat — either a random
# sample (samp=TRUE) or all of them — and return the matching rows of mat, transposed
select_ks <- function(cluster_designation, cluster, mat, samp = TRUE){
    yx <- cluster_designation[cluster_designation == cluster]
    kx <- sapply(strsplit(names(yx), '_'), getElement, 2) |> unique()
    yx_kx <- lapply(kx, function(x) yx[endsWith(names(yx), x)])

    set.seed(2025)

    if(samp){
        selected_yx_kx <- sapply(yx_kx, function(x){
            x <- x[sample.int(length(x), 1)]
            return(x)
        }) %>% names()
    } else {
        selected_yx_kx <- names(unlist(yx_kx))
    }
    

    tryCatch({
       mt <- mat[selected_yx_kx, ] |> t()
       return(mt)
    }, error = function(e){
        message(glue('Error: {e} at {cluster}'))
    })
}

# reshape a wide tf_tissue-by-subsample contribution matrix to long format. average=FALSE keeps only the
# top-N contributing tf_tissue rows per subsample; average=TRUE instead averages contribution across
# subsamples per tf_tissue (topn is ignored in that branch)
create_cluster_data <- function(dt_yselected, topn = 20, average = FALSE){
    # convert to data.table and to long format

    if(average == FALSE){
        dt_yselected <- as.data.table(dt_yselected, keep.rownames = 'tf_tissue') %>% 
            tidyr::pivot_longer(cols = !tf_tissue, names_to = 'subsample', values_to = 'contribution') %>%
            tidyr::separate(tf_tissue, c('tf', 'tissue'), sep = '_', remove = F) %>%
            dplyr::group_by(subsample) %>%
            dplyr::arrange(desc(contribution)) %>%
            slice_head(n = topn) %>%
            dplyr::ungroup() %>%
            dplyr::group_by(tf_tissue) %>%
            dplyr::mutate(tf_tissue = factor(tf_tissue, levels = unique(.$tf_tissue)),
                tf_tissue = reorder(tf_tissue, -contribution))
    } else if(average){
        dt_yselected <- as.data.table(dt_yselected, keep.rownames = 'tf_tissue') %>% 
            tidyr::pivot_longer(cols = !tf_tissue, names_to = 'subsample', values_to = 'contribution') %>%
            tidyr::separate(tf_tissue, c('tf', 'tissue'), sep = '_', remove = F) %>%
            dplyr::group_by(tf_tissue) %>%
            dplyr::summarize(contribution = mean(contribution)) %>%
            dplyr::mutate(tf_tissue = factor(tf_tissue, levels = unique(.$tf_tissue)),
                tf_tissue = reorder(tf_tissue, -contribution))
    }

    return(dt_yselected)
}