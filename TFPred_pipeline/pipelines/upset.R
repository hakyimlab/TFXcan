

library(ComplexUpset)

dataset <- 'cistrome' # change this variable name 
TF <- 'FOXA1' #HOXB13' #'PPARG'
todays_date <- data_date <- '2023-01-11'

regions_dir <- glue('{project_dir}/files/defined_regions/{dataset}_{TF}')
save_dir <- glue('{regions_dir}/regions_data_{todays_date}')
binding_rdata <- readRDS(glue('{save_dir}/regions_information.RData'))
cistrome_dt <- binding_rdata$binding_matrix
cistrome_dt <- cistrome_dt[cistrome_dt$binding_class != 0, ]

cdt <- cistrome_dt[, -c(1:5)] |> as.data.frame()# |> t()
cdt <- cdt == 1
rownames(cdt) <- 1:nrow(cdt)

cdt <- cdt |> as.data.frame()

# change class to TF
colnames(cdt) <- paste0('TF ChIP ', 1:length(colnames(cdt)), sep='')

upset(cdt, colnames(cdt)[1:5], name='genre', min_size=0, width_ratio=0.1, wrap=TRUE, keep_empty_groups=F, mode='exclusive_intersection') +
ggtitle('Binding strength definition for 5 ChIP-seq experiments') +
theme(plot.title = element_text(size = 15, color='black'))

upset(cdt, colnames(cdt)[1:20], name='genre', max_size=500, width_ratio=0.1, wrap=TRUE, keep_empty_groups=F, sort_intersections_by=c('degree', 'cardinality'))


# ==== Freedman
dataset <- 'freedman' # change this variable name 
TF <- 'FOXA1' #HOXB13' #'PPARG'
todays_date <- data_date <- '2023-01-11'

regions_dir <- glue('{project_dir}/files/defined_regions/{dataset}_{TF}')
save_dir <- glue('{regions_dir}/regions_data_{todays_date}')
binding_rdata <- readRDS(glue('{save_dir}/regions_information.RData'))
freedman_dt <- binding_rdata$binding_matrix
freedman_dt <- freedman_dt[freedman_dt$binding_class != 0, ]

individuals <- data.table::fread(glue('{project_dir}/metadata/{dataset}_individuals.txt'), header=F)$V1
individuals <- individuals[-c(1)][1:5]


cdt <- freedman_dt[, individuals] |> as.data.frame()# |> t()
cdt <- cdt == 1
rownames(cdt) <- 1:nrow(cdt)
cdt <- cdt |> as.data.frame()
upset(cdt, colnames(cdt), name='individuals', min_size=50, width_ratio=0.1, wrap=TRUE, keep_empty_groups=T)