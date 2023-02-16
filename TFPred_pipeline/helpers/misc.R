

library(data.table)
library(glue)

plogs <- data.table::fread('/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/helpers/geuvadis_files_sizes.out')

plogs_save <- plogs[plogs$V1 == 18655, ]$V2 %>% gsub('_log.csv', '', .)
plogs_delete <- plogs[plogs$V1 != 18655, ]$V2

plots_delete <- plogs_delete[-length(plogs_delete)]

glue("rm -rf {gsub('_log.csv', '', plots_delete)}")

paste('rm -rf ', gsub('_log.csv', '', plots_delete), collapse = NULL)

glue("rm -rf {plots_delete}")

# write psave to afile

utils::write.table(plogs_save, '/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/metadata/completed_geuvadis_240.txt', row.names=F, quote=F, col.names=F)

predixcan_ids <- data.frame(cbind(plogs_save, plogs_save))
utils::write.table(predixcan_ids, '/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/metadata/predixcan_geuvadis_240.txt', row.names=F, quote=F, col.names=F, sep='\t')
