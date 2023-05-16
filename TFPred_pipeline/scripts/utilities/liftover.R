# liftover is pretty quick

arguments <- commandArgs(trailingOnly=TRUE)

library(glue)
library(parallel)

liftover_cmd <- arguments[1]
chain_file <- arguments[2]
input_bed <- arguments[3]
output_bed <- arguments[4]

peaks_dir <- '/lus/grand/projects/covid-ct/imlab/data/freedman_data/peak_files'
chain_file <- '/lus/grand/projects/covid-ct/imlab/data/liftover_files/hg19ToHg38.over.chain.gz'

output_dir <- '/lus/grand/projects/covid-ct/imlab/data/freedman_data/peaks_liftover'
if(!dir.exists(output_dir)){dir.create(output_dir, recursive=T)}

unmapped_dir <- '/lus/grand/projects/covid-ct/imlab/data/freedman_data/unmapped_liftover'
if(!dir.exists(unmapped_dir)){dir.create(unmapped_dir, recursive=T)}

#files_to_liftover <- list.files(peaks_dir, full.names=T)

#fnames <- gsub('.bed', '', basename(files_to_liftover))

liftover_cmd <- '/lus/grand/projects/covid-ct/imlab/users/temi/software/liftOver'

print(glue('[INFO] Found {length(fnames)} files to liftover.'))

for(i in 1:length(fnames)){
    print(glue('[INFO] Lifting over: {i} - {fnames[i]}.bed'))
    system(glue('{liftover_cmd} {peaks_dir}/{fnames[i]}.bed {chain_file} {output_dir}/{fnames[i]}.bed {unmapped_dir}/{fnames[i]}.bed'))
}

# parallel::mclapply(1:length(fnames), function(i){
#     print(glue('[INFO] Lifting over: {i} - {fnames[i]}.bed'))
#     system(glue('{liftover_cmd} {peaks_dir}/{fnames[i]}.bed {chain_file} {output_dir}/{fnames[i]}.bed {unmapped_dir}/{fnames[i]}.bed'))
# })
