

current_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(current_path)

current_date <- Sys.Date()

download_path <- '/lus/theta-fs0/projects/covid-ct/imlab/data/freedman_data'
download_file <- 

if(!dir.exists('./files')){
    dir.create('./files')
} else {
    print('./files directory exists.')
}

file_url <- 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE161nnn/GSE161948/suppl/GSE161948_RAW.tar'

file_name <- paste0('./files/freedman_human_raw_bed_files_', current_date, '.tar')
download.file(url=file_url, destfile=file_name)