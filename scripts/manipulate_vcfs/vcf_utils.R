suppressPackageStartupMessages(library("optparse"))

option_list <- list(
    make_option(c("--command"), help='The command to run.'),
    make_option("--file_prefix", default=NULL, help = "Prefix of files before the extension."),
    make_option("--output_prefix", default=NULL, help="Full path plus prefix of output file."),
    make_option("--files_folder", default=NULL, help='Full path to output folder'),
    make_option("--output_folder", default=NULL, help='Full path to output folder'),
    make_option("--combine", default='no', help='Should you combine the SNP annotation and genotype files for all chromosomes.')
)

opt <- parse_args(OptionParser(option_list=option_list))

# some of these are meant to be used by some functions and not others
# file_prefix <- args[1]
# output_prefix <- args[2]

create_genotype_dosage_and_snp_annot_files <- function(file_prefix, output_prefix){

    library(glue)
    library(tidyverse)
    library(data.table)

    snp_annot <- data.table::fread(glue::glue("{file_prefix}.bim")) %>% 
        setnames(.,names(.), c("chr", "snp", "CM", "pos", "alt_vcf", "ref_vcf")) %>%
        dplyr::mutate(rsid = paste(chr, pos, ref_vcf, alt_vcf, 'b38', sep=':')) %>%
        dplyr::mutate(maf = 0.01) %>%  
        dplyr::mutate(varID = str_replace_all(rsid,":","_")) %>%
        dplyr::select(chr, pos, varID, ref_vcf, alt_vcf, maf, rsid)
    data.table::fwrite(snp_annot, file=glue("{output_prefix}.snp_annot.txt.gz"), sep='\t', quote=F, row.names=F, compress='gzip')

    genotype <- data.table::fread(glue::glue("{file_prefix}.traw")) %>% 
        tidyr::unite('varID', CHR, POS, COUNTED, ALT, sep = '_', remove=FALSE) %>%
        dplyr::mutate(varID = paste(varID, 'b38', sep='_')) %>%
        dplyr::select(-c(CHR,`(C)M`,POS, COUNTED, ALT, SNP)) %>% 
        setnames(.,names(.),gsub("0_", "", colnames(.))) %>% 
        dplyr::mutate(across(-varID, ~ case_when(
                . == 0 ~ 2, 
                . == 2 ~ 0,
                TRUE ~ .
        )))
    # genotype[genotype == 0] <- 2
    # genotype[genotype == 2] <- 0
    data.table::fwrite(genotype, file=glue("{output_prefix}.geno.txt.gz"), sep='\t', quote=F, row.names=F, compress='gzip')
}

merge_genotype_dosage_and_snp_annot_files <- function(files_folder, output_folder, combine){

    library(glue)
    library(tidyverse)
    library(data.table)

    valid_chr <- paste0('chr', c(1:22), sep='')
    geno_files <- file.path(files_folder, glue('*{valid_chr}*.geno.txt.gz')) |> Sys.glob() |> unique()
    snpannot_files <- file.path(files_folder, glue('*{valid_chr}*.snp_annot.txt.gz')) |> Sys.glob() |> unique()
    print(glue('INFO - Found {length(geno_files)} genotype .txt.gz files to merge.')) 
    print(glue('INFO - Found {length(snpannot_files)} snp annotation .txt.gz files to merge.'))

    #print(geno_files) ; print(snpannot_files)

    if(length(geno_files) == 0 | length(snpannot_files) == 0){
        stop(glue('ERROR - Either there are no geno .txt.gz files snp_annot .txt.gz files or both in {files_folder}.'))
    } else {

        all_geno_files <- purrr::map(geno_files, function(gfile){
            #ff <- glue('{files_folder}/ALL.{chrom}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.geno.txt')
            if(file.exists(gfile)){
                dt <- data.table::fread(gfile) %>% dplyr::mutate(varID=as.character(varID))
                return(dt)
            }
        }, .progress=TRUE)
        geno_file <- dplyr::bind_rows(all_geno_files)
        data.table::fwrite(geno_file, file=glue("{output_folder}/all_chrs.geno.txt"), sep='\t', quote=F, row.names=F)
        data.table::fwrite(geno_file, file=glue("{output_folder}/all_chrs.geno.txt.gz"), sep='\t', quote=F, row.names=F, compress='gzip')

        rm('all_geno_files')

        all_snpannot_files <- purrr::map(snpannot_files, function(sfile){
            #ff <- glue('{files_folder}/ALL.{chrom}.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.snp_annot.txt')
            if(file.exists(sfile)){
                dt <- data.table::fread(sfile) %>% dplyr::mutate(chr=as.character(chr))
                return(dt)
            }
        }, .progress=TRUE)
        snpannot_file <- dplyr::bind_rows(all_snpannot_files)
        data.table::fwrite(snpannot_file, file=glue("{output_folder}/all_chrs.snp_annot.txt"), sep='\t', quote=F, row.names=F)
        data.table::fwrite(snpannot_file, file=glue("{output_folder}/all_chrs.snp_annot.txt.gz"), sep='\t', quote=F, row.names=F, compress='gzip')
        
        rm('all_snpannot_files')

        if(combine == 'yes'){
            metadata='all_chrs'
            print("INFO - combining snp and geno files into one.")
            merged_dt <- snpannot_file %>% 
                dplyr::select(chr, varID, pos, ref_vcf, alt_vcf, maf) %>% 
                dplyr::inner_join(geno_file, by=c('varID' = 'varID'))

            ff1 <- glue("{output_folder}/{metadata}.text_dosages.txt.gz")
            data.table::fwrite(merged_dt, file=ff1, sep='\t', quote=F, row.names=F, col.names=F, compress = 'gzip')

            print("INFO - writing out samples file.")
            data.frame(colnames(geno_file)[-1], colnames(geno_file)[-1]) %>%
                data.table::fwrite(file=glue("{output_folder}/samples.text_dosages.txt"), sep='\t', quote=F, row.names=F, col.names=F)

        }
    }
}


main <- function(command){
    run <- switch(as.character(command),
        'create_genotype_dosage_and_snp_annot_files' = create_genotype_dosage_and_snp_annot_files(opt$file_prefix, opt$output_prefix),
        'merge_genotype_dosage_and_snp_annot_files' = merge_genotype_dosage_and_snp_annot_files(files_folder = opt$files_folder, output_folder = opt$output_folder, combine = opt$combine),
        stop('The command does not exist.')
    )
    return(run)
}

main(opt$command)