#!/bin/bash

rscript=${1}
utility_functions=${2}
vcf_files_arr=${3}
IFS=,$'\n' read -d '' -r -a vcf_arr < ${vcf_files_arr}
for vfile in "${vcf_arr[@]}"; do

    pattern=`echo ${vfile} | rev | cut -d '/' -f 1 | cut -d '.' -f 3- | rev`

    echo "Pattern is =" ${pattern}

    plink_results_pattern=${4}/${pattern}
    formatted_files_folder=${5}

    if [ ! -f ${plink_results_pattern}.traw ] && [ ! -f ${plink_results_pattern}.bim ]; then
        printf '\n%s\n' "INFO - preparing ${pattern}: plink portion"
        plink2 \
            --vcf ${vfile} \
            --geno 0.01 \
            --mind 0.01 \
            --make-bed \
            --maf 0.05 \
            --hwe 1e-6 \
            --out ${plink_results_pattern}

        plink2 \
            --bfile ${plink_results_pattern} \
            --recode A-transpose \
            --out ${plink_results_pattern}
    fi

    if [ ! -f ${formatted_files_folder}/${pattern}.geno.txt.gz ] && [ ! -f ${formatted_files_folder}/{pattern}.snp_annot.txt.gz ]; then
        printf '\n%s\n' "INFO - preparing ${pattern}: R portion"
        ${rscript} ${utility_functions} --file_prefix ${plink_results_pattern} --output_prefix ${formatted_files_folder}/${pattern} --command "create_genotype_dosage_and_snp_annot_files"
    fi
done
# Rscript -e "source('/beagle3/haky/users/temi/projects/TFXcan/scripts/ss.R'); summ(4, 5)"

# function create_file_formats(){
# }

# export -f create_file_formats

# "$@"