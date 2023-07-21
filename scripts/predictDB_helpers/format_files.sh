#!/bin/bash

vcf_files_arr=${1}
IFS=,$'\n' read -d '' -r -a vcf_arr < ${vcf_files_arr}
for vfile in "${vcf_arr[@]}"; do

    printf '\n%s\n' "INFO - preparing for ${vfile}: plink portion"

    pattern=`echo ${vfile} | rev | cut -d '/' -f 1 | cut -d '.' -f 3- | rev`

    echo "Pattern is =" ${pattern}

    plink_results_pattern=${2}/${pattern}
    formatted_files_folder=${3}

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

        printf '\n%s\n' "INFO - preparing ${pattern}: R portion"
        /home/temi/miniconda3/envs/r-env/bin/Rscript /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/scripts/predictDB_helpers/format_files.R ${plink_results_pattern} ${formatted_files_folder}/${pattern}
    else
        printf '\n%s\n' "INFO - for ${pattern}, .traw and .bim files exist."
        printf '\n%s\n' "INFO - preparing ${pattern}: R portion"
        /home/temi/miniconda3/envs/r-env/bin/Rscript /lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/scripts/predictDB_helpers/format_files.R ${plink_results_pattern} ${formatted_files_folder}/${pattern}
    fi
done

# function create_file_formats(){
# }

# export -f create_file_formats

# "$@"