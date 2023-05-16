#!/bin/bash

function split_vcf_per_chromosome (){

    j=${1}
    input_vcf=${2}
    output_folder=${3}

    echo "INFO - splitting off chr${j}"

    bcftools filter "${input_vcf}" -r "chr${j}" --output-type z --output "${output_folder}/chr${j}_CWAS_SNPs_2023-05-15_phased.vcf.gz" && tabix -p vcf "${output_folder}/chr${j}_CWAS_SNPs_2023-05-15_phased.vcf.gz"

    echo "INFO - Finished with chr${j}"
}

function rename_chromosomes () {

    from_vcf=${1}
    to_vcf=${2}
    chromosome_names_txt=${3}

    printf "\nRENAMING - ${from_vcf}\n"
    bcftools annotate --rename-chrs ${chromosome_names_txt} --threads 10 -Oz -o ${to_vcf} ${from_vcf}
    printf "INDEXING...\n"
    tabix -p vcf ${to_vcf}
    printf "DONE\n" && exit 0;
}

export -f split_vcf_per_chromosome
export -f rename_chromosomes

"$@"