#!/bin/bash


function predict_cwas_weights () {

    weights_db=${1}
    #vcf_file_pattern=${2} # full path to the vcf file
     # full path to the weights db
    predixcan_dir=${2} # full path to the predixcan dir
    predict_exe=/lus/grand/projects/TFXcan/imlab/users/temi/software/MetaXcan/software/Predict.py
    vcf_file_pattern='/lus/grand/projects/TFXcan/imlab/data/baca_cwas/imputed_vcf_hg38_snps_only/chr*.dose.vcf.gz'

    m_name=$( echo ${weights_db} | rev | cut -d '_' -f 1 | rev )
    m_name=${m_name%.*}

    printf "\nINFO - running ${m_name}\n"

    # echo weights = $weights_db
    # echo vcf pattern = $vcf_file_pattern
    # echo predixcan dir = $predixcan_dir

    ${predict_exe} \
    --model_db_path ${weights_db} \
    --vcf_genotypes ${vcf_file_pattern} \
    --vcf_mode genotyped \
    --on_the_fly_mapping METADATA "{}_{}_{}_{}_b38" \
    --variant_mapping "${predixcan_dir}/tutorial/data/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz" id rsid \
    --prediction_output "${predixcan_dir}/output/cwas_individuals_imputed/${m_name}/baca_cwas_predict.txt" \
    --prediction_summary_output "${predixcan_dir}/output/cwas_individuals_imputed/${m_name}/baca_cwas_summary.txt" \
    --verbosity 9 \
    --throw

    #--text_sample_ids "/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/TFPred_pipeline/metadata/predixcan_geuvadis_240.txt"

    printf "INFO - finished ${m_name}"
}

export -f predict_cwas_weights

"$@"