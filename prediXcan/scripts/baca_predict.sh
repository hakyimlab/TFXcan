
temi_dir="/lus/grand/projects/covid-ct/imlab/users/temi"
predixcan_dir="/lus/grand/projects/covid-ct/imlab/users/temi/projects/prediXcan"
conda activate imlabtools

for model_db in `ls ${temi_dir}/projects/TFXcan/baca_cwas/db_folder/`; do
    m_name=$( echo ${model_db} | cut -d '_' -f 3 )
    m_name=${m_name%.*}
    echo ${m_name}

    ${temi_dir}/software/MetaXcan/software/Predict.py \
    --model_db_path "${temi_dir}/projects/TFXcan/baca_cwas/db_folder/baca_cwas_${m_name}.db" \
    --vcf_genotypes ${temi_dir}/projects/TFXcan/geuvadis_individuals/vcf_files/chr*_99_individuals.vcf.gz \
    --vcf_mode genotyped \
    --on_the_fly_mapping METADATA "{}_{}_{}_{}_b38" \
    --variant_mapping "${predixcan_dir}/tutorial/data/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz" id rsid \
    --prediction_output "${predixcan_dir}/output/${m_name}/baca_cwas_predict.txt" \
    --prediction_summary_output "${predixcan_dir}/output/${m_name}/baca_cwas_summary.txt" \
    --verbosity 9 \
    --throw 
done




























${temi_dir}/software/MetaXcan/software/Predict.py \
--model_db_path ${temi_dir}/projects/TFXcan/baca_cwas/db_folder/baca_cwas_top1.qtl.db \
--vcf_genotypes ${temi_dir}/projects/TFXcan/geuvadis_individuals/vcf_files/chr22_99_individuals.vcf.gz \
--vcf_mode genotyped \
--model_db_snp_key varID \
--on_the_fly_mapping METADATA "{}_{}_{}_{}_b38" \
--prediction_output ${predixcan_dir}/output/baca_cwas_top1.qtl_predict.txt \
--prediction_summary_output ${predixcan_dir}/output/baca_cwas_top1.qtl_summary.txt \
--verbosity 9 \
--throw


${temi_dir}/software/MetaXcan/software/Predict.py \
--model_db_path ${temi_dir}/projects/TFXcan/baca_cwas/db_folder/baca_cwas_lasso.db \
--vcf_genotypes ${temi_dir}/projects/TFXcan/geuvadis_individuals/vcf_files/chr*_99_individuals.vcf.gz \
--vcf_mode genotyped \
--on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" \
--variant_mapping ${predixcan_dir}/tutorial/data/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz id rsid \
--prediction_output ${predixcan_dir}/output/predixcan_test.txt \
--prediction_summary_output $${predixcan_dir}/output/predixcan_test_summary.txt \
--verbosity 9 \
--throw


temi_dir="/lus/grand/projects/covid-ct/imlab/users/temi"
predixcan_dir="/lus/grand/projects/covid-ct/imlab/users/temi/projects/prediXcan"

DATA=${predixcan_dir}/tutorial/data
RESULTS=${predixcan_dir}/tutorial/results

python3 ${temi_dir}/software/MetaXcan/software/Predict.py \
--model_db_path $DATA/models/gtex_v8_en/en_Whole_Blood.db \
--vcf_genotypes $DATA/1000G_hg38/ALL.chr*.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz \
--vcf_mode genotyped \
--variant_mapping $DATA/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz id rsid \
--on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" \
--prediction_output $RESULTS/vcf_1000G_hg38_en/Whole_Blood__predict.txt \
--prediction_summary_output $RESULTS/vcf_1000G_hg38_en/Whole_Blood__summary.txt \
--verbosity 9 \
--throw