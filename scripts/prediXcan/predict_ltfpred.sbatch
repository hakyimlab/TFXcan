#!/bin/bash

# Author:
# Description:
# Usage:
# Date:

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --job-name=split_vcf_files
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/split_vcf_files.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/split_vcf_files.err
#SBATCH --time=06:00:00	
#SBATCH --partition=caslake


module load openmpi
module load parallel

slurm_workdir=${SLURM_SUBMIT_DIR}
SLURM_O_WORKDIR=${slurm_workdir}

mkdir -p ${SLURM_O_WORKDIR}
echo Working directory is $SLURM_O_WORKDIR
cd $SLURM_O_WORKDIR

echo Jobid: $SLURM_JOBID
echo Running on host `hostname`

printf "Starting to run\n"

source ~/.bashrc
conda activate compbio-tools

output_folder="//beagle3/haky/users/temi/projects/TFXcan/experiments/compare_predictors"
mkdir -p ${output_folder}
vcf_file="/beagle3/haky/data/CWAS/merged/merged_cwas_genotypes_chrRenamed.vcf.gz"

printf "Starting to run\n"

source ~/.bashrc
temi_dir="/lus/grand/projects/covid-ct/imlab/users/temi"
predixcan_dir="/lus/grand/projects/covid-ct/imlab/users/temi/projects/prediXcan"
conda activate imlabtools

# arguments
db_name="gtex_v7_EUR_AR_Prostate_filtered_signif"
model_db="/beagle3/haky/users/temi/projects/l-TFPred/results/filtered_db/${db_name}.db"
predixcan_exec="/beagle3/haky/users/temi/software/MetaXcan/software/Predict.py"

${predixcan_exec} \
--model_db_path "${model_db}" \
--vcf_genotypes /lus/grand/projects/covid-ct/imlab/data/GEUVADIS/vcf_snps_only/ALL.chr*.shapeit2_integrated_SNPs_v2a_27022019.GRCh38.phased.vcf.gz \
--vcf_mode genotyped \
--on_the_fly_mapping METADATA "{}_{}_{}_{}_b38" \
--variant_mapping "${predixcan_dir}/tutorial/data/gtex_v8_eur_filtered_maf0.01_monoallelic_variants.txt.gz" id rsid \
--prediction_output "${predixcan_dir}/output/geuvadis_240/${m_name}/baca_cwas_predict.txt" \
--prediction_summary_output "${predixcan_dir}/output/geuvadis_240/${m_name}/baca_cwas_summary.txt" \
--verbosity 9 \
--text_sample_ids "/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/TFPred_pipeline/metadata/predixcan_geuvadis_240.txt"
--throw & sleep 1
    fi
done
wait

printf "INFO - Finished with all models."