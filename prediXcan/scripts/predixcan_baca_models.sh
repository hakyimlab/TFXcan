#!/bin/bash
#PBS -l select=1:system=polaris
#PBS -l walltime=06:00:00,filesystems=home:grand
#PBS -A covid-ct
#PBS -q preemptable   
#PBS -N predixcan_baca_geuvadis_240
#PBS -k doe
#PBS -o /grand/projects/covid-ct/imlab/users/temi/projects/prediXcan/log/predixcan_baca_geuvadis_240.out
#PBS -e /grand/projects/covid-ct/imlab/users/temi/projects/prediXcan/log/predixcan_baca_geuvadis_240.err

echo Working directory is $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

echo Jobid: $PBS_JOBID
echo Running on host `hostname`
echo Running on nodes `cat $PBS_NODEFILE`

NNODES=`wc -l < $PBS_NODEFILE`
NRANKS=1
NDEPTH=48
NTHREADS=2

#NTOTRANKS=$(( NNODES * NRANKS ))
NTOTRANKS=$(( NNODES * NRANKS_PER_NODE ))
echo "NUM_OF_NODES=${NNODES}  TOTAL_NUM_RANKS=${NTOTRANKS}  RANKS_PER_NODE=${NRANKS}  THREADS_PER_RANK=${NTHREADS}"
echo "PBS_JOBID = " $PBS_JOBID
printf "Starting to run\n"

source ~/.bashrc
temi_dir="/lus/grand/projects/covid-ct/imlab/users/temi"
predixcan_dir="/lus/grand/projects/covid-ct/imlab/users/temi/projects/prediXcan"
conda activate imlabtools

for model_db in `ls ${temi_dir}/projects/TFXcan/baca_cwas/db_folder/`; do
    m_name=$( echo ${model_db} | cut -d '_' -f 3 )
    m_name=${m_name%.*}
    if [[ ${m_name} != 'top1' ]]; then
        continue
    else
        echo ${m_name}

        ${temi_dir}/software/MetaXcan/software/Predict.py \
        --model_db_path "${temi_dir}/projects/TFXcan/baca_cwas/db_folder/baca_cwas_${m_name}.db" \
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

printf "[INFO] Finished with all models."