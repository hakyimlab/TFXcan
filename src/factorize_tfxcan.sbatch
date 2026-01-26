#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G    
#SBATCH --job-name=factorize_tfxcan
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/factorize_tfxcan.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/factorize_tfxcan.err
#SBATCH --time=06:00:00
#SBATCH --partition=caslake

date
source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/tfxcan-env
module load parallel

gwas_summary_statistic=${1}
tfxcan_summary_statistic=${2}
output_basename=${3}
base_output_directory=${4}

output_directory=${base_output_directory}/${output_basename}/results
temporary_directory=${base_output_directory}/${output_basename}/temp

mkdir -p ${output_directory}
mkdir -p ${temporary_directory}


if [ ! -f ${temporary_directory}/${output_basename}.random_subsets.txt ]; then
    Rscript /beagle3/haky/users/temi/projects/TFXcan/src/prepare_summary_matrices.R \
        --gwas_summary_statistic ${gwas_summary_statistic} \
        --tfxcan_summary_statistic ${tfxcan_summary_statistic} \
        --output_basename ${output_basename} \
        --output_directory ${output_directory} \
        --temporary_directory ${temporary_directory} \
        --number_resamples 1000 \
        --number_batches 10
    status=$?
    echo "Exit status of preparing data is: $status"
else
    echo "${temporary_directory}/${output_basename}.random_subsets.txt exists. Skipping..."
fi


if [ ! -f ${output_directory}/${output_basename}.*.*Iters.rds.gz ]; then
    if [ ${status} -eq 0 ]; then
        echo "INFO - Running consensus flashier..."
        batch_list=${temporary_directory}/${output_basename}.random_subsets.txt
        priorL='ebnm_point_exponential'
        priorF='ebnm_point_exponential'
        greedy_Kmax=40

        parallel -a ${batch_list} --jobs 20 \
            Rscript /beagle3/haky/users/temi/projects/TFXcan/src/repeat_flash.R \
            --data ${output_directory}/${output_basename}.zratios.matrices.rds \
            --splits ${temporary_directory}/${output_basename}.random_subsets.{}.rds \
            --batch {} \
            --output_basename ${temporary_directory}/${output_basename} \
            --priorL ${priorL} \
            --priorF ${priorF} \
            --greedy_Kmax ${greedy_Kmax}
    fi

    status=$?
    echo "Exit status of repeated factorization is: $status"

    if [ $status -eq 0 ]; then
        echo "INFO - Compiling the flashier results..."
        Rscript /beagle3/haky/users/temi/projects/TFXcan/src/gather_repeats.R \
            --input_basename ${temporary_directory}/${output_basename} \
            --output_basename ${output_directory}/${output_basename} \
            --batch_list ${temporary_directory}/${output_basename}.random_subsets.txt
    fi

    status=$?
    echo "Exit status of merging is: $status"
else 
    echo "`ls ${output_directory}/${output_basename}.*.*Iters.rds.gz` exists. Skipping..."
fi


# resFile=`ls ${output_directory}/${output_basename}.loci_assignments.*_norm.txt.gz`
# {opt$output_basename}.loci_assignments.{length(rFlash_results)}samples.{opt$normalization}_norm.txt.gz
if [ ! -f ${output_directory}/${output_basename}.loci_assignments.*_norm.txt.gz ]; then
    echo "INFO - Clustering the flashier results..."
    input_file=`ls ${output_directory}/${output_basename}.*.*Iters.rds.gz`
    echo "INFO - Found ${input_file}..."
    Rscript /beagle3/haky/users/temi/projects/TFXcan/src/cluster_programs.R \
        --flash_results ${input_file} \
        --output_basename ${output_directory}/${output_basename}
    status=$?
    echo "Exit status of Clustering is: $status"
else 
    echo "${output_directory}/${output_basename}.loci_assignments.*_norm.txt.gz exists. Skipping..."
fi

# echo "Exit status of extracting epigenome is: $status"

# sbatch /beagle3/haky/users/temi/projects/TFXcan/src/factorize_tfxcan.sbatch /beagle3/haky/users/temi/projects/TFXcan-snakemake/data/prostate_cancer_risk_2024-09-30/collection/prostate_cancer_risk.filteredGWAS.topSNPs.txt.gz /beagle3/haky/users/temi/projects/Enpact/misc/reruns/prca_tfxcan/prostate_cancer_risk.enpactScores.2025-10-21.complete.spredixcan.tsv TEMO /scratch/beagle3/temi/test /scratch/beagle3/temi/temp

# sbatch /beagle3/haky/users/temi/projects/TFXcan/src/factorize_tfxcan.sbatch /beagle3/haky/users/shared_data/from_david/9Jan2026/asthma_children_hg38_chr17_35000000_45000000_Format.filteredGWAS.topSNPs.txt.gz /beagle3/haky/users/shared_data/from_david/9Jan2026/asthma_children_hg38_chr17_35000000_45000000_Format.enpactScores.2025-11-18.spredixcan.txt asthma_children_hg38_chr17_35000000_45000000 /scratch/beagle3/temi