#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH ---cpus-per-task=6
#SBATCH --mem-per-cpu=12G    
#SBATCH --job-name=liftover
#SBATCH --account=pi-haky
#SBATCH --output=/project/haky/users/temi/projects/Enpact/logs/liftover.out
#SBATCH --error=/project/haky/users/temi/projects/Enpact/logs/liftover.err
#SBATCH --time=02:00:00	
#SBATCH --partition=caslake


date

module load parallel

job_log='/lus/grand/projects/TFXcan/imlab/users/temi/projects/TFXcan/logs/liftover_beds.log' 

liftover_cmd='/lus/grand/projects/TFXcan/imlab/users/temi/software/liftOver'
chain_file='/lus/grand/projects/covid-ct/imlab/data/liftover_files/hg19ToHg38.over.chain.gz'
to_folder='/lus/grand/projects/TFXcan/imlab/data/baca_cwas/liftover_hg38'
mkdir -p ${to_folder}

# gather the input bed files
input_beds=( `find /lus/grand/projects/TFXcan/imlab/data/baca_cwas/sorted_bed -type f` )

echo ${input_beds[@]}

function lift_hg19_hg38(){
    input_bed=${1}
    liftover_cmd=${2}
    chain_file=${3}
    to_folder=${4}
    # get filename from input bed
    fname=$( echo ${input_bed} | rev | cut -d '/' -f 1 | rev )
    printf "\nINFO - started with ${fname}\n"
    ${liftover_cmd} ${input_bed} ${chain_file} ${to_folder}/${fname} ${to_folder}/unmapped_${fname}
    printf "\nINFO - finished with ${fname}\n"
}

export -f lift_hg19_hg38

parallel -j 12 --joblog ${job_log} "lift_hg19_hg38 {1} ${liftover_cmd} ${chain_file} ${to_folder}" ::: ${input_beds[@]}

#${mpiexec} -n ${NRANKS} --ppn ${NRANKS} --depth ${NDEPTH} --cpu-bind depth --env OMP_NUM_THREADS="${NTHREADS}" "${liftover_cmd}" "${input_bed}" "${chain_file}" "${output_bed}" "${unmapped_bed}"

status=$?

printf "\nINFO - Finished lifting over all files\n"
echo "Exit status of lifting over is: $status"

# qsub -v 'data_file=/lus/grand/projects/covid-ct/imlab/users/temi/projects/TFXcan/modeling_pipeline/data/train-test-val/kawakami/data_2022-12-12/kawakami_aggByCenter_FOXA1_old.csv.gz,metainfo=old' train_enet_model_pbs.sh