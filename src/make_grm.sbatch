#!/bin/bash
# Author: Temi
# Date: Wednesday June 18 2025
# Description: build a genome-wide GRM (genetic relationship matrix) for the Baca cohort from the merged
#   plink bfiles, using GCTA's --make-grm. Paths are hardcoded, not passed in as args.
# Usage: sbatch make_grm.sbatch

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4G    
#SBATCH --job-name=make_grm
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/Enpact/logs/make_grm.out
#SBATCH --error=/beagle3/haky/users/temi/projects/Enpact/logs/make_grm.err
#SBATCH --time=02:00:00
#SBATCH --partition=caslake

date

module load parallel
module load plink/1.90b6.21

# earlier step: merge per-chromosome plink files into one bfile set (left here for reference, not run)
# plink --merge-list "/beagle3/haky/users/temi/projects/Enpact/misc/mergelist.txt" --make-bed --out /project2/haky/Data/baca_cwas/vcfs/hg38/merged_files/bfiles/baca_individuals.plink.merged


gcta_executable="/beagle3/haky/users/temi/software/gcta-1.94.4-linux-kernel-3-x86_64/gcta64"

${gcta_executable} --bfile "/project2/haky/Data/baca_cwas/vcfs/hg38/merged_files/bfiles/baca_individuals.plink.merged" --make-grm --out "/project2/haky/Data/baca_cwas/vcfs/hg38/merged_files/grm/baca_individuals.plink.merged"

# earlier version of this script: took input_directory/output_directory/input_basenames as args and built
# one GRM per basename in parallel (this is basically what make_grm_single.sbatch does now)
# input_directory=${1}
# output_directory=${2}
# input_basenames=${3}

# IFS=',' read -r -a input_array <<< "${input_basenames}"

# gcta_executable="/beagle3/haky/users/temi/software/gcta-1.94.4-linux-kernel-3-x86_64/gcta64"

# parallel -j 5 "
#   ${gcta_executable} --bfile ${input_directory}/{1} --make-grm --out ${output_directory}/{1}" ::: ${input_array[@]}

status=$?
printf "INFO - All jobs are done. Exit status is $status.\n"



# for i in {1..22}
# do
# echo "/project2/haky/Data/baca_cwas/vcfs/hg38/plink_geno/chr$i.dose" >> "/beagle3/haky/users/temi/projects/Enpact/misc/mergelist.txt"
# done

# /project2/haky/Data/baca_cwas/vcfs/hg38/merged_files/bfiles

# plink --merge-list "/beagle3/haky/users/temi/projects/Enpact/misc/mergelist.txt" --make-bed --out /project2/haky/Data/baca_cwas/vcfs/hg38/merged_files/bfiles/baca_individuals.plink.merged