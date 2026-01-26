#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G    
#SBATCH --job-name=find_motifs_given_bed
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/find_motifs_given_bed.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/find_motifs_given_bed.err
#SBATCH --time=04:00:00	
#SBATCH --partition=caslake


date

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools
export PATH=$PATH:/beagle3/haky/users/temi/software/homer/bin

perl '/beagle3/haky/users/temi/software/homer/bin/findMotifsGenome.pl' ${1} hg38 ${2} -size ${3} -bg ${4}



