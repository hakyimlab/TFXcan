#!/bin/bash
# Author: Temi
# Date: Monday January 26 2026
# Description: de novo motif discovery on a BED file with HOMER's findMotifsGenome.pl (hg38)
# Usage: sbatch find_motifs_given_bed.sbatch <bed_file> <output_dir> <promoter_size> <background_file>

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

# $1=bed file, $2=output dir, $3=promoter/region size, $4=background bed file
perl '/beagle3/haky/users/temi/software/homer/bin/findMotifsGenome.pl' ${1} hg38 ${2} -size ${3} -bg ${4}



