#!/bin/bash
# Author: Temi
# Description: annotates an existing BED file with known-motif hits using HOMER's annotatePeaks.pl (given motif files),
#   as opposed to find_motifs_given_bed.sbatch which does de novo motif discovery
# Usage: sbatch annotate_bed_with_motifs.sbatch <input_bed> <output_bed_with_motifs> <output_annotation_file> <motif_file1> [motif_file2 ...]

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --job-name=annotate_bed_with_motifs
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/annotate_bed_with_motifs.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/annotate_bed_with_motifs.err
#SBATCH --time=00:10:00
#SBATCH --partition=caslake

date

source ~/.bashrc
conda activate /beagle3/haky/users/shared_software/TFXcan-pipeline-tools
export PATH=$PATH:/beagle3/haky/users/temi/software/homer/bin

input_bed_file=${1} #  /beagle3/haky/users/temi/projects/TFXcan/experiments/motifs/gwas_loci.512bp.bed
output_bed_file=${2} # /beagle3/haky/users/temi/projects/TFXcan/experiments/motifs/prca_gwas_motifs_instances.bed
output_annotation_file=${3} # /beagle3/haky/users/temi/projects/TFXcan/experiments/motifs/prca_gwas_annotation_with_motifs.txt
motif_files=("$@") # /beagle3/haky/users/temi/software/homer/data/knownTFs/motifs/are.motif /beagle3/haky/users/temi/software/homer/data/knownTFs/motifs/ar-half.motif

perl /beagle3/haky/users/temi/software/homer/bin/annotatePeaks.pl ${input_bed_file} hg38 \
    -m "${motif_files[@]}" \
    -noann \
    -nmotifs \
    -mscore \
    -mbed ${output_bed_file} \
    > ${output_annotation_file}


# input_bed_file=${1} #  /beagle3/haky/users/temi/projects/TFXcan/experiments/motifs/gwas_loci.512bp.bed
# output_directory=${2} # /beagle3/haky/users/temi/projects/TFXcan/experiments/motifs/prca_gwas_motifs_instances.bed
# outfile=${3}
# motif_files=${4} # /beagle3/haky/users/temi/software/homer/data/knownTFs/motifs/are.motif /beagle3/haky/users/temi/software/homer/data/knownTFs/motifs/ar-half.motif


# perl /beagle3/haky/users/temi/software/homer/bin/findMotifsGenome.pl ${input_bed_file} hg38 ${output_directory} -find "${motif_files[@]}" > ${outfile}



