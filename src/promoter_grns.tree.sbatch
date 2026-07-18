#!/bin/bash
# Author: Temi
# Description: Runs promoter_grns.tree.R once per gene in gene_list, in parallel via GNU parallel —
#   builds promoter GRN fit stats (oob R^2) for every gene against one TF's predicted binding.
#   Paths below are hardcoded to the GATA1 gene-modules run; edit them for a different TF/dataset.
# Usage: sbatch promoter_grns.tree.sbatch   (no positional args — edit the hardcoded paths inside)

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=8G    
#SBATCH --job-name=promoter_grns_tree
#SBATCH --account=pi-haky
#SBATCH --output=/beagle3/haky/users/temi/projects/TFXcan/logs/promoter_grns_tree.out
#SBATCH --error=/beagle3/haky/users/temi/projects/TFXcan/logs/promoter_grns_tree.err
#SBATCH --time=24:00:00
#SBATCH --partition=caslake

date

source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/renviron
module load parallel/

bmat="/beagle3/haky/users/temi/projects/TFXcan/experiments/gene_modules/files/predictedTFbinding.GATA1.tsv.gz"
emat="/beagle3/haky/users/temi/projects/TFXcan/experiments/gene_modules/files/geuvadis.lcl.gene_expression.tsv.gz"
outfile_bname="/scratch/beagle3/temi/gene_modules/promoterGRN/GATA1/temp/temp"
gene_list="/beagle3/haky/users/temi/projects/TFXcan/experiments/gene_modules/files/genes_tss_locus.tsv"

# -a ${gene_list} -a ${gene_list} 
# run
parallel -j 100 --joblog /beagle3/haky/users/temi/projects/TFXcan/logs/promoter_grns_tree.job "Rscript /beagle3/haky/users/temi/projects/TFXcan/src/promoter_grns.tree.R --gene {} --binding_mat ${bmat} --expression_mat ${emat} --output_file_basename ${outfile_bname}" :::: <(cat ${gene_list})   #:::: <(head -10 ${gene_list})

status=$?
echo "Exit status of job is: $status"