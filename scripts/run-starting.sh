#!/bin/bash

#SBATCH --job-name=STARTING
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=4:00:00
#SBATCH --mem=32G

#SBATCH --error=/project2/haky/temi/projects/TFXcan/log/starting.err
#SBATCH --output=/project2/haky/temi/projects/TFXcan/log/starting.out


module load R/4.2.0
echo "Loaded R==="

Rscript /project2/haky/temi/projects/TFXcan/scripts/starting.R

echo "Done ==="