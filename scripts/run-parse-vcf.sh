#!/bin/bash

#SBATCH --job-name=PARSE_VCF
#SBATCH --partition=broadwl
#SBATCH --nodes=60
#SBATCH --ntasks-per-node=10
#SBATCH --time=24:00:00
#SBATCH --mem=32G

#SBATCH --error=/project2/haky/temi/projects/TFXcan/log/parse_vcf.err
#SBATCH --output=/project2/haky/temi/projects/TFXcan/log/parse_vcf.out


echo "Starting to parse VCF files ===\n"

python /project2/haky/temi/projects/TFXcan/scripts/parse-vcf.py

echo "\nDone with parsing VCF files ==="