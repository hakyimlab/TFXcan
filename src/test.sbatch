#!/bin/bash
#SBATCH --job-name=test
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --mem-per-cpu=1G
#SBATCH --error=logs/lEnpact_training.err
#SBATCH --output=logs/lEnpact_training.log
#SBATCH --account=pi-haky

date
source ~/.bashrc
conda activate /beagle3/haky/users/temi/software/conda_envs/predictdb-env

echo "Running test"