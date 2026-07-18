#!/bin/bash

# Author: Temi
# Date: Wednesday June 18 2025
# Description: placeholder/scratch script -- just activates a conda env and echoes "Running test"; not part of the real pipeline
# Usage: sbatch test.sbatch

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