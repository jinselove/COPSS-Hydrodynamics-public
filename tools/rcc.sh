#!/bin/bash
#SBATCH --job-name=f_2
#SBATCH --output=f_2
#SBATCH --error=f_2
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jyli@uchicago.edu

bash run.sh

