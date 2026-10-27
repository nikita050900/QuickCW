#!/bin/bash

#SBATCH --job-name=1553113_broad_UL
#SBATCH --output=1553113_broad_UL__woJ1713_N1e9_thin10.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_1553113_broad.py

echo "Run complete."

