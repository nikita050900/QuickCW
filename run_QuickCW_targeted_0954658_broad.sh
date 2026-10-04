#!/bin/bash

#SBATCH --job-name=0954658_broad
#SBATCH --output=0954658_broad_UL_wJ1713_N1e9_thin10.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_0954658_broad.py

echo "Run complete."

