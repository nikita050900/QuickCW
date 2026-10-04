#!/bin/bash

#SBATCH --job-name=3C273_broad
#SBATCH --output=3C273_broad_UL_woJ1713_N1e9_thin10.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_3C273_broad.py

echo "Run complete."

