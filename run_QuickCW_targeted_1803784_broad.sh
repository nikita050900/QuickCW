#!/bin/bash

#SBATCH --job-name=1803784_broad
#SBATCH --output=1803784_broad_UL_wJ1713_N1e9_thin10.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_1803784_broad.py

echo "Run complete."

