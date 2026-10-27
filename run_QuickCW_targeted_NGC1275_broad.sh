#!/bin/bash

#SBATCH --job-name=NGC1275_broad
#SBATCH --output=NGC1275_broad_UL_wJ1713_N1e9_thin10.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_NGC1275_broad.py

echo "Run complete."

