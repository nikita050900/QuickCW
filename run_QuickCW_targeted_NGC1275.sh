#!/bin/bash

#SBATCH --job-name=NGC1275_targeted
#SBATCH --output=NGC1275_targeted_671E8_UL_wJ1713_N1e9_thin10.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_NGC1275.py

echo "Run complete."

