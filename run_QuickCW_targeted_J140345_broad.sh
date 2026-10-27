#!/bin/bash

#SBATCH --job-name=J140345_broad
#SBATCH --output=J140345_broad_UL_wJ1713_N1e9_thin10.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_J140345_broad.py

echo "Run complete."

