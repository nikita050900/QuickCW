#!/bin/bash

#SBATCH --job-name=1553113_targeted_288_2
#SBATCH --output=1553113_targeted_288E8_UL_woJ1713_N1e9_thin10.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_1553+113_2.88E_8_2.py

echo "Run complete."

