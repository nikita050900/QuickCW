#!/bin/bash

#SBATCH --job-name=1928738_targeted
#SBATCH --output=1928738_targeted_104E8_UL_wJ1713_N1e9_thin10.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_1928+738_1.04E_8.py

echo "Run complete."

