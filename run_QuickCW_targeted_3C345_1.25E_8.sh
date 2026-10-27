#!/bin/bash

#SBATCH --job-name=3C345_125E_8
#SBATCH --output=3C345_125E_8_newMcmax_UL.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_3C345_1.25E_8.py

echo "Run complete."

