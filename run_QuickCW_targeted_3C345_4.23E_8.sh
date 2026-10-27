#!/bin/bash

#SBATCH --job-name=3C345_423E_8
#SBATCH --output=3C345_423E_8_newMcmax_UL.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_3C345_4.23E_8.py

echo "Run complete."

