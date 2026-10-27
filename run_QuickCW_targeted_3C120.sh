#!/bin/bash

#SBATCH --job-name=3C120
#SBATCH --output=3C120_newMcmax_UL.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_3C120.py

echo "Run complete."

