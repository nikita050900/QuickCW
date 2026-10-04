#!/bin/bash

#SBATCH --job-name=3C345_746E_9
#SBATCH --output=3C345_746E_9_newMcmax_UL.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_3C345_7.46E_9.py

echo "Run complete."

