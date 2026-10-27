#!/bin/bash

#SBATCH --job-name=FBQS_J1159+2914
#SBATCH --output=FBQS_J1159+2914.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted.py

echo "Run complete."

