#!/bin/bash

#SBATCH --job-name=J050925
#SBATCH --output=J050925_newMcmax_UL.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_J050925.py

echo "Run complete."

