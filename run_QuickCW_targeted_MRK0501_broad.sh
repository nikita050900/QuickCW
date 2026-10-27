#!/bin/bash

#SBATCH --job-name=MRK0501_broad
#SBATCH --output=MRK0501_broad_UL_woJ1713_N1e9_thin10.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_MRK0501_broad.py

echo "Run complete."

