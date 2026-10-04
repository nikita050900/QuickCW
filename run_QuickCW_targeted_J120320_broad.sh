#!/bin/bash

#SBATCH --job-name=J120320_broad_UL
#SBATCH --output=J12032061_broad_UL_woJ1713_N1e9_thin10.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_J120320_broad.py

echo "Run complete."

