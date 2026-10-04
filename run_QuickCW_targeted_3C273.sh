#!/bin/bash

#SBATCH --job-name=3C273_targeted
#SBATCH --output=3C273_targeted_186E8_detect_woJ1713_N1e8_thin10.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_3C273.py

echo "Run complete."

