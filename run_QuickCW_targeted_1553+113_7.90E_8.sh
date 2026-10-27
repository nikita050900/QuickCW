#!/bin/bash

#SBATCH --job-name=1553113_targeted_1
#SBATCH --output=1553113_targeted_790E8_detect_woJ1713_N1e8_thin10.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_1553+113_7.90E_8.py

echo "Run complete."

