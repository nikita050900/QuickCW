#!/bin/bash

#SBATCH --job-name=0954658_targeted
#SBATCH --output=0954658_targeted_334E8_UL_wJ1713_N1e9_thin10.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_0954+658.py

echo "Run complete."

