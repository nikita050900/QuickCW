#!/bin/bash

#SBATCH --job-name=1928+738_167E_8
#SBATCH --output=1928+738_167E_8_newMcmax_UL.out
#SBATCH -p sbs0016
#SBATCH --mem-per-cpu=64G
#SBATCH --ntasks=1

which python

python runQuickMCMC_targeted_1928+738_1.67E_8.py

echo "Run complete."

