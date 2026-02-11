#!/bin/sh
#SBATCH -p i8cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t "00:30:00"
#SBATCH -e error.txt
#==========================================
# Step 5: Predict energies with AENET predict.x
#==========================================

predict.x predict.in > predict.out
