#!/bin/sh
#SBATCH -p i8cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t "00:30:00"
#SBATCH -e error.txt
#==========================================
# Step 4: Train neural network potential with AENET train.x
#==========================================

train.x train.in > train.out
