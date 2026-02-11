#!/bin/sh
#SBATCH -p i8cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t "00:30:00"
#SBATCH -e error.txt
#==========================================
# Step 3: Generate fingerprints with AENET generate.x
#==========================================

generate.x generate.in > generate.out
