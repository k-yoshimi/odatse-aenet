#!/bin/sh
#SBATCH -p i8cpu
#SBATCH -N 1
#SBATCH -n 8
#SBATCH -t "01:00:00"
#SBATCH -e error.txt
#==========================================
# Step 1 & 2: Generate structures, run QE relaxation, and create training data
#
# Adjust partition (-p), node count, and QE module path for your cluster.
#==========================================

# Load Quantum ESPRESSO (adjust for your environment)
# source /path/to/espresso/espressovars.sh

python3 structure_make.py

for i in $(seq 0 19)
do
    cd "directory_${i}"
    srun -n 8 pw.x < n2_dimer.pwi > n2_dimer.pwo
    python3 ../teach_data_make.py --input n2_dimer.pwo --output-dir teach_data
    cd ../
done
