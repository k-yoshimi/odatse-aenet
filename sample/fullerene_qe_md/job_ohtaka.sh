#!/bin/sh
#SBATCH -p F16cpu
#SBATCH -N 16
#SBATCH -n 64
#SBATCH -t "20:00:00"
#SBATCH --mail-type=BEGIN,END
#SBATCH -e error.txt
#==========================================
# Quantum ESPRESSO NVT MD on ISSP Ohtaka cluster
#
# Adjust the partition (-p), node count (-N), and task count (-n)
# for your system. Load the appropriate QE module before running.
#==========================================

source /home/issp/materiapps/oneapi_compiler_classic-2023.0.0--openmpi-4.1.5/espresso/espressovars-7.1-1.sh

python3 ase_qe.py
