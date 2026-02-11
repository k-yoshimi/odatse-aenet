#!/bin/bash
#==========================================
# AENET Training Pipeline - Full Execution Script
#
# Requirements: Quantum ESPRESSO (pw.x), AENET (generate.x, train.x, predict.x),
#               ASE, matplotlib
#
# Configuration: Adjust variables below for your environment.
#==========================================
set -eu

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# --- Configuration ---
NPROCS_QE=8                              # Number of MPI processes for QE
QE_CMD="mpirun -np ${NPROCS_QE} pw.x"   # QE execution command
NSTRUCTS=20                              # Number of random structures

echo "=== AENET Training Pipeline ==="
echo "Working directory: $(pwd)"
echo "QE command: ${QE_CMD}"
echo ""

#------------------------------------------
# Step 1: Generate structures & run QE relaxation
#------------------------------------------
echo "--- Step 1: Generate structures and run QE ---"
cd "$SCRIPT_DIR/relax"

python3 structure_make.py

for i in $(seq 0 $((NSTRUCTS - 1))); do
    dir="directory_${i}"
    echo "  Running QE in ${dir} ..."
    cd "$dir"
    ${QE_CMD} < n2_dimer.pwi > n2_dimer.pwo
    cd ..
done
echo ""

#------------------------------------------
# Step 2: Convert QE output to XSF training data
#------------------------------------------
echo "--- Step 2: Create training data (XSF) ---"
cd "$SCRIPT_DIR/relax"

for i in $(seq 0 $((NSTRUCTS - 1))); do
    dir="directory_${i}"
    cd "$dir"
    python3 ../teach_data_make.py --input n2_dimer.pwo --output-dir teach_data
    cd ..
done
echo ""

#------------------------------------------
# Step 3: Generate fingerprints
#------------------------------------------
echo "--- Step 3: Generate fingerprints ---"
cd "$SCRIPT_DIR/generate"
generate.x generate.in > generate.out
echo ""

#------------------------------------------
# Step 4: Train neural network potential
#------------------------------------------
echo "--- Step 4: Train ANN potential ---"
cd "$SCRIPT_DIR/train"
ln -sf ../generate/N2.train .
ln -sf ../generate/N.fingerprint.stp .
train.x train.in > train.out
echo ""

#------------------------------------------
# Step 5: Predict energies
#------------------------------------------
echo "--- Step 5: Predict with trained potential ---"
cd "$SCRIPT_DIR/predict"
ln -sf ../train/N.5t-5t.ann .
predict.x predict.in > predict.out
echo ""

#------------------------------------------
# Step 6: Plot results
#------------------------------------------
echo "--- Step 6: Plot distance-energy curve ---"
cd "$SCRIPT_DIR/predict"
python3 plot_distance_energy.py
echo ""

echo "=== Done ==="
echo "Key outputs:"
echo "  generate/N2.train          - Fingerprint data"
echo "  train/N.5t-5t.ann          - Trained ANN potential"
echo "  predict/predict.out        - Prediction results"
echo "  predict/distance_energy_plot.png - Distance-energy plot"
