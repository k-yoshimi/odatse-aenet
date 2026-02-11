#!/bin/bash
#==========================================
# AENET PAMC Tutorial - Full Execution Script
#
# Requirements: ODAT-SE (odatse-aenet), AENET (predict.x),
#               MPI, matplotlib
#
# Before running, ensure:
#   - predict.x is in PATH
#   - N.5t-5t.ann is in this directory
#   - template.xsf is in this directory
#
# Configuration: Adjust variables below for your environment.
#==========================================
set -eu

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# --- Configuration ---
NPROCS=2   # Number of MPI processes for PAMC

echo "=== AENET PAMC Tutorial ==="
echo "Working directory: $(pwd)"
echo "MPI processes: ${NPROCS}"
echo ""

# --- Prerequisite check ---
for f in N.5t-5t.ann template.xsf input.toml predict.in; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required file '$f' not found in $(pwd)"
        exit 1
    fi
done

#------------------------------------------
# Step 1: Run PAMC optimization
#------------------------------------------
echo "--- Step 1: Run PAMC ---"
export PYTHONUNBUFFERED=1
export OMPI_MCA_rmaps_base_oversubscribe=true
time mpiexec -np ${NPROCS} odatse-aenet input.toml
echo ""

#------------------------------------------
# Step 2: Summarize results
#------------------------------------------
echo "--- Step 2: Summarize results ---"
python3 summarize_results.py -i input.toml --idnum 0
echo ""

#------------------------------------------
# Step 3: Plot histograms
#------------------------------------------
echo "--- Step 3: Plot histograms ---"
python3 batch_plot.py
echo ""

echo "=== Done ==="
echo "Outputs:"
echo "  output/           - PAMC output directories"
echo "  *_eachT/          - Per-temperature summary files"
echo "  *.png             - Histogram plots"
