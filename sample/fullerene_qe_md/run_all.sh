#!/bin/bash
#==========================================
# Fullerene QE MD Tutorial - Full Execution Script
#
# Requirements: Quantum ESPRESSO (pw.x), ASE, matplotlib
#
# Before running, ensure:
#   - pw.x is in PATH (or adjust QE_COMMAND below)
#   - C.cif is in this directory
#   - C.pbe-n-kjpaw_psl.1.0.0.UPF is in this directory
#==========================================
set -eu

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== Fullerene QE MD Tutorial ==="
echo "Working directory: $(pwd)"
echo ""

# --- Prerequisite check ---
for f in C.cif C.pbe-n-kjpaw_psl.1.0.0.UPF; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required file '$f' not found in $(pwd)"
        exit 1
    fi
done

echo "--- Step 1: Run QE NVT MD simulation ---"
python3 ase_qe.py
echo ""

echo "--- Step 2: Plot kinetic energy ---"
python3 time_energy_plot.py
echo ""

echo "=== Done ==="
echo "Outputs:"
echo "  md.log                     - MD log"
echo "  benzene_optimization.traj  - ASE trajectory"
echo "  time_energy_plot_QE.png    - Ekin plot"
