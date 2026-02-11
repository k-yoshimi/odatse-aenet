#!/bin/bash
#==========================================
# Fullerene ELSES MD Tutorial - Full Execution Script
#
# Requirements: ELSES (elses, elses-xml-generate), ASE,
#               ElsesCalculator, matplotlib
#               imagemagick (for animation)
#
# Before running, ensure:
#   - elses and elses-xml-generate are in PATH
#   - C.cif is in this directory
#==========================================
set -eu

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== Fullerene ELSES MD Tutorial ==="
echo "Working directory: $(pwd)"
echo ""

if [ ! -f "C.cif" ]; then
    echo "ERROR: Required file 'C.cif' not found in $(pwd)"
    exit 1
fi

echo "--- Step 1: Run ELSES NVT MD simulation ---"
python3 test_elses.py
echo ""

echo "--- Step 2: Generate training data (XSF) ---"
python3 teach_data_make.py
echo ""

echo "--- Step 3: Plot kinetic energy ---"
python3 time_energy_plot.py
echo ""

echo "--- Step 4: Create trajectory animation ---"
if command -v convert >/dev/null 2>&1; then
    python3 animation.py --traj benzene_optimization.traj --output animation.gif \
        --natoms 480 --dt 10.0 --nframes 10
else
    echo "SKIP: imagemagick not found. Install it to generate animation.gif"
fi
echo ""

echo "=== Done ==="
echo "Outputs:"
echo "  md.log                     - MD log"
echo "  benzene_optimization.traj  - ASE trajectory"
echo "  teach_data/                - XSF training data"
echo "  time_energy_plot_QE.png    - Ekin plot"
echo "  animation.gif              - Trajectory animation (if imagemagick available)"
