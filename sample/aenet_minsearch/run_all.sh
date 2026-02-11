#!/bin/bash
#==========================================
# AENET Minsearch Tutorial - Full Execution Script
#
# Requirements: ODAT-SE (odatse-aenet), AENET (predict.x)
#
# Before running, ensure:
#   - predict.x is in PATH
#   - N.5t-5t.ann is in this directory
#   - template.xsf is in this directory
#==========================================
set -eu

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Auto-detect odatse-aenet: check PATH, then user site-packages bin
ODATSE_CMD="odatse-aenet"
if ! command -v odatse-aenet > /dev/null 2>&1; then
    USER_BIN="$(python3 -m site --user-base 2>/dev/null)/bin"
    if [ -x "${USER_BIN}/odatse-aenet" ]; then
        ODATSE_CMD="${USER_BIN}/odatse-aenet"
    else
        echo "ERROR: odatse-aenet not found. Install with: pip install odatse-aenet" >&2
        exit 1
    fi
fi

echo "=== AENET Minsearch Tutorial ==="
echo "Working directory: $(pwd)"
echo "odatse-aenet: ${ODATSE_CMD}"
echo ""

#------------------------------------------
# Step 0: Copy ANN potential from training output if not present
#------------------------------------------
if [ ! -f "N.5t-5t.ann" ]; then
    TRAIN_ANN="$SCRIPT_DIR/../aenet_training/train/N.5t-5t.ann"
    if [ -f "$TRAIN_ANN" ]; then
        echo "--- Step 0: Copying N.5t-5t.ann from training output ---"
        cp "$TRAIN_ANN" .
        echo ""
    else
        echo "ERROR: N.5t-5t.ann not found. Run the aenet_training tutorial first." >&2
        exit 1
    fi
fi

# --- Prerequisite check ---
for f in N.5t-5t.ann template.xsf input.toml predict.in; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required file '$f' not found in $(pwd)"
        exit 1
    fi
done

#------------------------------------------
# Step 1: Run minsearch (Nelder-Mead optimization)
#------------------------------------------
echo "--- Step 1: Run minsearch ---"
export PYTHONUNBUFFERED=1
time ${ODATSE_CMD} input.toml
echo ""

#------------------------------------------
# Step 2: Show results
#------------------------------------------
echo "--- Step 2: Optimization results ---"
if [ -f "output/res.txt" ]; then
    echo "=== res.txt ==="
    cat output/res.txt
    echo ""
fi
if [ -f "output/SimplexData.txt" ]; then
    echo "=== SimplexData.txt (last 5 lines) ==="
    tail -5 output/SimplexData.txt
    echo ""
fi

echo "=== Done ==="
echo "Outputs:"
echo "  output/res.txt         - Optimization result (optimal parameters)"
echo "  output/SimplexData.txt - Simplex iteration history"
