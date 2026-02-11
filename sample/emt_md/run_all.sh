#!/bin/bash
#==========================================
# EMT MD Tutorial - Full Execution Script
#
# Requirements: ASE
#==========================================
set -eu

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== EMT MD Tutorial ==="
echo "Working directory: $(pwd)"
echo ""

echo "--- Step 1: Run NVT MD simulation ---"
python3 run_md.py
echo ""

echo "=== Done ==="
echo "Output: md.traj"
