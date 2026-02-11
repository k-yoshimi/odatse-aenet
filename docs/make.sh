#!/bin/sh

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Build Japanese documentation
echo "Building Japanese documentation..."
python3 -m sphinx -b html ja/source ja/build/html
echo "  -> docs/ja/build/html/"

# Build English documentation
echo "Building English documentation..."
python3 -m sphinx -b html en/source en/build/html
echo "  -> docs/en/build/html/"

echo "Done."
