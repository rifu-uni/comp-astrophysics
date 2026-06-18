#!/bin/bash
# Run the full Alfven wave pipeline
set -e

HERE="$(cd "$(dirname "$0")/.." && pwd)"
PLUTO_DIR="$(cd "$HERE/../../../tools/pluto" && pwd)"

echo "=== 1. Python Solver ==="
cd "$HERE"
pixi run python scripts/alfven_solver.py

echo ""
echo "=== 2. PLUTO CP_Alfven ==="
BUILD="$HERE/pluto/build"
mkdir -p "$BUILD"
cp "$HERE/pluto/pluto.ini" "$BUILD/"
cp "$HERE/pluto/definitions.h" "$BUILD/"
cp "$PLUTO_DIR/Test_Problems/MHD/CP_Alfven/init.c" "$BUILD/"
cd "$BUILD"
python3 "$PLUTO_DIR/setup.py" --auto-update
make -j4
export PLUTO_DIR="$PLUTO_DIR"
./pluto -i pluto.ini

echo ""
echo "=== 3. pyPLUTO Analysis ==="
cd "$HERE"
pixi run python pluto/plot_pluto.py

echo ""
echo "=== 4. Comparison ==="
cd "$HERE"
pixi run python scripts/compare.py

echo ""
echo "=== Done ==="
ls -lh "$HERE"/*.png "$HERE"/*.gif "$HERE"/*.npz "$HERE"/pluto/build/data.*.vtk 2>/dev/null
