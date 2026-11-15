#!/usr/bin/env bash
# Strong and weak scaling for the LJ engines (tpt-dynamics-sim examples/lj_scaling.rs).
# Usage (from the repository root, Git Bash or any POSIX shell):
#   bash tpt-dynamics-sim/benchmarks/run_scaling.sh > docs/benchmarks/scaling-raw.txt
# Build first: cargo build --release -p tpt-dynamics-sim --examples
set -euo pipefail
BIN=target/release/examples/lj_scaling
STEPS=10
WARMUP=2
THREADS="1 2 4 8 16 20"

echo "# strong scaling: fixed 1,000,000 atoms (100x100x100 simple cubic), domain engine"
for t in $THREADS; do
  "$BIN" --mode domain --side 100,100,100 --threads "$t" --steps "$STEPS" --warmup "$WARMUP"
done

echo "# strong scaling: fixed 1,000,000 atoms, single-domain engine (rayon pool of t threads)"
for t in $THREADS; do
  "$BIN" --mode single --side 100,100,100 --threads "$t" --steps "$STEPS" --warmup "$WARMUP"
done

echo "# single-domain profile on the final configuration (1 thread and 20 threads)"
"$BIN" --mode single --side 100,100,100 --threads 1 --steps "$STEPS" --warmup "$WARMUP" --profile
"$BIN" --mode single --side 100,100,100 --threads 20 --steps "$STEPS" --warmup "$WARMUP" --profile

echo "# weak scaling: 50,000 atoms per thread, box 50 x 50 x (20 t) cells, domain engine"
for t in $THREADS; do
  "$BIN" --mode domain --side "50,50,$((20 * t))" --threads "$t" --steps "$STEPS" --warmup "$WARMUP"
done
