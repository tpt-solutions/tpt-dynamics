#!/usr/bin/env python3
"""Write a small synthetic event CSV for the `tptdyn analyze` and `histogram` examples.

Standard library only. Each row is one particle candidate with columns:

    event, pt, eta, phi, weight

`pt` is exponentially distributed (GeV), `eta` is uniform in [-2.5, 2.5], `phi` is
uniform in [-pi, pi), and `weight` is 1.0 for every row. The output is deterministic
for a given seed.

Usage:

    python make_events.py                       # writes events.csv, 2000 rows
    python make_events.py -n 50000 -o big.csv   # more rows, another name
"""

from __future__ import annotations

import argparse
import csv
import math
import random
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("-n", "--rows", type=int, default=2000, help="number of rows")
    parser.add_argument("-o", "--output", type=Path, default=Path("events.csv"))
    parser.add_argument("--seed", type=int, default=2024, help="random seed")
    parser.add_argument("--mean-pt", type=float, default=15.0, help="mean pt in GeV")
    args = parser.parse_args()

    if args.rows < 1:
        parser.error("--rows must be at least 1")

    rng = random.Random(args.seed)
    with args.output.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["event", "pt", "eta", "phi", "weight"])
        for i in range(args.rows):
            pt = rng.expovariate(1.0 / args.mean_pt)
            eta = rng.uniform(-2.5, 2.5)
            phi = rng.uniform(-math.pi, math.pi)
            writer.writerow([i // 4, f"{pt:.6f}", f"{eta:.6f}", f"{phi:.6f}", "1.0"])

    print(f"wrote {args.rows} rows to {args.output}")


if __name__ == "__main__":
    main()
