"""Write a synthetic event table for the tptdyn analysis examples.

Usage: python make_events.py [--rows N] [--seed S] [--out events.csv]

Uses only the Python standard library. Columns: event, pt, eta, phi, n_jets.
The transverse momentum is exponentially falling (a rough shape, not physics),
pseudorapidity is uniform in [-2.5, 2.5], azimuth is uniform in [-pi, pi].
"""

import argparse
import csv
import math
import random


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--rows", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=2024)
    parser.add_argument("--out", default="events.csv")
    args = parser.parse_args()

    rng = random.Random(args.seed)
    with open(args.out, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["event", "pt", "eta", "phi", "n_jets"])
        for event in range(args.rows):
            pt = rng.expovariate(1 / 25.0)
            eta = rng.uniform(-2.5, 2.5)
            phi = rng.uniform(-math.pi, math.pi)
            n_jets = min(8, int(rng.expovariate(0.9)))
            writer.writerow([event, f"{pt:.4f}", f"{eta:.4f}", f"{phi:.4f}", n_jets])
    print(f"wrote {args.rows} rows to {args.out}")


if __name__ == "__main__":
    main()
