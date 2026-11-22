#!/usr/bin/env python3
"""Fill a 1D histogram from a CSV file with the tpt_dynamics Python API.

Reads a CSV with a header (for example one written by make_events.py), keeps the
rows with pt above a threshold using a lazy frame, fills a Histogram1D from the
pt column, and prints the bin contents as a text bar chart.

Requires the tpt_dynamics package (build it with `maturin develop` in
bindings/tpt-dynamics-py) and numpy:

    python python_histogram.py events.csv --pt-min 20 --bins 40 --lo 20 --hi 100
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from tpt_dynamics import Histogram1D, LazyFrame, col


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("csv", type=Path, help="input CSV with a 'pt' column")
    parser.add_argument("--pt-min", type=float, default=20.0, help="keep pt > this value")
    parser.add_argument("--bins", type=int, default=25)
    parser.add_argument("--lo", type=float, default=20.0)
    parser.add_argument("--hi", type=float, default=100.0)
    args = parser.parse_args()

    # Nothing is read until collect(). The filter runs in the lazy engine.
    selected = LazyFrame.scan_csv(str(args.csv)).filter(col("pt") > args.pt_min).collect()
    pt = np.asarray(selected["pt"], dtype=np.float64)

    hist = Histogram1D(args.bins, args.lo, args.hi)
    hist.fill(pt)

    values, edges = hist.to_numpy()
    print(f"selected {len(pt)} rows with pt > {args.pt_min}")
    print(f"entries {hist.entries}, underflow {hist.underflow}, overflow {hist.overflow}")
    peak = max(float(values.max()), 1.0)
    for lo_edge, hi_edge, content in zip(edges[:-1], edges[1:], values):
        bar = "#" * round(40 * content / peak)
        print(f"[{lo_edge:7.2f}, {hi_edge:7.2f}) {content:8.0f} {bar}")


if __name__ == "__main__":
    main()
