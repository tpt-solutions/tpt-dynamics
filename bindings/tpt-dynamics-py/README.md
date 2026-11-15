# tpt-dynamics (Python)

Python bindings for [tpt-dynamics](../../README.md): lazy data frames, physics
histograms and molecular dynamics runs, implemented in Rust (PyO3) and exposed
as the `tpt_dynamics` package. Licensed under MIT OR Apache-2.0.

The package ships as an abi3 wheel (Python 3.11 and later, one wheel per
platform). It depends only on NumPy.

## Install (from source)

```bash
pip install maturin
maturin develop --release        # from bindings/tpt-dynamics-py/
```

## Histograms

```python
import numpy as np
import tpt_dynamics as td

h = td.Histogram1D(50, 0.0, 200.0)            # 50 uniform bins over [0, 200)
h.fill(np.random.normal(91.2, 5.0, 10_000))   # numpy array of values
h.fill(x, weights=w)                          # optional per-entry weights

h.contents()        # weighted bin contents (numpy float64)
h.errors()          # sqrt(sum of w^2) per bin
h.edges()           # bin edges, length nbins + 1
h.underflow, h.overflow, h.entries
h.mean(), h.std_dev(), h.quantile(0.5)

other = td.Histogram1D(50, 0.0, 200.0)
h.merge(other)      # identical binning required; additive
d = h.to_dict()     # plain dict of arrays and scalars

v = td.Histogram1D.variable([0.0, 10.0, 50.0, 200.0])

h2 = td.Histogram2D(40, 0.0, 200.0, 40, -5.0, 5.0)
h2.fill(mass, eta)
h2.contents()       # shape (40, 40), indexed [ix, iy]
h2.project_x()      # Histogram1D

# Hist-style interop (boost-histogram / hist naming)
h.axes[0].edges     # Axis.edges, numpy array; .centers and len() also available
h.values()          # bin contents (same as contents())
h.variances()       # per-bin sum of w^2 (errors() is its square root)
values, edges = h.to_numpy()
h3 = td.Histogram1D.from_numpy(values, edges)                 # variances = values
h4 = td.Histogram1D.from_numpy(values, edges, variances=var)  # weighted import

h2.values()         # shape (nbins_x, nbins_y)
h2.variances()
values, xedges, yedges = h2.to_numpy()
h5 = td.Histogram2D.from_numpy(values, xedges, yedges)
xaxis, yaxis = h2.axes
```

`from_numpy` builds a histogram from bin contents and edges, so values from
another library can be imported. Flow counts and fill counts start at zero, so
`mean()` and `std_dev()` return `None` until new fills arrive. Negative or
non-finite values, and variances that are negative or non-finite, raise
`ValueError`.

## Lazy data frames

A `LazyFrame` records a plan; `collect()` runs it and returns a dict of numpy
arrays. Filters and projections are pushed down, and the data is processed in
chunks.

```python
import numpy as np
import tpt_dynamics as td
from tpt_dynamics import col, Agg

lf = td.LazyFrame.scan_dict({
    "pt": pt,                        # float64 arrays
    "eta": eta,
    "phi": phi,
    "charge": charge,                # int64 arrays
})

result = (
    lf.filter((col("pt") > 20.0) & (col("eta").gt(-2.5)) & (col("eta") < 2.5))
      .with_column((col("pt") / 1000.0).alias("pt_gev"))
      .group_by(["charge"])
      .agg([Agg.count().alias("n"), Agg.mean(col("pt_gev")).alias("mean_pt")])
      .collect()
)
result["n"], result["mean_pt"], result["charge"]

# Four-vector helpers
mass = td.invariant_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2)
dr = td.delta_r(eta1, phi1, eta2, phi2)
```

Scan a CSV file (parsed in full when the frame is created; columns default to
float64, others can be given as `"i64"` or `"bool"`):

```python
lf = td.LazyFrame.scan_csv("events.csv", dtypes={"njets": "i64"})
```

Use `lf.explain()` to print the optimised plan.

## Molecular dynamics

```python
import tpt_dynamics as td

config = """
seed = 42
units = "reduced"
dt = 0.002
nsteps = 1000
temperature = 0.7
integrator = "velocity-verlet"

[neighbour]
cutoff = 2.5
skin = 0.3

[lennard_jones]
epsilon = 1.0
sigma = 1.0
modifier = "shift"
"""

td.validate_config(config)                      # raises ValueError if invalid
gro_text = open("water_box.gro").read()
energies = td.run_simulation(config, gro_text, format="gro", sample_every=10)
energies["time"], energies["total"], energies["temperature"]
```

`run_simulation` takes the structure as text (GRO or PDB), so no file is needed
for the input. It returns `step` (uint64) and float64 arrays `time`, `kinetic`,
`potential`, `total` and `temperature`. The output-file settings in the
configuration are ignored, because the arrays are returned instead.

## Errors

Invalid configurations, unknown columns, mismatched lengths and bad CSV values
raise `ValueError`. Out-of-range bin indices raise `IndexError`. Missing input
files raise `OSError`. Passing a non-numpy column raises `TypeError`.

## Development

```bash
python -m venv .venv && . .venv/bin/activate   # or .venv\Scripts\activate on Windows
pip install maturin numpy pytest
maturin develop
pytest tests -q
```

The Rust side is checked with `cargo test -p tpt-dynamics-py` and
`cargo clippy -p tpt-dynamics-py --all-targets -- -D warnings` from the
workspace root. See [ADR-007](../../docs/adr/ADR-007-pyo3-unsafe-exception.md)
for the unsafe-code exception that applies to this crate only.

For users coming from ROOT/PyROOT, see
[docs/migration/pyroot.md](../../docs/migration/pyroot.md).
