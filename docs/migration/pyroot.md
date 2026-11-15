# Migrating from PyROOT to `tpt_dynamics`

This page maps common PyROOT idioms to their `tpt_dynamics` equivalents. The
examples are the kind of analysis code that is typical of ROOT-based HEP work.
Every `tpt_dynamics` call shown here is covered by the test suite in
`bindings/tpt-dynamics-py/tests/test_basic.py`.

## Quick reference

| PyROOT | tpt_dynamics | Notes |
|--------|--------------|-------|
| `TH1F("h", "title", 50, 0, 200)` | `td.Histogram1D(50, 0.0, 200.0)` | Bins are `[lo, hi)`. Titles and names are not stored. |
| `TH1F(..., nbins, array_of_edges)` | `td.Histogram1D.variable(edges)` | `edges` is a sequence of strictly increasing floats. |
| `h.Fill(x)` | `h.fill(x_array)` | Takes a whole numpy array at once, which avoids a Python loop per entry. |
| `h.Fill(x, w)` | `h.fill(x_array, weights=w_array)` | |
| `h.GetBinContent(i)` | `h.bin(i)` | Bin index `i` is zero-based, as in numpy. |
| `h.GetBinError(i)` | `h.bin_error(i)` | `sqrt(sum of w^2)`, the same convention as ROOT. |
| `h.GetBinCenter(i)` | `h.bin_centers()[i]` | Returns the whole array. |
| `h.GetNbinsX()` | `h.nbins()` | |
| `h.Integral()` | `h.integral()` | Sum of in-range bins. |
| `h.GetEntries()` | `h.entries` | A property, not a method. |
| `h.GetMean()`, `h.GetStdDev()` | `h.mean()`, `h.std_dev()` | `None` for an empty histogram (ROOT returns 0). |
| `h.GetQuantiles(...)` | `h.quantile(p)` | One quantile per call. |
| `h.GetBinContent(i)` for all i | `h.contents()` | numpy array, length `nbins`. |
| `h.GetUnderflow()`, `h.GetOverflow()` | `h.underflow`, `h.overflow` | Properties. |
| `h1.Add(h2)` / `hadd` | `h1.merge(h2)` | Requires identical binning. Associative and commutative. |
| `TH2F(...)` | `td.Histogram2D(nx, xlo, xhi, ny, ylo, yhi)` | Contents are indexed `[ix, iy]`. |
| `h2.ProjectionX()` | `h2.project_x()` | Returns a `Histogram1D`. |
| `TTree` + `Draw` / `Project` | `td.LazyFrame.scan_dict(...)` or `scan_csv(...)`, then `.filter(...).collect()` | Lazy, chunked, with predicate pushdown. |
| `df.Filter("pt > 20")` | `.filter(td.col("pt") > 20.0)` | Expressions are built in Python; there are no string cuts. |
| `df.Define("pt2", "pt*pt")` | `.with_column((td.col("pt") * td.col("pt")).alias("pt2"))` | |
| `TLorentzVector` invariant mass | `td.invariant_mass(pt1, eta1, phi1, m1, pt2, eta2, phi2, m2)` | Columnar: each argument is an expression or a scalar. |
| `ROOT.Math.VectorUtil::DeltaR` | `td.delta_r(eta1, phi1, eta2, phi2)` | |
| `hist.SaveAs("out.root")` | `h.to_dict()` | There is no ROOT file format. Save the dict or the numpy arrays (for example with `np.savez`). Binary histogram serialisation exists in the Rust crate but is not yet exposed to Python. |

## Example 1: a filled histogram with errors

PyROOT:

```python
import ROOT
h = ROOT.TH1F("mll", "m_ll;m [GeV];events", 60, 60.0, 120.0)
for m in masses:
    h.Fill(m)
for i in range(1, h.GetNbinsX() + 1):
    c, e = h.GetBinContent(i), h.GetBinError(i)
    print(h.GetBinCenter(i), c, e)
```

tpt_dynamics:

```python
import numpy as np
import tpt_dynamics as td

h = td.Histogram1D(60, 60.0, 120.0)   # ROOT bins 1..60 are 0..59 here
h.fill(np.asarray(masses))
centers = h.bin_centers()
for c, content, err in zip(centers, h.contents(), h.errors()):
    print(c, content, err)
```

Note the index base: ROOT `GetBinContent(1)` is the first bin, and
`tpt_dynamics` uses zero-based indices. Histograms here have no ROOT-style
"bin 0 is underflow" convention; underflow and overflow are separate
properties.

## Example 2: summing histograms from several files (`hadd`)

PyROOT:

```python
total = ROOT.TH1F("total", "", 50, 0, 200)
for name in ["run1.root", "run2.root"]:
    f = ROOT.TFile.Open(name)
    total.Add(f.Get("pt"))
```

tpt_dynamics: histograms are built per worker, then merged in memory. Histograms
are not yet written to or read from files from Python; if you need to pass them
between processes, send `to_dict()` output.

```python
import numpy as np
import tpt_dynamics as td

parts = []
for chunk in chunks:                       # each chunk is a numpy array
    h = td.Histogram1D(50, 0.0, 200.0)
    h.fill(chunk)
    parts.append(h)

total = td.Histogram1D(50, 0.0, 200.0)
for p in parts:
    total.merge(p)
```

## Example 3: a selection on a tree-like dataset

PyROOT:

```python
tree.Draw("pt >> h(50, 0, 200)", "pt > 20 && abs(eta) < 2.5 && njets >= 2")
```

tpt_dynamics:

```python
import numpy as np
import tpt_dynamics as td
from tpt_dynamics import col

data = {
    "pt": pt,           # numpy float64
    "eta": eta,         # numpy float64
    "njets": njets,     # numpy int64
}
sel = (col("pt") > 20.0) & (col("eta").ge(-2.5)) & (col("eta").le(2.5)) & (col("njets") >= 2)

selected = td.LazyFrame.scan_dict(data).filter(sel).collect()

h = td.Histogram1D(50, 0.0, 200.0)
h.fill(selected["pt"])
```

`abs(eta) < 2.5` is written as two bounds, because the expression language has
no `abs` yet. Only the expression functions listed in the reference
(`invariant_mass`, `delta_r`, `transverse_momentum`) are available as built-ins.

## Example 4: a 2D distribution with its projections

PyROOT:

```python
h2 = ROOT.TH2F("eta_phi", "", 40, -2.5, 2.5, 32, -3.2, 3.2)
for e, p in zip(eta, phi):
    h2.Fill(e, p)
proj = h2.ProjectionX("px")
```

tpt_dynamics:

```python
h2 = td.Histogram2D(40, -2.5, 2.5, 32, -3.2, 3.2)
h2.fill(eta, phi)
proj = h2.project_x()
grid = h2.contents()        # shape (40, 32)
```

## Plotting

`tpt_dynamics` has no `Draw()`. Plot the arrays with matplotlib:

```python
import matplotlib.pyplot as plt

h = td.Histogram1D(50, 0.0, 200.0)
h.fill(pt)
plt.stairs(h.contents(), h.edges(), fill=False)
plt.errorbar(h.bin_centers(), h.contents(), yerr=h.errors(), fmt="none")
plt.xlabel("pT [GeV]")
plt.ylabel("events")
```

## Simulation (no ROOT equivalent)

`run_simulation` runs a molecular dynamics simulation from a TOML string and
returns energies as numpy arrays. See the package README for an example.
