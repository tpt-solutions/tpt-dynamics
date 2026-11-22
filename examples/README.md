# Examples

Small, self-contained inputs and scripts. None of them need the full benchmark
data. Scripts run from any directory. Output files are written where you run them.

| Path | What it shows | How to run |
| --- | --- | --- |
| [`lj-fluid/lj-small.toml`](lj-fluid/) | A 64-atom Lennard-Jones run for `tptdyn simulate`, with its structure file `lj-64.gro` | `cd examples/lj-fluid && tptdyn simulate lj-small.toml --output-dir out` |
| [`tptdyn init water-box`](../tpt-dynamics-cli/templates/water-box/) | 512 rigid three-site water molecules in a periodic box, with PME and SETTLE. The force field is not bundled: supply your own `water.toml` (see the template's README) | `tptdyn init water-box my-water && cd my-water`, then add `water.toml` and run `tptdyn validate simulate.toml` |
| [`make_events.py`](make_events.py) | Writes a synthetic event CSV (`event,pt,eta,phi,weight`). Standard library only | `python examples/make_events.py -o events.csv` |
| [`python_histogram.py`](python_histogram.py) | Fills a `Histogram1D` from a lazy-frame selection with the `tpt_dynamics` Python API | `python examples/python_histogram.py events.csv` |

## Command-line analysis

```sh
python examples/make_events.py -o events.csv

# Keep pt > 20 and write pt.json (a histogram file) into ./out
tptdyn analyze events.csv --filter 'pt > 20' --hist pt,50,0,100 --out-dir out

# Or fill one histogram directly from one column
tptdyn histogram events.csv --column pt --bins 50 --range 0,100 --out pt.json
```

`analyze` prints a per-column summary (count, min, max, mean, standard deviation)
before it writes the histograms.

## Python

The Python example needs the `tpt_dynamics` package. Build it from the bindings
directory first (see the root README) and install `numpy`.

```sh
python examples/python_histogram.py events.csv --pt-min 20 --bins 40 --lo 20 --hi 100
```

The generated CSV and the simulation output are ignored by git (`*.csv` in
`.gitignore`), so generate them locally.
