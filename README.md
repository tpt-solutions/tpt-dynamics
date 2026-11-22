# tpt-dynamics

A pure-Rust, zero-FFI engine for large-scale physical simulation and data analysis.
It is meant to replace two workflows:

* **ROOT**: analysis of very large particle-physics event data (lazy data frames,
  physics histograms, columnar I/O). Crate: `tpt-dynamics-data`.
* **GROMACS**: classical molecular dynamics of large biomolecular systems (neighbour
  lists, integrators, thermostats, constraints, PME). Crate: `tpt-dynamics-sim`.

## Status

**Version `0.1.0-dev`. Nothing is published yet.** There is no release on GitHub,
no crates.io package and no PyPI package. Interfaces can still change before 0.1.0.
Progress is tracked in [`todo.md`](todo.md), and the design specification is
[`spec.txt`](spec.txt). Release history will go in [`CHANGELOG.md`](CHANGELOG.md).

## Install from source

You need a Rust toolchain. The minimum supported version is 1.95, declared as
`rust-version` in [`Cargo.toml`](Cargo.toml).

```sh
git clone https://github.com/tpt-solutions/tpt-dynamics
cd tpt-dynamics
cargo build --release -p tpt-dynamics-cli
./target/release/tptdyn --version      # Windows: .\target\release\tptdyn.exe
```

Python bindings (PyO3 and maturin, abi3 for Python 3.11 and later). Run this from
a virtual environment:

```sh
python -m pip install maturin numpy pytest
cd bindings/tpt-dynamics-py
maturin develop --release
python -c "import tpt_dynamics; print(tpt_dynamics.version())"
```

Prebuilt binaries: `install.sh` (Linux and macOS) and `install.ps1` (Windows) download
the archive for a GitHub Release, check its SHA-256, and install `tptdyn`. **They work
only after the first release exists.** Until then, build from source. The Homebrew and
winget files under [`packaging/`](packaging/) are drafts with placeholder values.

A Docker image and a devcontainer are also provided: see [`Dockerfile`](Dockerfile)
and [`.devcontainer/`](.devcontainer/).

## Five-minute quickstart

The examples use the `tptdyn` binary built above. Run them from the repository root,
or put `target/release` on your `PATH`.

**1. Scaffold a project.** `init` writes a starter project from a template (`lj-fluid`
or `hep-analysis`) into a directory.

```sh
tptdyn init lj-fluid lj-run          # creates lj-run/simulate.toml and a README
```

**2. Check the configuration.** `validate` checks the file without running anything.
Configuration errors exit with code 3.

```sh
tptdyn validate lj-run/simulate.toml
```

**3. Run it.** `simulate` writes the energy log (`energy.csv`) and the trajectory
(`traj.bin`) where the configuration says. `--report` also writes `report.html`.

```sh
tptdyn simulate lj-run/simulate.toml --report
```

**4. Report on a run.** `report` writes `report.html` for a run directory. It is the
same report that `--report` produces, and you can rerun it on an old run.

```sh
tptdyn report lj-run
```

**5. Analyse events.** `analyze` reads a CSV with a header row. It prints a summary of
each column, keeps the rows that match each `--filter`, and writes one histogram per
`--hist` option. The `hep-analysis` template includes a generator for a sample file.

```sh
python examples/make_events.py -o events.csv   # writes 2000 synthetic rows
tptdyn analyze events.csv --filter 'pt > 20' --hist pt,50,0,100 --out-dir out
```

This writes `out/pt.json`. The `histogram` subcommand fills one histogram from one
column:

```sh
tptdyn histogram events.csv --column pt --bins 50 --range 0,100 --out pt.json
```

Run `tptdyn --help` or `tptdyn <subcommand> --help` for every option. Exit codes are
listed at the end of `--help`: `0` success, `1` runtime error, `2` usage error and
`3` input or configuration error.

## Python quickstart

```python
import numpy as np
from tpt_dynamics import Histogram1D, LazyFrame, col

# Nothing is read until collect(). The filter runs in the lazy engine.
selected = LazyFrame.scan_csv("events.csv").filter(col("pt") > 20.0).collect()

hist = Histogram1D(40, 20.0, 100.0)
hist.fill(np.asarray(selected["pt"], dtype=np.float64))
values, edges = hist.to_numpy()
```

The full script is [`examples/python_histogram.py`](examples/python_histogram.py).

## Crate map

| Crate | Purpose |
|-------|---------|
| `tpt-dynamics-core` | Physical constants, unit conversions, typed quantities and `Vec3`. `no_std`. |
| `tpt-dynamics-data` | Physics histograms, lazy data frames, CSV, columnar and trajectory I/O. The ROOT replacement. |
| `tpt-dynamics-sim` | Molecular dynamics engine: force fields, neighbour lists, integrators, thermostats, constraints. The GROMACS replacement. |
| `tpt-dynamics-compute` | Pair kernels (scalar and SIMD), Ewald real-space terms and PME. |
| `tpt-dynamics-cli` | The `tptdyn` command-line tool (crate `tpt-dynamics-cli`, binary `tptdyn`). |
| `bindings/tpt-dynamics-py` | Python bindings (PyO3), published as the `tpt_dynamics` package. Not published yet. |
| `embedded-demo` | Fixed-point event filter on `thumbv7em-none-eabihf`. A separate workspace, not a member of the root one. |

Checking the `no_std` core for a bare-metal target needs
`rustup target add thumbv7em-none-eabihf` and
`cargo build -p tpt-dynamics-core --no-default-features --target thumbv7em-none-eabihf`.

## Principles

* **Permissive dependencies only.** Every dependency must be MIT, BSD, Zlib, ISC,
  Unicode, or dual `MIT OR Apache-2.0`. `cargo deny` enforces this. See
  [ADR-005](docs/adr/ADR-005-licence-policy.md).
* **No FFI.** No `*-sys` crates, no `bindgen`, no C or C++ wrappers in the normal
  dependency tree. See [ADR-001](docs/adr/ADR-001-pure-rust-md.md).
* **Safe Rust.** Unsafe code is denied by default. Each exception is an ADR and is
  audited. See [ADR-006](docs/adr/ADR-006-unsafe-code-policy.md).

## Documentation

* [`docs/tutorials/`](docs/tutorials/): [Lennard-Jones fluid](docs/tutorials/lj-fluid.md)
  and [dilepton mass analysis](docs/tutorials/hep-analysis.md).
* [`docs/adr/`](docs/adr/): architecture decision records.
* [`docs/benchmarks/`](docs/benchmarks/README.md): performance measurements and the
  regression baseline.
* [`docs/migration/pyroot.md`](docs/migration/pyroot.md): moving analysis code from
  PyROOT to tpt-dynamics.
* [`docs/site/index.html`](docs/site/index.html): a static index of the Markdown
  documentation, generated by `python docs/site/build.py`. A workflow also deploys it
  to GitHub Pages from `master`.
* [`docs/release-checklist.md`](docs/release-checklist.md): the release gates.
* [`examples/`](examples/): sample configurations and scripts.
* [`CONTRIBUTING.md`](CONTRIBUTING.md).

## Licence

Copyright © 2026 TPT Solutions. Dual-licensed under either of

* [MIT licence](LICENSE-MIT)
* [Apache Licence, Version 2.0](LICENSE-APACHE)

at your option. The [NOTICE](NOTICE) file must be carried with any redistribution of
the Apache-2.0 form of this work.

Unless you explicitly state otherwise, any contribution intentionally submitted for
inclusion in this work shall be dual licensed as above, without any additional terms
or conditions.
