# tpt-dynamics-cli

The `tptdyn` command-line tool. Part of [tpt-dynamics](../README.md).
Licensed under MIT OR Apache-2.0.

## Global options

| Option | Meaning |
| --- | --- |
| `--log-level LEVEL` | `off`, `error`, `warn` (default), `info`, `debug`, `trace`. Diagnostics go to stderr. |
| `-q, --quiet` | Suppress progress lines and the end-of-run summary. |
| `--config FILE` | CLI defaults (TOML). Keys: `log_level`, `quiet`, `threads`, `output_dir`. Unknown keys are rejected. Flags override the file. |
| `--threads N` | Worker threads for rayon (default: all logical cores). |
| `-h, --help`, `-V, --version` | Help and version. |

## Subcommands

### `tptdyn init TEMPLATE [DIR] [--force]`

Writes a runnable project from a built-in template. Templates:

- `lj-fluid`: a Lennard-Jones fluid on a 4 x 4 x 4 cubic lattice. Writes
  `simulate.toml` (energy log every 10 steps, trajectory every 50) and a README.
- `hep-analysis`: `make_events.py`, which writes a synthetic event CSV using only
  the Python standard library, plus `analyze.txt`, a command sheet of `analyze`
  and `histogram` examples.

`DIR` defaults to the current directory and is created if missing. No file is
overwritten unless `--force` is given. If any target exists, nothing is written
and the command exits with code 3.

### `tptdyn validate CONFIG.toml`

Parses a simulation configuration the way `simulate` does, builds the starting
system and engine, and prints a summary without stepping: atom count, box, number
of steps, expected energy-log rows and an estimated trajectory size in MB. Exit
code 0 if the configuration is valid, 3 if not.

Errors name the file. An unknown key also gets a suggestion when a close key
exists (edit distance at most 2):

```text
tptdyn: error: config 'run.toml': invalid configuration: TOML parse error at line 17, column 1
...
unknown field `energy_intervall`, expected one of `energy_interval`, ...
  hint: did you mean `energy_interval`?
```

### `tptdyn simulate CONFIG.toml [--output-dir DIR] [--report]`

Runs the MD engine described by a TOML file. The file is the engine's
`SimConfig` plus a `[system]` table naming the starting structure:

```toml
seed = 7
dt = 0.002
nsteps = 50
temperature = 0.7

[system]
structure = "start.gro"   # .gro or .pdb, or:
# cubic_cells = 3         # cubic_cells^3 identical atoms on a lattice,
# spacing = 1.6           # with spacing (nm) between sites,
# mass = 1.0              # and mass (amu, default 1.0)

[neighbour]
cutoff = 2.0
skin = 0.3

[output]
energy_interval = 10
energy_file = "energy.csv"          # CSV log
trajectory_interval = 10
trajectory_file = "traj.bin"        # binary trajectory
```

Relative paths in the file resolve against the directory containing it.
`--output-dir` writes the energy and trajectory files into that directory
under their configured names.

Progress: when stderr is a terminal and `--quiet` is not given, one line is
redrawn at most four times a second:

```text
step 500/1000 ( 50.0%)  ETA 10s       8.64 ns/day
```

ETA is the remaining steps at the average rate so far. ns/day is simulated time
per wall-clock day. When stderr is not a terminal (a file or a pipe), a plain
`step i/n (p%)` line is printed about every tenth of the run instead.

`--report` writes `report.html` (see `tptdyn report`) next to the energy log
after the run. It needs `energy_file` in `[output]`.

### `tptdyn report RUN_DIR [--trajectory FILE] [--out FILE]`

Reads `RUN_DIR/energy.csv`, the log written by `simulate` with columns
`step,time,kinetic,potential,total,temperature,pressure` (columns are found by
name), and writes `RUN_DIR/report.html`, or the file given with `--out`.

The page is self-contained. It has inline SVG plots of energy drift, temperature
and pressure against time, and a summary table: mean and sample standard
deviation of T and P, the energy drift in percent of the initial total energy
(final value and largest excursion), and, with `--trajectory`, the frame and atom
count of a plain trajectory. A relative `--trajectory` path is taken from
`RUN_DIR`. The page has no scripts, no external stylesheets and no network access.

### `tptdyn analyze INPUT [options]`

Reads an event table and prints the row count and, for each column, the min, max,
mean and standard deviation. Every column must be numeric: integer columns become
`i64`, others `f64`. The input format is detected from the file's first bytes (see
[Input formats](#input-formats)).

- `--filter 'COLUMN OP VALUE'`: keep rows where the comparison holds.
  `OP` is one of `>`, `>=`, `<`, `<=`, `==`, `!=`. Repeat for AND.
- `--where 'EXPR'`: keep rows where the expression holds (see
  [Expression syntax](#expression-syntax)). Repeat for AND.
- `--hist COLUMN,BINS,LO,HI`: uniform histogram of a (filtered) column.
  Repeatable.
- `--out-dir DIR` (default `.`) and `--format json|bin` (default `json`):
  histogram files are written as `COLUMN.json` or `COLUMN.bin`.
- `--delimiter CHAR` and `--no-header` (CSV only).
- `--trajectory FILE` analyses a trajectory instead (see below).

Filters and expressions run through the lazy-frame engine of `tpt-dynamics-data`.
CSV is processed in chunks. `--filter` and `--where` can be mixed; every condition
must hold.

#### Expression syntax

`--where` takes one expression over column names and numbers:

| Element | Examples |
| --- | --- |
| Numbers | `2`, `0.5`, `.5`, `1e-3` (arithmetic is `f64`) |
| Columns | identifiers: letters, digits and `_`, not starting with a digit |
| Arithmetic | `+ - * /` and unary `-`; `*` and `/` bind tighter than `+` and `-` |
| Comparisons | `< <= > >= == !=`; they do not chain (`a < b < c` is an error) |
| Logic | `not`, `and`, `or` (lowercase). Precedence: `not` over `and` over `or` |
| Grouping | parentheses, up to 128 levels deep |

Examples:

```text
--where 'pt > 20 and eta < 2.5 and eta > -2.5'
--where '(pt > 20 or n_jets >= 2) and eta < 2.5'
--where 'not (x < 0) and y * 0.5 > 1.2'
```

Division by zero is not an error. It follows IEEE 754: `x / 0` gives `±inf` and
`0 / 0` gives NaN. Comparisons with NaN are false, so `<`, `<=`, `>`, `>=` and `==`
drop such rows, while `!=` keeps them. A syntax error exits with code 2 and reports
the column. An unknown column exits with code 3.

#### Input formats

| Format | Recognised by | Notes |
| --- | --- | --- |
| CSV | default, when no magic bytes match | streamed in chunks; `--delimiter`, `--no-header` |
| Parquet | leading `PAR1` | `tpt-dynamics-data` `parquet` feature: flat columns, SNAPPY and GZIP codecs |
| Columnar | leading `TPTC` | `tpt-dynamics-data` columnar table; `u32` columns are read as `i64` |
| Row store | leading `TPTROWS1` | `tpt-dynamics-data` row store, read through its lazy source |

The binary formats are loaded whole, as one table. A trajectory file (`TPTTRJ01`
or `TPTCTRJ1`) is rejected, with a pointer to `--trajectory`.

### `tptdyn analyze --trajectory FILE`

Analyses a trajectory and writes one CSV row per frame: RMSD against the first
frame and radius of gyration. Plain (`TPTTRJ01`) and chunked (`TPTCTRJ1`)
trajectories are both accepted. Use `--out FILE` for the CSV and `--structure` for
masses from a `.gro` or `.pdb` file.

### `tptdyn histogram INPUT.csv --column NAME|INDEX --bins N --range LO,HI --out FILE`

Streams one column into a uniform histogram without loading the whole file.
`--column` is a name when the file has a header, otherwise a zero-based index.
`--format json|bin` selects the output format. By default a `.bin` file name
gives binary and anything else gives JSON. Both formats are read back with
`tpt_dynamics_data::serialise`.

### `tptdyn completions SHELL`

Writes a completion script for `bash`, `zsh`, `fish`, `powershell` or `elvish`
to stdout, for example `tptdyn completions bash > tptdyn.bash`.

### `tptdyn man [--dir DIR]`

Writes `tptdyn.1` and `tptdyn-<subcommand>.1` into `DIR`. Without `--dir`, all
pages are written to stdout.

## Exit codes

| Code | Meaning |
| --- | --- |
| 0 | Success (including `--help` and `--version`). |
| 1 | Runtime error: a simulation step failed, or output could not be written. |
| 2 | Usage error: unknown option, missing or malformed argument, bad filter, expression or range. |
| 3 | Input or configuration error: missing or unreadable file, malformed CSV or TOML, unknown column, existing file without `--force`, invalid simulation setup. |

Errors are printed to stderr as `tptdyn: error: <message>`.

## Testing

`cargo test -p tpt-dynamics-cli` runs the unit tests (expression parser, progress
formatting, report, input format detection, init templates) and the integration
tests in `tests/cli.rs`. The integration tests cover help and version, exit codes,
histogram bin counts, analyze filters and expressions, binary input formats, init,
validate, report, the `--report` flag, completions, man pages and short
simulations.
