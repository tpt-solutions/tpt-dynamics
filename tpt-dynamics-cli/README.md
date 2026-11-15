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

### `tptdyn simulate CONFIG.toml [--output-dir DIR]`

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
under their configured names. Progress (`step i/n (p%)`) goes to stderr.

### `tptdyn analyze INPUT.csv [options]`

Loads a CSV file with a header (or `--no-header`). Every column must be
numeric: integer columns become `i64`, others `f64`. Prints the row count and
per-column min, max, mean and standard deviation to stdout.

- `--filter 'COLUMN OP VALUE'`: keep rows where the comparison holds.
  `OP` is one of `>`, `>=`, `<`, `<=`, `==`, `!=`. Repeat for AND.
  Filters run through the lazy-frame engine. There is no general expression
  syntax yet.
- `--hist COLUMN,BINS,LO,HI`: uniform histogram of a (filtered) column.
  Repeatable.
- `--out-dir DIR` (default `.`) and `--format json|bin` (default `json`):
  histogram files are written as `COLUMN.json` or `COLUMN.bin`.
- `--delimiter CHAR`, `--no-header`.

### `tptdyn histogram INPUT.csv --column NAME|INDEX --bins N --range LO,HI --out FILE`

Streams one column into a uniform histogram without loading the whole file.
`--column` is a name when the file has a header, otherwise a zero-based index.
`--format json|bin` selects the output format; by default a `.bin` file name
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
| 2 | Usage error: unknown option, missing or malformed argument, bad filter or range. |
| 3 | Input or configuration error: missing or unreadable file, malformed CSV or TOML, unknown column, invalid simulation setup. |

Errors are printed to stderr as `tptdyn: error: <message>`.

## Testing

`cargo test -p tpt-dynamics-cli` runs the integration tests in `tests/cli.rs`,
which cover help and version, exit codes, histogram bin counts, analyze
filters, completions, man pages and a 27-atom, 50-step simulation.
