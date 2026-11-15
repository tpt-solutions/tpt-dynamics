# tpt-dynamics

A pure-Rust, zero-FFI engine for large-scale physical simulation and data
analysis. It aims to replace two workflows:

* **ROOT** — analysis of very large particle-physics event data (lazy data frames,
  physics histograms, columnar I/O). Crate: `tpt-dynamics-data`.
* **GROMACS** — classical molecular dynamics of large biomolecular systems
  (neighbour lists, integrators, force fields, domain decomposition). Crate:
  `tpt-dynamics-sim`.

**Status: early development (`0.1.0-dev`).** Nothing is published yet. The
project plan and progress are tracked in [`todo.md`](todo.md); the design
specification is [`spec.txt`](spec.txt).

| Crate | Purpose | Status |
|-------|---------|--------|
| `tpt-dynamics-core` | Constants, unit conversions, typed quantities, `Vec3` (`no_std`) | Usable |
| `tpt-dynamics-data` | 1D physics histograms (more to come) | 1D histograms done |
| `tpt-dynamics-sim` | Molecular dynamics engine | Placeholder |
| `tpt-dynamics-compute` | SIMD / GPU kernels | Placeholder |
| `tpt-dynamics-cli` | `tptdyn` command-line tool | Placeholder (`--help`, `--version`) |
| `bindings/tpt-dynamics-py` | Python bindings (PyO3) | Placeholder |

## Principles

* **Permissive dependencies only.** Every dependency must be MIT, BSD, Zlib, ISC,
  Unicode, or dual `MIT OR Apache-2.0`. GPL, LGPL, AGPL, MPL and Apache-2.0-only
  crates are rejected by `cargo deny`. See
  [ADR-005](docs/adr/ADR-005-licence-policy.md).
* **No FFI.** No `*-sys` crates, no `bindgen`, no C/C++ wrappers.
* **Safe Rust.** `#![forbid(unsafe_code)]` by default. See
  [ADR-006](docs/adr/ADR-006-unsafe-code-policy.md).
* **`no_std` core.** `tpt-dynamics-core` builds for bare-metal targets (for example
  `thumbv7em-none-eabihf`) so it can run on FPGA controllers and DAQ hardware.

## Building

Requirements: a Rust toolchain. The version is pinned by
[`rust-toolchain.toml`](rust-toolchain.toml); the minimum supported version is
declared as `rust-version` in [`Cargo.toml`](Cargo.toml).

```sh
cargo build --workspace
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all -- --check
cargo deny check

# no_std check (requires the target to be installed)
rustup target add thumbv7em-none-eabihf
cargo build -p tpt-dynamics-core --no-default-features --target thumbv7em-none-eabihf
```

## Example

```rust
use tpt_dynamics_core::quantity::Energy;
use tpt_dynamics_data::Histogram1D;

let barrier = Energy::from_kcal_per_mol(3.2);
println!("barrier = {:.3} eV", barrier.ev());

let mut h = Histogram1D::uniform(50, 0.0, 10.0)?;
h.fill(2.5);
assert_eq!(h.integral(), 1.0);
# Ok::<(), tpt_dynamics_data::HistogramError>(())
```

## Documentation

* [`docs/adr/`](docs/adr/) — architecture decision records
* [`docs/dependencies.md`](docs/dependencies.md) — dependency licence audit
* [`CONTRIBUTING.md`](CONTRIBUTING.md) — how to contribute
* [`CHANGELOG.md`](CHANGELOG.md) — release history

## Licence

Dual-licensed under either of

* [MIT licence](LICENSE-MIT)
* [Apache Licence, Version 2.0](LICENSE-APACHE)

at your option.

Copyright © 2026 TPT Solutions.

The project [NOTICE](NOTICE) file must be carried with any redistribution of the
Apache-2.0 licensed form of this work, as the Apache-2.0 licence requires.

### Contributions

Unless you explicitly state otherwise, any contribution intentionally submitted
for inclusion in this work by you, as defined in the Apache-2.0 licence, shall be
dual licensed as above, without any additional terms or conditions.
