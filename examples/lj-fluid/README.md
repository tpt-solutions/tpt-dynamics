# Lennard-Jones fluid example

`lj-64.gro` is a 64-atom simple-cubic Lennard-Jones structure in reduced units
(density 0.8442, T* = 0.72). It was written by

```sh
cargo run --release -p tpt-dynamics-sim --example gen_lj_system -- 4 4 4 examples/lj-fluid/lj-64.gro
```

`lj-small.toml` runs 2000 velocity-Verlet steps and writes an energy log every 100
steps. Run it from this directory so the relative paths resolve:

```sh
cd examples/lj-fluid
tptdyn simulate lj-small.toml --output-dir out
```

The energy log is `out/lj-small-energy.csv`. For the full 10,000-atom benchmark see
`tpt-dynamics-sim/benchmarks/lj-10k-cli.toml`.
