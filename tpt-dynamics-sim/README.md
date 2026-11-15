# tpt-dynamics-sim

Classical molecular dynamics engine: structure-of-arrays system model, PDB and
GRO I/O, TOML run configuration, orthorhombic and triclinic periodic boundaries,
cell and Verlet neighbour lists (rayon-parallel build), Lennard-Jones forces
with cutoff/shift/switch, velocity Verlet and leap-frog behind traits.

Phase 5 additions: topology builder and synthetic-test force-field file format,
harmonic bonds/angles (Urey-Bradley), periodic and improper dihedrals,
exclusions and scaled 1-4 pairs, PME Coulomb, SHAKE / LINCS-style / SETTLE
constraints, Berendsen / v-rescale / Nosé–Hoover thermostats, Berendsen and
C-rescale barostats, steepest descent and L-BFGS minimisation, solvation and
ion placement, and checkpoint/restart.

Additions in this phase: rigid water presets (TIP3P, SPC/E; TIP4P geometry
only, its virtual-site force distribution is not implemented) with SETTLE, an
isotropic Parrinello–Rahman barostat, a streaming trajectory writer that feeds
the `tpt-dynamics-data` frame sink, a reproducibility mode (`deterministic =
true`), bounded line reading for PDB and GRO files (1 MiB lines, 10 million
atoms), and parser fuzz tests.

Not included: AMBER, CHARMM or GROMOS parameter data (only a synthetic test
file ships), anisotropic (tensor) pressure coupling, and validation against
GROMACS output. The Lennard-Jones equation-of-state comparison is recorded in
`docs/validation/lj-fluid.md` and is not yet a passing validation.

Run the SPC/E water box (216 waters, NVT with v-rescale) and write an energy log:

```text
cargo run --release -p tpt-dynamics-sim --example water_box -- energy.csv
```

Run the 10,000-atom Lennard-Jones example in release mode:

```text
cargo run --release -p tpt-dynamics-sim --example lj_fluid
```

Part of [tpt-dynamics](../README.md). Licensed under MIT OR Apache-2.0.
