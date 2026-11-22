# Changelog

All notable changes to `tpt-dynamics-sim` are documented in this file. The workspace
changelog is in the [root `CHANGELOG.md`](../CHANGELOG.md).

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the
crate uses [Semantic Versioning](https://semver.org/) once it reaches 0.1.0.

## [Unreleased]

### Added

- `system`: structure-of-arrays particle state with masses, types and topology.
- `pbc`: orthorhombic and triclinic periodic boxes in GROMACS lower-triangular form,
  with the minimum-image convention.
- `pdb`, `gro`: fixed-width structure readers and writers with triclinic-box support.
  Line length and atom count are bounded (1 MiB lines, 10 million atoms).
- `config`: TOML run configuration (`SimConfig::from_toml_str`). Unknown keys are
  rejected at every level, and invalid combinations give a `SimError::Config`.
- `topology`, `ffparams`: bond graph to angles, dihedrals, impropers, 1-4 pairs and
  exclusions; TOML force-field parameter files. The shipped parameter file is
  synthetic.
- `bonded`: harmonic bonds and angles with Urey–Bradley terms, periodic dihedrals and
  impropers.
- `nonbonded`: Lennard-Jones with exclusions and scaled 1-4 pairs, and PME Coulomb
  through `tpt-dynamics-compute`.
- `neighbour`: cell lists and Verlet lists with a skin, built in parallel with rayon.
- `forcefield`: the `ForceField` trait and Lennard-Jones with cutoff, shift or switch
  modifiers.
- `integrator`: the `Integrator` trait, with velocity Verlet and leap-frog.
- `constraints`: SHAKE, a LINCS-style matrix expansion, and SETTLE for rigid water.
- `thermostat`: Berendsen, stochastic v-rescale and Nosé–Hoover (chain length 1).
- `barostat`: isotropic Berendsen, C-rescale and Parrinello–Rahman pressure coupling.
- `minimize`: adaptive steepest descent and L-BFGS.
- `solvate`: water-grid placement, ion placement and rigid-water presets (TIP3P, SPC/E).
- `checkpoint`: versioned, CRC-32-protected binary checkpoints, with bitwise restart.
- `output`: observer hook, binary trajectory format, streaming trajectory writer and
  CSV energy log.
- `analysis`: RMSD, radius of gyration, hydrogen bonds and the radial distribution
  function.
- A reproducibility mode (`deterministic = true`) and a deterministic PRNG.
- Parser fuzz tests for PDB and GRO input.
- Examples: `water_box`, `lj_fluid`, `lj_eos_check`, `lj_scaling`, `gen_lj_system`
  and `ffenergy`.

### Known limitations

- TIP4P is not supported. Its virtual-site force distribution is not implemented, so
  the water builders refuse it.
- The `bonded` force field needs a system built through `ParamSet::assemble` or
  `assemble_chains`. GRO and PDB structures carry neither topology nor types, so
  bonded runs are built from the Rust API.
- `write_checkpoint` must be called after the run. `tptdyn simulate` does not call it.
- AMBER, CHARMM and GROMOS parameter data are not included, and anisotropic pressure
  coupling is not implemented.
- Validation against GROMACS output has not been done. The Lennard-Jones equation of
  state comparison is recorded in `docs/validation/lj-fluid.md` and is not yet a
  passing validation.
