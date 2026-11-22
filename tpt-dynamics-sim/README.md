# tpt-dynamics-sim

Classical molecular dynamics engine: structure-of-arrays system model, PDB and
GRO I/O, TOML run configuration, orthorhombic and triclinic periodic boundaries
(GRO and PDB files with triclinic boxes are read and written), cell and Verlet
neighbour lists (rayon-parallel build), Lennard-Jones forces with
cutoff/shift/switch, velocity Verlet and leap-frog behind traits.

Also included: topology builder and synthetic-test force-field file format,
harmonic bonds/angles (Urey-Bradley), periodic and improper dihedrals,
exclusions and scaled 1-4 pairs, PME Coulomb, SHAKE / LINCS / SETTLE
constraints, Berendsen / v-rescale / Nosé–Hoover thermostats, Berendsen,
C-rescale and isotropic Parrinello–Rahman barostats, steepest descent and
L-BFGS minimisation, solvation and ion placement, checkpoint/restart, rigid
water presets (TIP3P, SPC/E; TIP4P geometry only, see below), a streaming
trajectory writer that feeds the `tpt-dynamics-data` frame sink, a
reproducibility mode (`deterministic = true`), bounded line reading for PDB and
GRO files (1 MiB lines, 10 million atoms), and parser fuzz tests.

## Quick start

A run is built from a configuration and a system. An empty system cannot be run, so
real input comes from a PDB or GRO file (or from `System::push_atom`):

```rust
use tpt_dynamics_sim::config::SimConfig;
use tpt_dynamics_sim::pbc::SimBox;
use tpt_dynamics_sim::system::System;

let cfg = SimConfig::from_toml_str("nsteps = 10\ntemperature = 0.7\n").unwrap();
let sys = System::new(SimBox::new(3.0, 3.0, 3.0).unwrap());
// An empty system cannot be run; a real system is built from `pdb` or `gro`
// input or by pushing atoms with `System::push_atom`.
assert!(cfg.build_simulation(sys).is_err());
```

From the command line, the [`tptdyn`](../tpt-dynamics-cli/README.md) tool runs a
configuration against a structure file and writes the energy and trajectory outputs:
`tptdyn validate` checks a configuration without running it, and `tptdyn simulate`
runs it.

## Modules

| Module | Contents |
| --- | --- |
| `system` | Structure-of-arrays particle state: positions, velocities, forces, masses and types. |
| `pbc` | Orthorhombic and triclinic periodic boxes in GROMACS lower-triangular form, with the minimum-image convention. |
| `pdb`, `gro` | Fixed-width structure readers and writers, with bounded line lengths. |
| `config` | The TOML run configuration (`SimConfig`). Unknown keys are rejected. |
| `topology`, `ffparams` | Bond graph to angles, dihedrals, impropers, 1-4 pairs and exclusions; TOML force-field parameter files. The shipped parameter file is synthetic. |
| `bonded` | Harmonic bonds and angles (with Urey–Bradley), periodic dihedrals and impropers. |
| `nonbonded` | Lennard-Jones with exclusions and scaled 1-4 pairs, and PME Coulomb. |
| `neighbour` | Cell lists and Verlet lists with a skin, built in parallel with rayon. |
| `forcefield` | The `ForceField` trait and Lennard-Jones with cutoff, shift or switch. |
| `integrator` | The `Integrator` trait, with velocity Verlet and leap-frog. |
| `constraints` | SHAKE, a LINCS-style matrix expansion, and SETTLE for rigid water. |
| `thermostat`, `barostat` | Berendsen, stochastic v-rescale and Nosé–Hoover thermostats; Berendsen, C-rescale and Parrinello–Rahman barostats. |
| `minimize` | Adaptive steepest descent and L-BFGS. |
| `solvate` | Water-grid placement and ion placement. |
| `checkpoint` | Versioned, CRC-32-protected binary checkpoints. Restart is bitwise. |
| `output` | Observer hook, the binary trajectory format and the CSV energy log. |
| `analysis` | RMSD, radius of gyration, hydrogen bonds and the radial distribution function. |
| `simulation` | The run loop that ties the above together. |

Units follow GROMACS: nm, ps, amu, kJ/mol and K. Reduced Lennard-Jones units
(σ = ε = m = 1) work equally well.

## Run configuration

A run is described by a TOML file. Unknown keys are rejected at every level, and
invalid combinations give a `SimError::Config` message. The top-level keys are
`seed`, `dt`, `nsteps`, `units`, `temperature`, `integrator`, `deterministic`,
and the tables below.

| Table | Options |
| --- | --- |
| `[neighbour]` | `cutoff`, `skin` (nm) |
| `[lennard_jones]` | `epsilon`, `sigma`, `modifier` (`none`, `shift`, `switch`), `switch_start` |
| `[forcefield]` | `kind` = `"lennard-jones"` (default) or `"bonded"`. For `bonded`, `params` is the ffparams file path, resolved against the working directory. `electrostatics` = `"none"` (default) or `"pme"`, with optional `pme_tolerance` and `pme_order`. PME needs `kind = "bonded"`. |
| `[thermostat]` | `kind` = `"none"`, `"berendsen"`, `"v-rescale"`, `"nose-hoover"`; `tau` |
| `[barostat]` | `kind` = `"none"`, `"berendsen"`, `"c-rescale"`, `"parrinello-rahman"`; `pressure`, `compressibility`, `tau` |
| `[constraints]` | `kind` = `"none"`, `"shake"` (`tolerance`, `max_iter`), `"lincs"` (`order`, `iterations`), or `"settle"` (`water` = `"tip3p"` or `"spc-e"`, `molecules` = list of `[O, H1, H2]` indices). SHAKE and LINCS constrain every topology bond at its starting length. |
| `[minimize]` | `enabled`, `method` = `"steepest-descent"` (default) or `"lbfgs"`, `max_iterations`, `force_tolerance`. Runs before the first step. |
| `[checkpoint]` | `restart` (file restored after the build), `write` (file written by `SimConfig::write_checkpoint` after the run) |
| `[output]` | `energy_interval`, `energy_file`, `trajectory_interval`, `trajectory_file`, `trajectory_format` |

Combinations that are rejected include PME with the Lennard-Jones force field,
`minimize.enabled` together with `checkpoint.restart`, a bonded force field
without `params`, and SHAKE or LINCS on a system with no topology bonds.

The `bonded` force field needs a system whose topology and atom types come from
`ParamSet::assemble` or `assemble_chains`. A GRO or PDB structure has neither,
so bonded runs must be built from the Rust API. The configuration is not
run-loop-aware: `write_checkpoint` must be called after the run, and the
`tptdyn simulate` command does not call it.

TIP4P is not supported. Its virtual-site force distribution is not implemented,
so the water builders refuse it.

## Not included

AMBER, CHARMM or GROMOS parameter data (only a synthetic test file ships),
anisotropic (tensor) pressure coupling, and validation against GROMACS output.
The Lennard-Jones equation-of-state comparison is recorded in
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
