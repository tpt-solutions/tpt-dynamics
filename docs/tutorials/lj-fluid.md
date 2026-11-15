# Tutorial: Lennard-Jones fluid with the simulation crate

This tutorial explains the `lj_fluid` example in
[`tpt-dynamics-sim/examples/lj_fluid.rs`](../../tpt-dynamics-sim/examples/lj_fluid.rs).
It builds a 10,000-atom face-centred cubic (FCC) crystal in reduced units, equilibrates
it with a thermostat, then checks energy conservation in a microcanonical (NVE) run.
This is the Phase 2 milestone check.

Run it in release mode. Debug builds are much slower:

```text
cargo run --release -p tpt-dynamics-sim --example lj_fluid
```

## 1. Build the system

`build_fcc()` places four atoms per conventional cell on an FCC basis. The box is
rectangular, `10 x 10 x 25` cells, so the atom count is exactly
`4 * 10 * 10 * 25 = 10,000`. The lattice constant comes from the reduced density
`rho* = 0.8442`:

```rust
let a = (4.0 / DENSITY).cbrt();
let lengths = [CELLS[0] as f64 * a, CELLS[1] as f64 * a, CELLS[2] as f64 * a];
let mut sys = System::new(SimBox::new(lengths[0], lengths[1], lengths[2]).expect("valid box"));
```

Each atom is added with `push_atom(NewAtom { .. })`. In reduced units the mass is 1
and the charge is 0. Velocities start at zero and are set in the next step.

## 2. Initial velocities

```rust
maxwell_boltzmann_velocities(&mut system, TEMPERATURE, 12345, UnitSystem::Reduced);
```

The seed (`12345`) makes the run reproducible. The target temperature is `T* = 0.72`.

## 3. Force field, neighbour list and thermostat

```rust
let ff = LennardJones::single_type(1.0, 1.0, CUTOFF, Modifier::Shift).expect("valid LJ");
let eval = ForceEvaluator::new(Box::new(ff), SKIN).expect("valid neighbour list");
let thermostat = Berendsen::new(TEMPERATURE, 0.5, UnitSystem::Reduced).expect("valid thermostat");
```

- `Modifier::Shift` shifts the potential so it is zero at the cutoff `r_c = 2.5 sigma`,
  which removes the energy discontinuity there. `Modifier::None` gives plain truncation,
  and `Modifier::Switch` gives a smooth switch.
- `ForceEvaluator` owns the neighbour list. `SKIN = 0.3` is the Verlet skin, so the list
  is rebuilt only after atoms have moved far enough to matter.
- `Berendsen` couples the system to a heat bath with relaxation time `0.5`.

The force field, integrator and thermostat are passed as trait objects
(`Box<dyn ...>`), so they can be swapped at run time (see ADR-001).

## 4. Equilibrate, then run NVE

```rust
let mut sim = Simulation::new(
    system,
    eval,
    Box::new(VelocityVerlet),
    Some(Box::new(thermostat)),
    DT,
    UnitSystem::Reduced,
)
.expect("simulation setup");
let equil = sim.run(300, &mut []).expect("equilibration"); // thermostat on
sim.thermostat = None;
let nve = sim.run(300, &mut []).expect("NVE run");          // microcanonical
```

`sim.run` returns a summary with the `initial` and `last` states, the wall time,
`relative_energy_drift()` and the number of neighbour-list rebuilds. In a correct NVE
run the drift stays small. The example uses `DT = 0.004`.

## 5. What the example checks

The example asserts that the final energies and all positions are finite. It then
prints:

- the mean temperature and potential energy per atom after equilibration;
- the NVE relative energy drift, which should be small;
- the total momentum, which should stay near zero;
- milliseconds per step, the number to compare when tuning performance.

## Changing the run

- Density and temperature are the constants `DENSITY` and `TEMPERATURE`. Keep
  `rho* = 0.8442` and `T* = 0.72` to reproduce the liquid state used in validation.
- `CELLS` sets the system size. `[4, 4, 4]` gives 256 atoms, which runs in seconds and
  is useful for experiments.
- Setting the modifier to `Modifier::None` shows the energy error that an unshifted
  cutoff introduces.

## Current limits

- The example uses velocity Verlet. The crate also has `LeapFrog`, which the example
  does not use.
- The only thermostat is Berendsen. Nose-Hoover and v-rescale are Phase 5 items.
- The example is single-threaded. Parallel neighbour-list builds and domain
  decomposition are later phases.
- `TrajectoryWriter` and `EnergyLogger` exist in `tpt_dynamics_sim::output`. The example
  does not attach them.
