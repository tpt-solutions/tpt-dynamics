//! Milestone check: a stable 10,000-atom Lennard-Jones fluid.
//!
//! Builds an FCC crystal of 10,000 atoms (10 × 10 × 25 conventional cells) at
//! reduced density 0.8442 and T* = 0.72, equilibrates with a Berendsen
//! thermostat, then runs NVE and reports the energy drift and timing.
//!
//! Run in release mode for meaningful timings:
//!
//! ```text
//! cargo run --release -p tpt-dynamics-sim --example lj_fluid
//! ```

use std::time::Instant;

use tpt_dynamics_core::Vec3;
use tpt_dynamics_sim::forcefield::{LennardJones, Modifier};
use tpt_dynamics_sim::integrator::{ForceEvaluator, VelocityVerlet};
use tpt_dynamics_sim::pbc::SimBox;
use tpt_dynamics_sim::simulation::{maxwell_boltzmann_velocities, Simulation};
use tpt_dynamics_sim::system::{NewAtom, System, UnitSystem};
use tpt_dynamics_sim::thermostat::Berendsen;

const DENSITY: f64 = 0.8442;
const TEMPERATURE: f64 = 0.72;
const CELLS: [usize; 3] = [10, 10, 25];
const DT: f64 = 0.004;
const CUTOFF: f64 = 2.5;
const SKIN: f64 = 0.3;

fn build_fcc() -> System {
    let n_atoms = 4 * CELLS.iter().product::<usize>();
    let volume = n_atoms as f64 / DENSITY;
    // Lattice constant for the cubic cell at this density; the box is
    // rectangular so that the atom count is exactly 10,000.
    let a = (4.0 / DENSITY).cbrt();
    let lengths = [CELLS[0] as f64 * a, CELLS[1] as f64 * a, CELLS[2] as f64 * a];
    debug_assert!((lengths.iter().product::<f64>() - volume).abs() < 1e-6 * volume);
    let mut sys = System::new(SimBox::new(lengths[0], lengths[1], lengths[2]).expect("valid box"));
    let basis = [
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(0.5, 0.5, 0.0),
        Vec3::new(0.5, 0.0, 0.5),
        Vec3::new(0.0, 0.5, 0.5),
    ];
    for ix in 0..CELLS[0] {
        for iy in 0..CELLS[1] {
            for iz in 0..CELLS[2] {
                for b in basis {
                    let p = Vec3::new(ix as f64 + b.x, iy as f64 + b.y, iz as f64 + b.z) * a;
                    sys.push_atom(NewAtom {
                        name: "AR".into(),
                        residue_name: "LJ".into(),
                        residue_id: 1,
                        element: "Ar".into(),
                        position: p,
                        velocity: Vec3::ZERO,
                        mass: 1.0,
                        charge: 0.0,
                        atom_type: 0,
                    });
                }
            }
        }
    }
    sys
}

fn main() {
    let mut system = build_fcc();
    maxwell_boltzmann_velocities(&mut system, TEMPERATURE, 12345, UnitSystem::Reduced);
    let n = system.len();
    println!("atoms: {n}, box: {:?}, density: {DENSITY}", system.sim_box.lengths);

    let ff = LennardJones::single_type(1.0, 1.0, CUTOFF, Modifier::Shift).expect("valid LJ");
    let eval = ForceEvaluator::new(Box::new(ff), SKIN).expect("valid neighbour list");
    let thermostat =
        Berendsen::new(TEMPERATURE, 0.5, UnitSystem::Reduced).expect("valid thermostat");

    let t_setup = Instant::now();
    let mut sim = Simulation::new(
        system,
        eval,
        Box::new(VelocityVerlet),
        Some(Box::new(thermostat)),
        DT,
        UnitSystem::Reduced,
    )
    .expect("simulation setup");
    println!("setup: {:.3} s", t_setup.elapsed().as_secs_f64());

    // Phase 1: equilibrate with the thermostat.
    let equil = sim.run(300, &mut []).expect("equilibration");
    println!(
        "equilibration: {} steps in {:.2} s, T = {:.4}, U/N = {:.4}",
        equil.steps,
        equil.wall_seconds,
        equil.last.temperature,
        equil.last.potential / n as f64
    );

    // Phase 2: NVE production run; check energy conservation.
    sim.thermostat = None;
    let nve = sim.run(300, &mut []).expect("NVE run");
    assert!(nve.last.total.is_finite() && nve.last.pressure.is_finite(), "non-finite state");
    assert!(sim
        .system
        .positions
        .iter()
        .all(|p| p.x.is_finite() && p.y.is_finite() && p.z.is_finite()));
    println!(
        "NVE: {} steps in {:.2} s ({:.1} ms/step), T = {:.4}, P = {:.4}",
        nve.steps,
        nve.wall_seconds,
        1e3 * nve.wall_seconds / nve.steps as f64,
        nve.last.temperature,
        nve.last.pressure
    );
    println!(
        "NVE energy: E0 = {:.4}, E_end = {:.4}, relative drift = {:.3e}",
        nve.initial.total,
        nve.last.total,
        nve.relative_energy_drift()
    );
    println!(
        "neighbour list rebuilds during NVE: {}, total momentum: {:.2e}",
        nve.neighbour_rebuilds,
        sim.system.total_momentum().norm()
    );
}
