//! SPC/E water box: 216 rigid waters in an NVT run with stochastic velocity
//! rescaling, PME electrostatics and SETTLE constraints. Writes an energy log.
//!
//! ```text
//! cargo run --release -p tpt-dynamics-sim --example water_box -- <energy.csv>
//! ```
//!
//! The solvation grid places the waters on a 6 × 6 × 6 lattice at roughly
//! liquid density (about 1 g/cm³). The orientations are randomised so the
//! starting configuration has no net polarisation. The run is short (a few ps)
//! and is a stability check, not an equilibrated water simulation.

use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use tpt_dynamics_core::Vec3;
use tpt_dynamics_sim::forcefield::{LjType, Modifier};
use tpt_dynamics_sim::integrator::{ForceEvaluator, VelocityVerlet};
use tpt_dynamics_sim::nonbonded::{
    Electrostatics, Nonbonded, NonbondedSettings, COULOMB_FACTOR_MD,
};
use tpt_dynamics_sim::output::EnergyLogger;
use tpt_dynamics_sim::pbc::SimBox;
use tpt_dynamics_sim::rng::Rng;
use tpt_dynamics_sim::simulation::{maxwell_boltzmann_velocities, Simulation};
use tpt_dynamics_sim::solvate::place_water_grid_geometry;
use tpt_dynamics_sim::system::UnitSystem;
use tpt_dynamics_sim::thermostat::VRescale;
use tpt_dynamics_sim::water::WaterModel;

const TEMPERATURE: f64 = 298.15;
const DT: f64 = 0.002;
const CUTOFF: f64 = 0.8;
const SKIN: f64 = 0.1;
const EQUIL_STEPS: u64 = 2000;
const LOG_EVERY: u64 = 50;
const SEED: u64 = 2024;

/// Rotates the molecule `atoms` (O, H1, H2) by a uniformly random rotation about
/// its oxygen, using Rodrigues' formula with a random axis and angle.
fn randomise_orientation(pos: &mut [Vec3], atoms: [usize; 3], rng: &mut Rng) {
    let axis = {
        let v = Vec3::new(rng.gaussian(), rng.gaussian(), rng.gaussian());
        v * (1.0 / v.norm())
    };
    let angle = std::f64::consts::TAU * rng.uniform();
    let (s, c) = angle.sin_cos();
    let origin = pos[atoms[0]];
    for &a in &atoms[1..] {
        let v = pos[a] - origin;
        let rotated = v * c + axis.cross(v) * s + axis * (axis.dot(v) * (1.0 - c));
        pos[a] = origin + rotated;
    }
}

fn main() {
    let out: PathBuf = std::env::args()
        .nth(1)
        .map_or_else(|| PathBuf::from("water_box_energy.csv"), PathBuf::from);

    let model = WaterModel::SpcE;
    let params = model.params();
    let l = 1.86;
    let mut sys = place_water_grid_geometry(
        SimBox::new(l, l, l).expect("valid box"),
        0.305,
        params.r_oh,
        model.d_hh(),
    )
    .expect("water grid");
    let n_mol = sys.len() / 3;
    println!("waters: {n_mol}, atoms: {}, box: {l} nm", sys.len());

    let mut rng = Rng::seed_from_u64(SEED);
    for m in 0..n_mol {
        randomise_orientation(&mut sys.positions, [3 * m, 3 * m + 1, 3 * m + 2], &mut rng);
    }
    for m in 0..n_mol {
        sys.atom_types[3 * m] = 0;
        sys.charges[3 * m] = params.q_o;
        for h in 1..3 {
            sys.atom_types[3 * m + h] = 1;
            sys.charges[3 * m + h] = params.q_h;
        }
    }
    maxwell_boltzmann_velocities(&mut sys, TEMPERATURE, SEED, UnitSystem::Md);

    // Oxygen carries the Lennard-Jones term; hydrogens have epsilon = 0.
    let settings = NonbondedSettings {
        lj_types: vec![
            LjType { epsilon: params.epsilon_o, sigma: params.sigma_o },
            LjType { epsilon: 0.0, sigma: 1.0 },
        ],
        cutoff: CUTOFF,
        modifier: Modifier::Shift,
        electrostatics: Electrostatics::Pme { tolerance: 1e-4, order: 4 },
        coulomb_factor: COULOMB_FACTOR_MD,
        fudge_lj: 0.5,
        fudge_qq: 0.8,
    };
    // Every pair inside a water molecule is excluded (the geometry is rigid).
    let mut exclusions: Vec<Vec<usize>> = vec![Vec::new(); sys.len()];
    for m in 0..n_mol {
        let (a, b, c) = (3 * m, 3 * m + 1, 3 * m + 2);
        exclusions[a] = vec![b, c];
        exclusions[b] = vec![a, c];
        exclusions[c] = vec![a, b];
    }
    let ff = Nonbonded::new(settings, exclusions, Vec::new()).expect("force field");
    let eval = ForceEvaluator::new(Box::new(ff), SKIN).expect("neighbour list");

    let thermostat =
        VRescale::new(TEMPERATURE, 0.1, UnitSystem::Md, SEED ^ 0x5EED_0001).expect("thermostat");
    let settle = model
        .settle((0..n_mol).map(|m| [3 * m, 3 * m + 1, 3 * m + 2]).collect(), sys.len(), 1e-10)
        .expect("SETTLE");
    let mut sim = Simulation::new(
        sys,
        eval,
        Box::new(VelocityVerlet),
        Some(Box::new(thermostat)),
        DT,
        UnitSystem::Md,
    )
    .expect("simulation setup")
    .with_constraints(Box::new(settle))
    .expect("constraints");

    let file = BufWriter::new(File::create(&out).expect("create energy log"));
    let mut observers: Vec<Box<dyn tpt_dynamics_sim::output::Observer>> =
        vec![Box::new(EnergyLogger::new(file, LOG_EVERY).expect("log header"))];
    let summary = sim.run(EQUIL_STEPS, &mut observers).expect("NVT run");

    let e = summary.last;
    assert!(e.total.is_finite() && e.pressure.is_finite(), "non-finite state");
    println!(
        "NVT v-rescale, {} steps ({:.2} ps) in {:.1} s",
        summary.steps,
        summary.steps as f64 * DT,
        summary.wall_seconds
    );
    println!("T = {:.1} K, U_pot/N_mol = {:.2} kJ/mol", e.temperature, e.potential / n_mol as f64);
    println!("E_total start {:.1}, end {:.1} kJ/mol", summary.initial.total, e.total);
    println!("pressure (kJ/mol/nm^3, constraint virial omitted) = {:.1}", e.pressure);
    println!("neighbour rebuilds: {}, energy log: {}", summary.neighbour_rebuilds, out.display());
}
