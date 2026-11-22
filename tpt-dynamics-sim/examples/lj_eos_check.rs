//! Lennard-Jones equation-of-state check at T* = 1.0, ρ* = 0.8.
//!
//! The potential is the full 12-6 form truncated at 3σ, with the standard
//! analytic long-range corrections for energy and pressure (the convention of
//! the NIST TMMC tables, "3σ + standard long range corrections"). The
//! simulation is a 864-atom FCC crystal melted at T* = 1.0 and sampled with a
//! stochastic velocity-rescaling thermostat (NVT). Pressure and energy are
//! block-averaged, and the block standard error is reported.
//!
//! ```text
//! cargo run --release -p tpt-dynamics-sim --example lj_eos_check
//! ```

use std::time::Instant;

use tpt_dynamics_core::Vec3;
use tpt_dynamics_sim::forcefield::{LennardJones, Modifier};
use tpt_dynamics_sim::integrator::{ForceEvaluator, VelocityVerlet};
use tpt_dynamics_sim::pbc::SimBox;
use tpt_dynamics_sim::simulation::{maxwell_boltzmann_velocities, Simulation};
use tpt_dynamics_sim::system::{NewAtom, System, UnitSystem};
use tpt_dynamics_sim::thermostat::VRescale;

const RHO: f64 = 0.8;
const T: f64 = 1.0;
const RC: f64 = 3.0;
const CELLS: usize = 6; // 4 * 6^3 = 864 atoms
const DT: f64 = 0.002;
const EQUIL: u64 = 5000;
const DEFAULT_BLOCKS: u64 = 10;
const STEPS_PER_BLOCK: u64 = 1000;
const SAMPLE_EVERY: u64 = 10;

/// Analytic long-range corrections for the truncated (unshifted) LJ potential.
/// Returns `(u_tail per particle, p_tail)` in reduced units.
fn tail_corrections(rho: f64, rc: f64) -> (f64, f64) {
    let r3 = 1.0 / (rc * rc * rc);
    let r9 = r3 * r3 * r3;
    let u = (8.0 * std::f64::consts::PI / 3.0) * rho * (r9 / 3.0 - r3);
    let p = (16.0 * std::f64::consts::PI / 3.0) * rho * rho * (2.0 * r9 / 3.0 - r3);
    (u, p)
}

fn build_fcc() -> System {
    let a = (4.0 / RHO).cbrt();
    let l = a * CELLS as f64;
    let mut sys = System::new(SimBox::new(l, l, l).expect("valid box"));
    let basis = [
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(0.5, 0.5, 0.0),
        Vec3::new(0.5, 0.0, 0.5),
        Vec3::new(0.0, 0.5, 0.5),
    ];
    for ix in 0..CELLS {
        for iy in 0..CELLS {
            for iz in 0..CELLS {
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
    // Optional: `--blocks N` (production blocks, default 10) and `--seed-offset K`
    // (added to the velocity and thermostat seeds, default 0, which reproduces the
    // committed run). Use several offsets for independent replicas.
    let args: Vec<String> = std::env::args().collect();
    let arg_value = |flag: &str| args.iter().position(|a| a == flag).and_then(|i| args.get(i + 1));
    let blocks: u64 =
        arg_value("--blocks").map_or(DEFAULT_BLOCKS, |v| v.parse().expect("--blocks N"));
    let offset: u64 = arg_value("--seed-offset").map_or(0, |v| v.parse().expect("--seed-offset K"));

    let mut sys = build_fcc();
    let n = sys.len();
    maxwell_boltzmann_velocities(&mut sys, T, 42 + offset, UnitSystem::Reduced);
    let ff = LennardJones::single_type(1.0, 1.0, RC, Modifier::None).expect("LJ");
    let eval = ForceEvaluator::new(Box::new(ff), 0.3).expect("neighbour list");
    let thermostat = VRescale::new(T, 0.5, UnitSystem::Reduced, 7 + offset).expect("thermostat");
    let mut sim = Simulation::new(
        sys,
        eval,
        Box::new(VelocityVerlet),
        Some(Box::new(thermostat)),
        DT,
        UnitSystem::Reduced,
    )
    .expect("simulation setup");

    let started = Instant::now();
    sim.run(EQUIL, &mut []).expect("equilibration");

    let (u_tail, p_tail) = tail_corrections(RHO, RC);
    let mut block_p = Vec::new();
    let mut block_u = Vec::new();
    for _ in 0..blocks {
        let (mut p_sum, mut u_sum, mut count) = (0.0, 0.0, 0u64);
        for s in 0..STEPS_PER_BLOCK {
            sim.advance().expect("step");
            if s % SAMPLE_EVERY == 0 {
                let e = sim.energies();
                p_sum += e.pressure;
                u_sum += e.potential / n as f64;
                count += 1;
            }
        }
        block_p.push(p_sum / count as f64 + p_tail);
        block_u.push(u_sum / count as f64 + u_tail);
    }
    let elapsed = started.elapsed().as_secs_f64();

    let (p_mean, p_err) = mean_and_sem(&block_p);
    let (u_mean, u_err) = mean_and_sem(&block_u);
    println!("LJ fluid, N = {n}, rho* = {RHO}, T* = {T}, rc = {RC} sigma (truncated + LRC)");
    println!(
        "run: {EQUIL} equilibration + {} production steps in {elapsed:.1} s",
        blocks * STEPS_PER_BLOCK
    );
    println!("P*  = {p_mean:.4} +/- {p_err:.4}  (block SEM, {blocks} blocks)");
    println!("U*/N = {u_mean:.4} +/- {u_err:.4}");
    println!("tail corrections: U_tail/N = {u_tail:.4}, P_tail = {p_tail:.4}");
    assert!(p_mean.is_finite() && u_mean.is_finite(), "non-finite result");
}

fn mean_and_sem(x: &[f64]) -> (f64, f64) {
    let n = x.len() as f64;
    let mean = x.iter().sum::<f64>() / n;
    let var = x.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / (n - 1.0);
    (mean, (var / n).sqrt())
}
