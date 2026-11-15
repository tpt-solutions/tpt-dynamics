//! Reference validation of the Lennard-Jones engine against exact low-density
//! theory.
//!
//! No published tabulated value was verified for this test, so the reference is
//! derived here: for a dilute fluid the pressure obeys the virial expansion
//! `P = ρ k_B T (1 + B₂ ρ + O(ρ²))`, where the second virial coefficient
//!
//! `B₂(T) = −2π ∫₀^∞ (exp(−U(r)/T) − 1) r² dr`
//!
//! is computed here by quadrature for the same truncated-and-shifted potential
//! the engine uses. The simulated time-averaged pressure must match this within
//! a few per cent at a density where the `ρ²` term is negligible.

use tpt_dynamics_core::Vec3;
use tpt_dynamics_sim::forcefield::{LennardJones, Modifier};
use tpt_dynamics_sim::integrator::{ForceEvaluator, VelocityVerlet};
use tpt_dynamics_sim::pbc::SimBox;
use tpt_dynamics_sim::simulation::{maxwell_boltzmann_velocities, Simulation};
use tpt_dynamics_sim::system::{NewAtom, System, UnitSystem};
use tpt_dynamics_sim::thermostat::Berendsen;

const CUTOFF: f64 = 2.5;

/// Second virial coefficient of the shifted, truncated LJ(σ=ε=1) fluid at
/// reduced temperature `t`, by composite Simpson quadrature.
fn second_virial_b2(ff: &LennardJones, t: f64) -> f64 {
    let n = 200_000usize; // even
    let h = CUTOFF / n as f64;
    let integrand = |r: f64| {
        if r == 0.0 {
            return 0.0;
        }
        let (u, _) = ff.pair(0, 0, r);
        ((-u / t).exp() - 1.0) * r * r
    };
    let mut sum = integrand(0.0) + integrand(CUTOFF);
    for k in 1..n {
        let w = if k % 2 == 1 { 4.0 } else { 2.0 };
        sum += w * integrand(k as f64 * h);
    }
    let integral = sum * h / 3.0;
    -2.0 * std::f64::consts::PI * integral
}

fn fcc_system(cells: usize, density: f64, temperature: f64, seed: u64) -> System {
    // 4 atoms per conventional cell, `cells` cells per side, at reduced number density `density`.
    let n_atoms = (4 * cells * cells * cells) as f64;
    let sim_box_len = (n_atoms / density).cbrt();
    let a = sim_box_len / cells as f64;
    let mut sys = System::new(SimBox::new(sim_box_len, sim_box_len, sim_box_len).unwrap());
    let basis = [
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(0.5, 0.5, 0.0),
        Vec3::new(0.5, 0.0, 0.5),
        Vec3::new(0.0, 0.5, 0.5),
    ];
    for ix in 0..cells {
        for iy in 0..cells {
            for iz in 0..cells {
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
    maxwell_boltzmann_velocities(&mut sys, temperature, seed, UnitSystem::Reduced);
    sys
}

#[test]
fn low_density_pressure_matches_second_virial_theory() {
    let t = 1.5;
    let density = 0.02;
    let cells = 4; // 256 atoms
    let ff = LennardJones::single_type(1.0, 1.0, CUTOFF, Modifier::Shift).unwrap();
    let b2 = second_virial_b2(&ff, t);

    let sys = fcc_system(cells, density, t, 2024);
    let eval = ForceEvaluator::new(Box::new(ff), 0.3).unwrap();
    let thermostat = Berendsen::new(t, 0.1, UnitSystem::Reduced).unwrap();
    let mut sim = Simulation::new(
        sys,
        eval,
        Box::new(VelocityVerlet),
        Some(Box::new(thermostat)),
        0.003,
        UnitSystem::Reduced,
    )
    .unwrap();

    sim.run(1500, &mut []).unwrap(); // equilibrate at T
    assert!((sim.energies().temperature - t).abs() < 0.2);
    let samples: u32 = 40_000;
    let mut sum_p = 0.0;
    let mut sum_t = 0.0;
    for _ in 0..samples {
        sim.advance().unwrap();
        let e = sim.energies();
        sum_p += e.pressure;
        sum_t += e.temperature;
    }
    eprintln!("mean T = {}", sum_t / f64::from(samples));
    let measured = sum_p / f64::from(samples);
    // Compare the second-virial term P/(rho T) - 1 directly. Its size is about 3.6%
    // of P at this density, so a 10% tolerance tests B2 to a few per cent.
    let measured_term = measured / (density * t) - 1.0;
    let predicted_term = b2 * density;
    eprintln!("B2 = {b2:.5}, B2 rho = {predicted_term:.5}, measured = {measured_term:.5}");
    assert!(
        (measured_term - predicted_term).abs() < 0.10 * predicted_term.abs(),
        "measured {measured_term} vs predicted {predicted_term}"
    );
}

#[test]
fn second_virial_is_negative_at_moderate_temperature_and_small_at_high() {
    let ff = LennardJones::single_type(1.0, 1.0, CUTOFF, Modifier::Shift).unwrap();
    // Attractive well dominates at low T (B2 < 0); repulsion dominates at high T.
    assert!(second_virial_b2(&ff, 1.0) < 0.0);
    assert!(second_virial_b2(&ff, 10.0) > 0.0);
}

#[test]
fn nve_energy_drift_small_for_shifted_potential() {
    let t = 1.0;
    let ff = LennardJones::single_type(1.0, 1.0, CUTOFF, Modifier::Shift).unwrap();
    let sys = fcc_system(5, 0.7, t, 7);
    let eval = ForceEvaluator::new(Box::new(ff), 0.3).unwrap();
    let mut sim =
        Simulation::new(sys, eval, Box::new(VelocityVerlet), None, 0.001, UnitSystem::Reduced)
            .unwrap();
    let summary = sim.run(2000, &mut []).unwrap();
    assert!(
        summary.relative_energy_drift().abs() < 1e-4,
        "drift {}",
        summary.relative_energy_drift()
    );
}
