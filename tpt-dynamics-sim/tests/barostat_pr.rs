//! Parrinello–Rahman (isotropic) barostat: the box responds to a pressure
//! difference and the coupled LJ run stays finite.

use tpt_dynamics_core::Vec3;
use tpt_dynamics_sim::barostat::ParrinelloRahman;
use tpt_dynamics_sim::forcefield::{LennardJones, Modifier};
use tpt_dynamics_sim::integrator::{ForceEvaluator, VelocityVerlet};
use tpt_dynamics_sim::pbc::SimBox;
use tpt_dynamics_sim::simulation::{maxwell_boltzmann_velocities, Simulation};
use tpt_dynamics_sim::system::{NewAtom, System, UnitSystem};

/// Simple cubic LJ lattice of `n³` atoms at reduced density `density`.
fn lj_cubic(n: usize, density: f64) -> System {
    let a = (1.0 / density).cbrt();
    let l = a * n as f64;
    let mut sys = System::new(SimBox::new(l, l, l).unwrap());
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                let p =
                    Vec3::new((ix as f64 + 0.5) * a, (iy as f64 + 0.5) * a, (iz as f64 + 0.5) * a);
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
    sys
}

fn build(density: f64, p0: f64) -> Simulation {
    let mut sys = lj_cubic(6, density);
    maxwell_boltzmann_velocities(&mut sys, 1.0, 9, UnitSystem::Reduced);
    let ff = LennardJones::single_type(1.0, 1.0, 2.5, Modifier::Shift).unwrap();
    let eval = ForceEvaluator::new(Box::new(ff), 0.3).unwrap();
    let pr = ParrinelloRahman::new(p0, 1.0, 0.1).unwrap();
    Simulation::new(sys, eval, Box::new(VelocityVerlet), None, 0.002, UnitSystem::Reduced)
        .unwrap()
        .with_barostat(Box::new(pr))
}

#[test]
fn box_expands_when_pressure_exceeds_target() {
    // At density 0.9 and T* = 1 the pressure is well above a target of 0.
    let mut sim = build(0.9, 0.0);
    let l0 = sim.system.sim_box.lengths.x;
    for _ in 0..100 {
        sim.advance().unwrap();
    }
    let l1 = sim.system.sim_box.lengths.x;
    assert!(l1 > l0, "box did not expand: {l0} -> {l1}");
}

#[test]
fn box_shrinks_when_pressure_is_below_target() {
    let mut sim = build(0.6, 5.0);
    let l0 = sim.system.sim_box.lengths.x;
    for _ in 0..100 {
        sim.advance().unwrap();
    }
    assert!(sim.system.sim_box.lengths.x < l0);
}

#[test]
fn short_coupled_run_stays_finite_and_bounded() {
    let mut sim = build(0.9, 0.0);
    let l0 = sim.system.sim_box.lengths.x;
    let summary = sim.run(1000, &mut []).unwrap();
    assert!(summary.last.total.is_finite() && summary.last.pressure.is_finite());
    assert!(sim.system.positions.iter().all(|p| p.x.is_finite() && p.y.is_finite()));
    let l1 = sim.system.sim_box.lengths.x;
    assert!(l1 > 0.5 * l0 && l1 < 2.0 * l0, "box length {l0} -> {l1}");
}
