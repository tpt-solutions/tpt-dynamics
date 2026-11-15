//! SETTLE rigid-water presets: the constrained geometry must hold to tolerance
//! over a run, for each supported water model.

use tpt_dynamics_core::Vec3;
use tpt_dynamics_sim::forcefield::{ForceField, ForceOutput};
use tpt_dynamics_sim::integrator::{ForceEvaluator, VelocityVerlet};
use tpt_dynamics_sim::neighbour::VerletList;
use tpt_dynamics_sim::pbc::SimBox;
use tpt_dynamics_sim::simulation::{maxwell_boltzmann_velocities, Simulation};
use tpt_dynamics_sim::solvate::place_water_grid_geometry;
use tpt_dynamics_sim::system::UnitSystem;
use tpt_dynamics_sim::water::WaterModel;

/// No interactions: rigid molecules then move freely under the constraints only.
struct NoForce;

impl ForceField for NoForce {
    fn name(&self) -> &'static str {
        "none"
    }

    fn cutoff(&self) -> f64 {
        0.5
    }

    fn compute(
        &self,
        _positions: &[Vec3],
        _atom_types: &[u16],
        _charges: &[f64],
        _sim_box: &SimBox,
        _neighbours: &VerletList,
        _forces: &mut [Vec3],
    ) -> ForceOutput {
        ForceOutput::default()
    }
}

/// Runs `steps` steps of a free rigid-water box and returns the largest
/// deviation of any O–H or H–H distance from its target (nm).
fn run_rigid(model: WaterModel, steps: usize) -> f64 {
    let p = model.params();
    let mut sys =
        place_water_grid_geometry(SimBox::new(2.0, 2.0, 2.0).unwrap(), 1.0, p.r_oh, model.d_hh())
            .unwrap();
    let n_mol = sys.len() / 3;
    maxwell_boltzmann_velocities(&mut sys, 300.0, 4, UnitSystem::Md);
    let molecules: Vec<[usize; 3]> = (0..n_mol).map(|m| [3 * m, 3 * m + 1, 3 * m + 2]).collect();
    let settle = model.settle(molecules.clone(), sys.len(), 1e-12).unwrap();
    let eval = ForceEvaluator::new(Box::new(NoForce), 0.1).unwrap();
    let mut sim = Simulation::new(sys, eval, Box::new(VelocityVerlet), None, 0.002, UnitSystem::Md)
        .unwrap()
        .with_constraints(Box::new(settle))
        .unwrap();
    assert_eq!(sim.system.n_constraints, 3 * n_mol);
    let mut worst: f64 = 0.0;
    for _ in 0..steps {
        sim.advance().unwrap();
        let sb = sim.system.sim_box;
        let pos = &sim.system.positions;
        for m in &molecules {
            let oh1 = sb.min_image(pos[m[0]] - pos[m[1]]).norm();
            let oh2 = sb.min_image(pos[m[0]] - pos[m[2]]).norm();
            let hh = sb.min_image(pos[m[1]] - pos[m[2]]).norm();
            worst = worst.max((oh1 - p.r_oh).abs());
            worst = worst.max((oh2 - p.r_oh).abs());
            worst = worst.max((hh - model.d_hh()).abs());
        }
    }
    worst
}

#[test]
fn tip3p_geometry_is_kept_over_a_run() {
    let worst = run_rigid(WaterModel::Tip3p, 500);
    assert!(worst < 1e-8, "TIP3P worst deviation {worst} nm");
}

#[test]
fn spce_geometry_is_kept_over_a_run() {
    let worst = run_rigid(WaterModel::SpcE, 500);
    assert!(worst < 1e-8, "SPC/E worst deviation {worst} nm");
}

#[test]
fn tip4p_is_refused_as_a_rigid_model() {
    let sys = place_water_grid_geometry(SimBox::new(2.0, 2.0, 2.0).unwrap(), 1.0, 0.09572, 0.15139)
        .unwrap();
    assert!(WaterModel::Tip4p.settle(vec![[0, 1, 2]], sys.len(), 1e-12).is_err());
}
