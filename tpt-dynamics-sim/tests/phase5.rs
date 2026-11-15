//! End-to-end tests for the Phase 5 features: rigid water dynamics, barostat
//! response, triclinic NVE, and robustness of the configuration parser.

use tpt_dynamics_core::Vec3;
use tpt_dynamics_sim::barostat::Berendsen;
use tpt_dynamics_sim::config::SimConfig;
use tpt_dynamics_sim::constraints::Settle;
use tpt_dynamics_sim::forcefield::{ForceField, ForceOutput, LennardJones, Modifier};
use tpt_dynamics_sim::integrator::{ForceEvaluator, VelocityVerlet};
use tpt_dynamics_sim::neighbour::VerletList;
use tpt_dynamics_sim::pbc::SimBox;
use tpt_dynamics_sim::rng::Rng;
use tpt_dynamics_sim::simulation::{maxwell_boltzmann_velocities, Simulation};
use tpt_dynamics_sim::solvate::{place_water_grid, WATER_HH, WATER_OH};
use tpt_dynamics_sim::system::{NewAtom, System, UnitSystem};

/// A force field that exerts no force: rigid bodies then move freely.
struct ZeroForce;

impl ForceField for ZeroForce {
    fn name(&self) -> &'static str {
        "zero"
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

/// Cubic lattice of `n³` LJ atoms with small random displacements.
fn lj_lattice(n: usize, density: f64, seed: u64, tilt: [f64; 3]) -> System {
    let a = (1.0 / density).cbrt();
    let l = a * n as f64;
    let sb = SimBox::triclinic(l, l, l, tilt[0], tilt[1], tilt[2]).unwrap();
    let mut sys = System::new(sb);
    let mut rng = Rng::seed_from_u64(seed);
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                let p = Vec3::new(
                    (ix as f64 + 0.5) * a + 0.05 * a * (rng.uniform() - 0.5),
                    (iy as f64 + 0.5) * a + 0.05 * a * (rng.uniform() - 0.5),
                    (iz as f64 + 0.5) * a + 0.05 * a * (rng.uniform() - 0.5),
                );
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

fn lj_eval(cutoff: f64) -> ForceEvaluator {
    let ff = LennardJones::single_type(1.0, 1.0, cutoff, Modifier::Shift).unwrap();
    ForceEvaluator::new(Box::new(ff), 0.3).unwrap()
}

#[test]
fn rigid_water_keeps_its_geometry_under_dynamics() {
    // 2 x 2 x 2 grid of placeholder waters, 8 molecules.
    let mut sys = place_water_grid(SimBox::new(2.0, 2.0, 2.0).unwrap(), 1.0).unwrap();
    let n_mol = sys.len() / 3;
    maxwell_boltzmann_velocities(&mut sys, 1.0, 12, UnitSystem::Reduced);
    let molecules: Vec<[usize; 3]> = (0..n_mol).map(|m| [3 * m, 3 * m + 1, 3 * m + 2]).collect();
    let settle = Settle::new(molecules.clone(), sys.len(), WATER_OH, WATER_HH, 1e-12).unwrap();
    let eval = ForceEvaluator::new(Box::new(ZeroForce), 0.1).unwrap();
    let mut sim =
        Simulation::new(sys, eval, Box::new(VelocityVerlet), None, 0.002, UnitSystem::Reduced)
            .unwrap()
            .with_constraints(Box::new(settle))
            .unwrap();
    assert_eq!(sim.system.n_constraints, 3 * n_mol);
    for _ in 0..500 {
        sim.advance().unwrap();
        for m in &molecules {
            let p = &sim.system.positions;
            let oh1 = sim.system.sim_box.min_image(p[m[0]] - p[m[1]]).norm();
            let oh2 = sim.system.sim_box.min_image(p[m[0]] - p[m[2]]).norm();
            let hh = sim.system.sim_box.min_image(p[m[1]] - p[m[2]]).norm();
            assert!((oh1 - WATER_OH).abs() < 1e-9, "O-H1 {oh1}");
            assert!((oh2 - WATER_OH).abs() < 1e-9, "O-H2 {oh2}");
            assert!((hh - WATER_HH).abs() < 1e-9, "H-H {hh}");
        }
    }
    // Free rigid bodies with no forces conserve total momentum.
    assert!(sim.system.total_momentum().norm() < 1e-8);
}

#[test]
fn berendsen_barostat_expands_a_compressed_fluid() {
    // LJ at reduced density 0.9 and T = 1 has a large positive pressure, so
    // coupling to a target of 0.5 must increase the box length.
    let mut sys = lj_lattice(5, 0.9, 3, [0.0; 3]);
    let l0 = sys.sim_box.lengths.x;
    maxwell_boltzmann_velocities(&mut sys, 1.0, 4, UnitSystem::Reduced);
    let eval = lj_eval(2.0);
    let bar = Berendsen::new(0.5, 0.5, 0.05).unwrap();
    let mut sim =
        Simulation::new(sys, eval, Box::new(VelocityVerlet), None, 0.002, UnitSystem::Reduced)
            .unwrap()
            .with_barostat(Box::new(bar));
    for _ in 0..300 {
        sim.advance().unwrap();
    }
    let l1 = sim.system.sim_box.lengths.x;
    assert!(l1 > l0, "box did not expand: {l0} -> {l1}");
    assert!(sim.energies().total.is_finite());
    for p in &sim.system.positions {
        let w = sim.system.sim_box.wrap(*p);
        assert!((w - *p).norm() < 1e-9, "position outside box after scaling");
    }
}

#[test]
fn triclinic_box_conserves_energy_under_nve() {
    // A genuinely skewed cell: the neighbour list and the forces both use the
    // triclinic minimum-image convention.
    let mut sys = lj_lattice(4, 0.7, 9, [0.3, -0.2, 0.25]);
    maxwell_boltzmann_velocities(&mut sys, 0.8, 2, UnitSystem::Reduced);
    let mut sim = Simulation::new(
        sys,
        lj_eval(1.0),
        Box::new(VelocityVerlet),
        None,
        0.0005,
        UnitSystem::Reduced,
    )
    .unwrap();
    let summary = sim.run(300, &mut []).unwrap();
    assert!(
        summary.relative_energy_drift().abs() < 1e-3,
        "drift {}",
        summary.relative_energy_drift()
    );
}

#[test]
fn config_parser_survives_random_and_truncated_input() {
    let valid = r#"
seed = 42
dt = 0.004
nsteps = 100
units = "reduced"
temperature = 0.7
integrator = "leap-frog"

[neighbour]
cutoff = 2.5
skin = 0.3

[lennard_jones]
modifier = "switch"
switch_start = 2.0

[thermostat]
kind = "v-rescale"
tau = 0.5

[barostat]
kind = "c-rescale"
pressure = 1.0
compressibility = 4.5e-5
tau = 2.0
"#;
    assert!(SimConfig::from_toml_str(valid).is_ok());
    for cut in 0..valid.len() {
        let _ = SimConfig::from_toml_str(&valid[..cut]);
    }
    let alphabet = b"[]=\"\n .0123456789abcdefghijklmnopqrstuvwxyz-_#";
    let mut rng = Rng::seed_from_u64(77);
    for _ in 0..3000 {
        let len = (rng.next_u64() % 200) as usize;
        let text: String = (0..len)
            .map(|_| alphabet[(rng.next_u64() % alphabet.len() as u64) as usize] as char)
            .collect();
        let _ = SimConfig::from_toml_str(&text);
    }
}

#[test]
fn config_selects_new_thermostats_and_barostat() {
    let cfg = SimConfig::from_toml_str(
        "units = \"reduced\"\ntemperature = 1.0\nnsteps = 5\n[neighbour]\ncutoff = 1.0\nskin = 0.3\n[thermostat]\nkind = \"nose-hoover\"\ntau = 0.2\n[barostat]\nkind = \"berendsen\"\npressure = 1.0\ncompressibility = 0.1\ntau = 1.0\n",
    )
    .unwrap();
    let sys = lj_lattice(4, 0.7, 1, [0.0; 3]);
    let sim = cfg.build_simulation(sys).unwrap();
    assert_eq!(sim.thermostat.as_ref().unwrap().name(), "nose-hoover");
    assert_eq!(sim.barostat.as_ref().unwrap().name(), "berendsen");
}

#[test]
fn benchmark_config_parses_and_matches_lj_example_settings() {
    let cfg = SimConfig::from_toml_str(include_str!("../benchmarks/lj-1m.toml")).unwrap();
    assert_eq!(cfg.units, UnitSystem::Reduced);
    assert!((cfg.temperature - 0.72).abs() < 1e-12);
    assert!((cfg.neighbour.cutoff - 2.5).abs() < 1e-12);
}
