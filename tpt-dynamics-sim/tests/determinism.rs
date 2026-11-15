//! Reproducibility mode: with `deterministic = true`, two runs with the same
//! seed produce bitwise identical state.

use tpt_dynamics_core::Vec3;
use tpt_dynamics_sim::config::SimConfig;
use tpt_dynamics_sim::pbc::SimBox;
use tpt_dynamics_sim::system::{NewAtom, System};

/// Cubic lattice of `n³` LJ atoms with a small deterministic jitter.
fn lj_jittered(n: usize, density: f64) -> System {
    let a = (1.0 / density).cbrt();
    let l = a * n as f64;
    let mut sys = System::new(SimBox::new(l, l, l).unwrap());
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                let k = (ix * 31 + iy * 17 + iz * 7) as f64;
                let jitter = Vec3::new((k * 0.37).sin(), (k * 0.91).cos(), (k * 0.53).sin()) * 0.03;
                let base =
                    Vec3::new((ix as f64 + 0.5) * a, (iy as f64 + 0.5) * a, (iz as f64 + 0.5) * a);
                sys.push_atom(NewAtom {
                    name: "AR".into(),
                    residue_name: "LJ".into(),
                    residue_id: 1,
                    element: "Ar".into(),
                    position: base + jitter * a,
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

fn config(deterministic: bool) -> SimConfig {
    let text = format!(
        "seed = 7\n\
         units = \"reduced\"\n\
         temperature = 1.0\n\
         dt = 0.002\n\
         nsteps = 300\n\
         deterministic = {deterministic}\n\
         \n\
         [neighbour]\n\
         cutoff = 2.5\n\
         skin = 0.3\n\
         \n\
         [thermostat]\n\
         kind = \"v-rescale\"\n\
         tau = 0.5\n"
    );
    SimConfig::from_toml_str(&text).unwrap()
}

/// Final state as raw bits, so that the comparison is exact.
fn final_state(deterministic: bool) -> (Vec<u64>, Vec<u64>, u64) {
    let cfg = config(deterministic);
    let mut sim = cfg.build_simulation(lj_jittered(6, 0.8)).unwrap();
    assert_eq!(sim.eval.is_deterministic(), deterministic);
    sim.run(cfg.nsteps, &mut []).unwrap();
    let bits = |v: &[Vec3]| -> Vec<u64> {
        v.iter().flat_map(|p| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()]).collect()
    };
    (bits(&sim.system.positions), bits(&sim.system.velocities), sim.energies().total.to_bits())
}

#[test]
fn deterministic_mode_repeats_bitwise() {
    let a = final_state(true);
    let b = final_state(true);
    assert_eq!(a.0, b.0, "positions differ between runs");
    assert_eq!(a.1, b.1, "velocities differ between runs");
    assert_eq!(a.2, b.2, "energies differ between runs");
}

#[test]
fn non_deterministic_mode_matches_in_this_test() {
    // Not a guarantee: the parallel list build is collected in atom order, so
    // the runs agree in practice here. Only deterministic mode is guaranteed.
    assert_eq!(final_state(false), final_state(false));
}
