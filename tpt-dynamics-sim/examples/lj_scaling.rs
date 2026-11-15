//! Timing driver for the single-domain and domain-decomposed Lennard-Jones engines.
//!
//! ```text
//! cargo run --release -p tpt-dynamics-sim --example lj_scaling -- --mode domain --side 100,100,100 --threads 20 --steps 20 --warmup 2
//! cargo run --release -p tpt-dynamics-sim --example lj_scaling -- --mode single --side 100,100,100 --threads 1 --steps 20 --profile
//! ```
//!
//! Options:
//! - `--mode single|domain` (default `single`).
//! - `--side NX,NY,NZ` simple cubic lattice sites (default `100,100,100`, 1,000,000 atoms).
//! - `--threads T`: rayon pool size for `single`, number of domains for `domain`.
//! - `--steps S` timed steps; `--warmup W` untimed steps first.
//! - `--rebalance R`: domain mode, rebalance every R steps (0 = never).
//! - `--profile`: single mode only. Times one neighbour-list rebuild and one
//!   force evaluation on the final configuration, to split the step cost.
//!
//! Output is one `key=value` line per run, for scripts to collect.

use std::time::Instant;

use tpt_dynamics_core::Vec3;
use tpt_dynamics_sim::domain::{auto_dims, run_decomposed_nve, DecompOptions};
use tpt_dynamics_sim::forcefield::{ForceField, LennardJones, Modifier};
use tpt_dynamics_sim::integrator::{ForceEvaluator, VelocityVerlet};
use tpt_dynamics_sim::neighbour::VerletList;
use tpt_dynamics_sim::pbc::SimBox;
use tpt_dynamics_sim::simulation::{maxwell_boltzmann_velocities, Simulation};
use tpt_dynamics_sim::system::{NewAtom, System, UnitSystem};

const DT: f64 = 0.004;
const CUTOFF: f64 = 2.5;
const SKIN: f64 = 0.3;
const TEMPERATURE: f64 = 0.72;

struct Args {
    mode: String,
    side: [usize; 3],
    threads: usize,
    steps: u64,
    warmup: u64,
    rebalance: usize,
    profile: bool,
}

fn parse_args() -> Result<Args, String> {
    let mut a = Args {
        mode: "single".into(),
        side: [100, 100, 100],
        threads: 1,
        steps: 20,
        warmup: 2,
        rebalance: 0,
        profile: false,
    };
    let mut it = std::env::args().skip(1);
    while let Some(flag) = it.next() {
        let mut value = || it.next().ok_or_else(|| format!("{flag} needs a value"));
        match flag.as_str() {
            "--mode" => a.mode = value()?,
            "--side" => {
                let v = value()?;
                let parts: Vec<usize> = v
                    .split(',')
                    .map(|s| s.trim().parse().map_err(|_| format!("bad --side '{v}'")))
                    .collect::<Result<_, _>>()?;
                if parts.len() != 3 || parts.contains(&0) {
                    return Err(format!("--side needs three positive integers, got '{v}'"));
                }
                a.side = [parts[0], parts[1], parts[2]];
            }
            "--threads" => a.threads = value()?.parse().map_err(|_| "bad --threads")?,
            "--steps" => a.steps = value()?.parse().map_err(|_| "bad --steps")?,
            "--warmup" => a.warmup = value()?.parse().map_err(|_| "bad --warmup")?,
            "--rebalance" => a.rebalance = value()?.parse().map_err(|_| "bad --rebalance")?,
            "--profile" => a.profile = true,
            other => return Err(format!("unknown argument '{other}'")),
        }
    }
    if a.threads == 0 {
        return Err("--threads must be >= 1".into());
    }
    if a.mode != "single" && a.mode != "domain" {
        return Err(format!("--mode must be single or domain, got '{}'", a.mode));
    }
    Ok(a)
}

/// Simple cubic Lennard-Jones lattice at reduced density 0.8442 with
/// Maxwell-Boltzmann velocities (seed 12345), as in `gen_lj_system`.
fn build(side: [usize; 3]) -> System {
    let spacing = (1.0_f64 / 0.8442).cbrt();
    let l = [0, 1, 2].map(|a| side[a] as f64 * spacing);
    let mut sys = System::new(SimBox::new(l[0], l[1], l[2]).expect("valid box"));
    for ix in 0..side[0] {
        for iy in 0..side[1] {
            for iz in 0..side[2] {
                let p = Vec3::new(ix as f64, iy as f64, iz as f64) * spacing;
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
    maxwell_boltzmann_velocities(&mut sys, TEMPERATURE, 12345, UnitSystem::Reduced);
    sys
}

/// Mean seconds for one neighbour-list rebuild and one LJ force evaluation.
fn profile(sys: &System, ff: &LennardJones) -> Result<(f64, f64), String> {
    let mut list = VerletList::new(CUTOFF, SKIN).map_err(|e| e.to_string())?;
    let reps = 5;
    let t = Instant::now();
    for _ in 0..reps {
        list.rebuild(&sys.positions, &sys.sim_box).map_err(|e| e.to_string())?;
    }
    let build = t.elapsed().as_secs_f64() / reps as f64;
    let mut forces = vec![Vec3::ZERO; sys.len()];
    let t = Instant::now();
    for _ in 0..reps {
        forces.fill(Vec3::ZERO);
        let _ = ff.compute(
            &sys.positions,
            &sys.atom_types,
            &sys.charges,
            &sys.sim_box,
            &list,
            &mut forces,
        );
    }
    let force = t.elapsed().as_secs_f64() / reps as f64;
    Ok((build, force))
}

fn run() -> Result<(), String> {
    let args = parse_args()?;
    let ff =
        LennardJones::single_type(1.0, 1.0, CUTOFF, Modifier::Shift).map_err(|e| e.to_string())?;
    let mut sys = build(args.side);
    let n = sys.len();
    let lengths = sys.sim_box.lengths;
    let threads = args.threads;

    let (wall, per_step, potential, extra) = if args.mode == "single" {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|e| e.to_string())?;
        let (wall, per_step, potential, profile_out) = pool.install(|| -> Result<_, String> {
            let eval =
                ForceEvaluator::new(Box::new(ff.clone()), SKIN).map_err(|e| e.to_string())?;
            let mut sim = Simulation::new(
                sys.clone(),
                eval,
                Box::new(VelocityVerlet),
                None,
                DT,
                UnitSystem::Reduced,
            )
            .map_err(|e| e.to_string())?;
            sim.run(args.warmup, &mut []).map_err(|e| e.to_string())?;
            let s = sim.run(args.steps, &mut []).map_err(|e| e.to_string())?;
            let profile_out = if args.profile { Some(profile(&sim.system, &ff)?) } else { None };
            Ok((
                s.wall_seconds,
                s.wall_seconds / s.steps.max(1) as f64,
                s.last.potential / n as f64,
                profile_out,
            ))
        })?;
        let extra = match profile_out {
            Some((b, f)) => {
                format!(" profile_rebuild_ms={:.2} profile_force_ms={:.2}", 1e3 * b, 1e3 * f)
            }
            None => String::new(),
        };
        (wall, per_step, potential, extra)
    } else {
        let opts =
            DecompOptions { dims: auto_dims(threads, lengths), rebalance_every: args.rebalance };
        // Untimed warm-up segment, then the timed segment on the same state.
        run_decomposed_nve(&mut sys, &ff, DT, args.warmup as usize, &opts)
            .map_err(|e| e.to_string())?;
        let s = run_decomposed_nve(&mut sys, &ff, DT, args.steps as usize, &opts)
            .map_err(|e| e.to_string())?;
        let extra = format!(
            " dims={}x{}x{} rebalances={} migrated={} ghosts_sent={} imbalance={:.3} t_migrate={:.3} t_ghost={:.3} t_force={:.3} t_kick={:.3}",
            s.dims[0], s.dims[1], s.dims[2], s.rebalances, s.migrated, s.ghosts_sent, s.imbalance,
            s.phase_seconds[0], s.phase_seconds[1], s.phase_seconds[2], s.phase_seconds[3]
        );
        (s.wall_seconds, s.wall_seconds / args.steps.max(1) as f64, s.potential / n as f64, extra)
    };

    println!(
        "mode={} atoms={} threads={} steps={} wall_s={:.4} ms_per_step={:.3} u_per_atom={:.6}{}",
        args.mode,
        n,
        threads,
        args.steps,
        wall,
        1e3 * per_step,
        potential,
        extra
    );
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e}");
        std::process::exit(2);
    }
}
