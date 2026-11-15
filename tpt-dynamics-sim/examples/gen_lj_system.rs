//! Writes a deterministic Lennard-Jones configuration as a `.gro` file.
//!
//! The atoms sit on a simple cubic lattice of `nx × ny × nz` sites. The lattice
//! constant follows from the reduced number density, and the velocities are
//! Maxwell-Boltzmann samples from a fixed seed with the centre-of-mass drift
//! removed. The same arguments always produce the same file.
//!
//! ```text
//! cargo run --release -p tpt-dynamics-sim --example gen_lj_system -- 20 20 25 lj-10k.gro
//! cargo run --release -p tpt-dynamics-sim --example gen_lj_system -- 100 100 100 lj-1m.gro
//! ```
//!
//! Optional trailing arguments: density (default 0.8442), temperature
//! (default 0.72), seed (default 12345).

use std::fs::File;
use std::io::BufWriter;
use std::process::ExitCode;

use tpt_dynamics_core::Vec3;
use tpt_dynamics_sim::gro::write_gro;
use tpt_dynamics_sim::pbc::SimBox;
use tpt_dynamics_sim::simulation::maxwell_boltzmann_velocities;
use tpt_dynamics_sim::system::{NewAtom, System, UnitSystem};

/// Builds the simple cubic Lennard-Jones configuration in reduced units.
pub fn simple_cubic_lj(dims: [usize; 3], density: f64, temperature: f64, seed: u64) -> System {
    let spacing = (1.0 / density).cbrt();
    let lengths =
        Vec3::new(dims[0] as f64 * spacing, dims[1] as f64 * spacing, dims[2] as f64 * spacing);
    let mut sys = System::new(SimBox::new(lengths.x, lengths.y, lengths.z).expect("valid box"));
    for ix in 0..dims[0] {
        for iy in 0..dims[1] {
            for iz in 0..dims[2] {
                let p = Vec3::new(ix as f64, iy as f64, iz as f64) * spacing;
                sys.push_atom(NewAtom {
                    name: "AR".into(),
                    residue_name: "AR".into(),
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
    maxwell_boltzmann_velocities(&mut sys, temperature, seed, UnitSystem::Reduced);
    sys
}

fn parse<T: std::str::FromStr>(
    args: &[String],
    i: usize,
    default: T,
    name: &str,
) -> Result<T, String> {
    match args.get(i) {
        None => Ok(default),
        Some(s) => s.parse().map_err(|_| format!("cannot parse {name} from '{s}'")),
    }
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.len() < 4 {
        eprintln!("usage: gen_lj_system NX NY NZ OUT.gro [density] [temperature] [seed]");
        return ExitCode::from(2);
    }
    let run = || -> Result<(), String> {
        let dims = [
            parse(&args, 0, 0_usize, "NX")?,
            parse(&args, 1, 0_usize, "NY")?,
            parse(&args, 2, 0_usize, "NZ")?,
        ];
        if dims.contains(&0) {
            return Err("lattice dimensions must be >= 1".into());
        }
        let density = parse(&args, 4, 0.8442_f64, "density")?;
        let temperature = parse(&args, 5, 0.72_f64, "temperature")?;
        let seed = parse(&args, 6, 12345_u64, "seed")?;
        let sys = simple_cubic_lj(dims, density, temperature, seed);
        let file = File::create(&args[3]).map_err(|e| format!("cannot create {}: {e}", args[3]))?;
        write_gro(BufWriter::new(file), &sys, "tpt-dynamics LJ lattice")
            .map_err(|e| e.to_string())?;
        println!(
            "wrote {} atoms, box {:.4} x {:.4} x {:.4} nm, to {}",
            sys.len(),
            sys.sim_box.lengths.x,
            sys.sim_box.lengths.y,
            sys.sim_box.lengths.z,
            args[3]
        );
        Ok(())
    };
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}
