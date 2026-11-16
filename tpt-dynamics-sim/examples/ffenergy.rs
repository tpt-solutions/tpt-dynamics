//! Prints the potential energy, positions and forces of a residue chain built from
//! an ffparams file. Used by `validation/amber_validate.py`, which compares the
//! output with OpenMM. Lennard-Jones and bonded terms only (electrostatics off).
//!
//! ```text
//! cargo run --release -p tpt-dynamics-sim --example ffenergy -- \
//!     <params.toml> <ACE,ALA,NME> <bond_length_nm> <box_nm>
//! ```
//!
//! Output lines: `potential <kJ/mol>`, `atoms <n>`, then `n` position lines and
//! `forces` followed by `n` force lines (kJ/mol/nm).

use std::process::ExitCode;

use tpt_dynamics_core::Vec3;
use tpt_dynamics_sim::ffparams::{FfOptions, ParamSet};
use tpt_dynamics_sim::forcefield::{ForceField, Modifier};
use tpt_dynamics_sim::neighbour::VerletList;
use tpt_dynamics_sim::nonbonded::Electrostatics;
use tpt_dynamics_sim::pbc::SimBox;
use tpt_dynamics_sim::system::UnitSystem;

const CUTOFF: f64 = 4.0;

fn run(args: &[String]) -> Result<(), String> {
    if args.len() != 4 {
        return Err("usage: ffenergy <params.toml> <ACE,ALA,NME> <bond_length_nm> <box_nm>".into());
    }
    let text = std::fs::read_to_string(&args[0]).map_err(|e| format!("{}: {e}", args[0]))?;
    let params = ParamSet::from_toml_str(&text).map_err(|e| e.to_string())?;
    let residues: Vec<&str> = args[1].split(',').map(str::trim).filter(|s| !s.is_empty()).collect();
    let bond_length: f64 = args[2].parse().map_err(|_| "bad bond length".to_string())?;
    let box_len: f64 = args[3].parse().map_err(|_| "bad box length".to_string())?;
    let sim_box = SimBox::new(box_len, box_len, box_len).map_err(|e| e.to_string())?;
    let origin = Vec3::new(box_len / 2.0, box_len / 2.0, box_len / 2.0);
    let asm = params
        .assemble_chains(&[(&residues, origin)], bond_length, sim_box)
        .map_err(|e| e.to_string())?;
    let system = asm.system;

    let opts = FfOptions {
        units: UnitSystem::Md,
        cutoff: CUTOFF,
        modifier: Modifier::None,
        electrostatics: Electrostatics::None,
    };
    let ff = params.force_field(&system, &opts).map_err(|e| e.to_string())?;
    let mut neighbours = VerletList::new(CUTOFF, 0.2).map_err(|e| e.to_string())?;
    neighbours.rebuild(&system.positions, &system.sim_box).map_err(|e| e.to_string())?;
    let mut forces = vec![Vec3::ZERO; system.len()];
    let out = ff.compute(
        &system.positions,
        &system.atom_types,
        &system.charges,
        &system.sim_box,
        &neighbours,
        &mut forces,
    );

    println!("potential {:?}", out.potential);
    println!("atoms {}", system.len());
    for p in &system.positions {
        println!("{:?} {:?} {:?}", p.x, p.y, p.z);
    }
    println!("forces");
    for f in &forces {
        println!("{:?} {:?} {:?}", f.x, f.y, f.z);
    }
    Ok(())
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match run(&args) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}
