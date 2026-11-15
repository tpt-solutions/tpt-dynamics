//! `tptdyn analyze --trajectory`: per-frame RMSD against the first frame and
//! radius of gyration, streamed from a trajectory file into a CSV.
//!
//! Frames are read one at a time, so only the first frame and the current one are
//! held in memory. Masses come from `--structure` (a `.pdb` or `.gro` file with the
//! same atoms); without it every atom has equal weight.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Write};
use std::path::Path;

use tpt_dynamics_core::Vec3;
use tpt_dynamics_data::trajectory::{TrajectoryError, TrajectoryReader};
use tpt_dynamics_sim::analysis::{radius_of_gyration, rmsd};
use tpt_dynamics_sim::{gro, pdb};

use crate::cli::AnalyzeArgs;
use crate::error::{input, runtime, sim_err, usage, CliError};

/// Runs `analyze --trajectory`. The caller has checked that `--trajectory` is set.
pub(crate) fn run(args: &AnalyzeArgs, quiet: bool) -> Result<(), CliError> {
    let traj_path = args.trajectory.as_deref().ok_or_else(|| usage("--trajectory is missing"))?;
    let out_path = args.out.as_deref().ok_or_else(|| usage("--trajectory needs --out FILE"))?;

    let masses = args.structure.as_deref().map(read_masses).transpose()?;

    let file = File::open(traj_path)
        .map_err(|e| input(format!("cannot open trajectory '{}': {e}", traj_path.display())))?;
    let mut reader = TrajectoryReader::open(BufReader::new(file)).map_err(traj_err)?;
    let n_atoms = reader.n_atoms();
    if let Some(m) = &masses {
        if m.len() != n_atoms {
            return Err(input(format!(
                "structure has {} atoms but the trajectory has {n_atoms}",
                m.len()
            )));
        }
    }
    let weights = masses.as_deref();

    let out_file = File::create(out_path)
        .map_err(|e| runtime(format!("cannot create '{}': {e}", out_path.display())))?;
    let mut out = BufWriter::new(out_file);
    let write_err = |e: io::Error| runtime(format!("cannot write '{}': {e}", out_path.display()));
    writeln!(out, "step,time_ps,rmsd_nm,rg_nm").map_err(write_err)?;

    let mut reference: Option<Vec<Vec3>> = None;
    let mut rows = 0usize;
    for frame in reader.frames() {
        let frame = frame.map_err(traj_err)?;
        let positions: Vec<Vec3> =
            frame.positions.iter().map(|p| Vec3::new(p[0], p[1], p[2])).collect();
        let first = reference.get_or_insert_with(|| positions.clone());
        let deviation =
            rmsd(first.as_slice(), &positions, weights).map_err(|e| analysis_err(&e))?;
        let rg = radius_of_gyration(&positions, weights).map_err(|e| analysis_err(&e))?;
        writeln!(out, "{},{},{},{}", frame.step, frame.time, deviation, rg).map_err(write_err)?;
        rows += 1;
    }
    out.flush().map_err(write_err)?;

    if !quiet {
        println!(
            "analysed {rows} frame(s) of {n_atoms} atom(s) from '{}'; wrote '{}'",
            traj_path.display(),
            out_path.display()
        );
    }
    Ok(())
}

/// Reads the per-atom masses from a `.pdb` or `.gro` structure file.
fn read_masses(path: &Path) -> Result<Vec<f64>, CliError> {
    let ext = path.extension().and_then(|e| e.to_str()).map(str::to_ascii_lowercase);
    let file = File::open(path)
        .map_err(|e| input(format!("cannot open structure '{}': {e}", path.display())))?;
    let reader = BufReader::new(file);
    let system = match ext.as_deref() {
        Some("pdb") => pdb::read_pdb(reader),
        Some("gro") => gro::read_gro(reader),
        _ => {
            return Err(usage(format!(
                "--structure must be a .pdb or .gro file, got '{}'",
                path.display()
            )))
        }
    }
    .map_err(sim_err)?;
    Ok(system.masses)
}

fn traj_err(e: TrajectoryError) -> CliError {
    match e {
        TrajectoryError::Io(_) => runtime(format!("trajectory read failed: {e}")),
        other => input(format!("cannot read trajectory: {other}")),
    }
}

fn analysis_err(e: &tpt_dynamics_sim::analysis::AnalysisError) -> CliError {
    input(format!("analysis failed: {e}"))
}
