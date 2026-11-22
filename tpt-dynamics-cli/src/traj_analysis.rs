//! `tptdyn analyze --trajectory`: per-frame RMSD against the first frame and
//! radius of gyration, streamed from a trajectory file into a CSV.
//!
//! Both trajectory formats are accepted, detected from the file's magic bytes:
//! the chunked format of `tpt-dynamics-data` (`TPTCTRJ1`, streamed one frame at a
//! time) and the simple format written by `tptdyn simulate` (`TPTTRJ01`, which
//! the engine reads whole into memory).
//!
//! Only the first frame and the current one are held in memory for chunked input.
//! Masses come from `--structure` (a `.pdb` or `.gro` file with the same atoms);
//! without it every atom has equal weight.

use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::Path;

use tpt_dynamics_core::Vec3;
use tpt_dynamics_data::trajectory::{self, TrajectoryError, TrajectoryReader};
use tpt_dynamics_sim::analysis::{radius_of_gyration, rmsd};
use tpt_dynamics_sim::output::{read_trajectory, TRAJ_MAGIC};
use tpt_dynamics_sim::{gro, pdb};

use crate::cli::AnalyzeArgs;
use crate::error::{input, runtime, sim_err, usage, CliError};

/// One frame's worth of data, as the analysis needs it.
struct Row {
    step: u64,
    time: f64,
    positions: Vec<Vec3>,
}

/// Runs `analyze --trajectory`. The caller has checked that `--trajectory` is set.
pub(crate) fn run(args: &AnalyzeArgs, quiet: bool) -> Result<(), CliError> {
    let traj_path = args.trajectory.as_deref().ok_or_else(|| usage("--trajectory is missing"))?;
    let out_path = args.out.as_deref().ok_or_else(|| usage("--trajectory needs --out FILE"))?;

    let masses = args.structure.as_deref().map(read_masses).transpose()?;
    let weights = masses.as_deref();

    let mut file = File::open(traj_path)
        .map_err(|e| input(format!("cannot open trajectory '{}': {e}", traj_path.display())))?;
    let mut magic = [0u8; 8];
    let is_chunked = file.read_exact(&mut magic).is_ok() && &magic == trajectory::MAGIC;
    let is_plain = !is_chunked && &magic == TRAJ_MAGIC;
    if !is_chunked && !is_plain {
        return Err(input(format!(
            "'{}' is not a tptdyn trajectory (unknown magic bytes)",
            traj_path.display()
        )));
    }
    file.seek(SeekFrom::Start(0))
        .map_err(|e| runtime(format!("cannot read trajectory '{}': {e}", traj_path.display())))?;

    let out_file = File::create(out_path)
        .map_err(|e| runtime(format!("cannot create '{}': {e}", out_path.display())))?;
    let mut out = BufWriter::new(out_file);
    let write_err = |e: io::Error| runtime(format!("cannot write '{}': {e}", out_path.display()));
    writeln!(out, "step,time_ps,rmsd_nm,rg_nm").map_err(write_err)?;

    let (n_atoms, rows) = if is_chunked {
        let mut reader = TrajectoryReader::open(BufReader::new(file)).map_err(traj_err)?;
        let n_atoms = reader.n_atoms();
        check_atoms(weights, n_atoms)?;
        let rows = reader.frames().map(|f| {
            f.map_err(traj_err).map(|f| Row {
                step: f.step,
                time: f.time,
                positions: f.positions.iter().map(|p| Vec3::new(p[0], p[1], p[2])).collect(),
            })
        });
        (n_atoms, write_rows(&mut out, weights, rows, out_path)?)
    } else {
        let traj = read_trajectory(BufReader::new(file)).map_err(sim_err)?;
        let n_atoms = traj.n_atoms;
        check_atoms(weights, n_atoms)?;
        let rows = traj.frames.into_iter().map(|f| {
            Ok(Row {
                step: f.step,
                time: f.time,
                positions: f
                    .positions
                    .iter()
                    .map(|p| Vec3::new(f64::from(p[0]), f64::from(p[1]), f64::from(p[2])))
                    .collect(),
            })
        });
        (n_atoms, write_rows(&mut out, weights, rows, out_path)?)
    };
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

/// Checks that the structure file (if any) matches the trajectory's atom count.
fn check_atoms(weights: Option<&[f64]>, n_atoms: usize) -> Result<(), CliError> {
    match weights {
        Some(m) if m.len() != n_atoms => {
            Err(input(format!("structure has {} atoms but the trajectory has {n_atoms}", m.len())))
        }
        _ => Ok(()),
    }
}

/// Computes RMSD and Rg for each row and writes one CSV line per row. Returns the
/// number of rows written.
fn write_rows<W: Write, I>(
    out: &mut W,
    weights: Option<&[f64]>,
    rows: I,
    out_path: &Path,
) -> Result<usize, CliError>
where
    I: IntoIterator<Item = Result<Row, CliError>>,
{
    let write_err = |e: io::Error| runtime(format!("cannot write '{}': {e}", out_path.display()));
    let mut reference: Option<Vec<Vec3>> = None;
    let mut count = 0usize;
    for row in rows {
        let row = row?;
        let first = reference.get_or_insert_with(|| row.positions.clone());
        let deviation =
            rmsd(first.as_slice(), &row.positions, weights).map_err(|e| analysis_err(&e))?;
        let rg = radius_of_gyration(&row.positions, weights).map_err(|e| analysis_err(&e))?;
        writeln!(out, "{},{},{},{}", row.step, row.time, deviation, rg).map_err(write_err)?;
        count += 1;
    }
    Ok(count)
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
