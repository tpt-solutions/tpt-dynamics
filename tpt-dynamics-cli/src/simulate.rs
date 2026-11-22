//! `tptdyn simulate`: builds the system and engine from a TOML file and runs it,
//! writing the energy log and trajectory named in `[output]`.
//!
//! The file is the engine's `SimConfig` plus one extra table, `[system]`, that
//! says where the starting structure comes from:
//!
//! ```toml
//! [system]
//! structure = "start.gro"   # .gro or .pdb file, or:
//! # cubic_cells = 3         # a cubic lattice of cubic_cells^3 identical atoms,
//! # spacing = 1.6           # nm between sites (required with cubic_cells)
//! # mass = 1.0              # amu for lattice atoms (default 1.0)
//! # params = "water.toml"    # ffparams file, with molecule, grid and spacing:
//! # molecule = "SOL"        # copies of one ffparams molecule template
//! # grid = [8, 8, 8]        # copies along x, y, z
//! # spacing = 0.31          # nm between grid sites; the box is grid * spacing
//! ```
//!
//! The molecule form needs `[forcefield] kind = "bonded"`. Its `params` file also
//! supplies the force field, so `[forcefield] params` is set from it.
//!
//! Relative paths in the file resolve against the directory that contains it.
//! With `--output-dir`, the energy and trajectory files are written into that
//! directory under their configured file names.

use std::fs::{self, File};
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::time::Instant;

use tpt_dynamics_core::Vec3;
use tpt_dynamics_sim::config::SimConfig;
use tpt_dynamics_sim::ffparams::ParamSet;
use tpt_dynamics_sim::gro::read_gro;
use tpt_dynamics_sim::output::{Observation, Observer};
use tpt_dynamics_sim::pbc::SimBox;
use tpt_dynamics_sim::pdb::read_pdb;
use tpt_dynamics_sim::simulation::{RunSummary, Simulation};
use tpt_dynamics_sim::system::{NewAtom, System};

use crate::error::{input, sim_err, usage, CliError};
use crate::progress::{self, ProgressObserver};
use crate::report;

/// A parsed configuration: where the starting system comes from, and the engine settings.
pub(crate) struct Loaded {
    /// The `[system]` table.
    pub(crate) spec: SystemSpec,
    /// Everything else, as the engine's configuration.
    pub(crate) cfg: SimConfig,
}

/// Where the starting system comes from.
#[derive(Debug, PartialEq)]
pub(crate) enum SystemSpec {
    /// A `.gro` or `.pdb` structure file.
    Structure(PathBuf),
    /// A simple cubic lattice of identical atoms.
    Cubic { cells: usize, spacing: f64, mass: f64 },
    /// Copies of one ffparams molecule template on a cubic grid; the box is `grid * spacing`.
    Molecules { params: PathBuf, molecule: String, grid: [usize; 3], spacing: f64 },
}

/// Keys allowed in the `[system]` table.
const SYSTEM_KEYS: &[&str] =
    &["structure", "cubic_cells", "spacing", "mass", "params", "molecule", "grid"];

/// Keys of `[output]` that name files.
const OUTPUT_FILE_KEYS: &[&str] = &["energy_file", "trajectory_file"];

/// Runs one simulation described by the TOML file at `config_path`.
pub(crate) fn run(
    config_path: &Path,
    output_dir: Option<&Path>,
    quiet: bool,
    write_report: bool,
) -> Result<(), CliError> {
    let Loaded { spec, cfg } = load_config(config_path, output_dir)?;
    if write_report && cfg.output.energy_file.is_none() {
        return Err(input("--report needs [output] energy_file to be set"));
    }

    let system = build_system(&spec)?;
    let mut sim = cfg.build_simulation(system).map_err(sim_err)?;
    let mut observers = cfg.open_observers(sim.system.len()).map_err(sim_err)?;
    // A live progress line replaces the plain 10% lines when stderr is a terminal.
    let live = progress::should_show(quiet);
    if live {
        observers.push(Box::new(ProgressObserver::new(cfg.nsteps)));
    }
    log::info!("simulating {} atoms for {} steps", sim.system.len(), cfg.nsteps);

    let summary = integrate(&mut sim, &mut observers, cfg.nsteps, !quiet && !live)?;
    if let Some(path) = &cfg.output.energy_file {
        log::info!("energy log written to {}", path.display());
    }
    if let Some(path) = &cfg.output.trajectory_file {
        log::info!("trajectory written to {}", path.display());
    }
    if !quiet {
        eprintln!(
            "done: {} steps in {:.2} s, {} neighbour rebuilds, relative energy drift {:.3e}",
            summary.steps,
            summary.wall_seconds,
            summary.neighbour_rebuilds,
            summary.relative_energy_drift()
        );
    }
    if let Some(energy) = &cfg.output.energy_file {
        if write_report {
            let dir =
                energy.parent().filter(|p| !p.as_os_str().is_empty()).unwrap_or(Path::new("."));
            let out = dir.join("report.html");
            report::write_report(energy, cfg.output.trajectory_file.as_deref(), &out)?;
            if !quiet {
                eprintln!("report written to {}", out.display());
            }
        }
    }
    Ok(())
}

/// Reads and parses a simulation configuration, without building anything.
///
/// Output file names are redirected into `output_dir` when given. Errors are input
/// errors whose message starts with the file name.
pub(crate) fn load_config(
    config_path: &Path,
    output_dir: Option<&Path>,
) -> Result<Loaded, CliError> {
    let shown = config_path.display();
    let text = fs::read_to_string(config_path)
        .map_err(|e| input(format!("cannot read config '{shown}': {e}")))?;
    let base = config_path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));

    let mut table: toml::Table =
        text.parse().map_err(|e| input(format!("config '{shown}': {e}")))?;
    let spec = take_system(&mut table, base)?;
    if let SystemSpec::Molecules { params, .. } = &spec {
        // The ffparams file of the molecule grid is also the force field.
        let ff = table
            .entry("forcefield")
            .or_insert_with(|| toml::Value::Table(toml::Table::new()))
            .as_table_mut()
            .ok_or_else(|| input("[forcefield] must be a table"))?;
        if ff.get("kind").and_then(|v| v.as_str()) != Some("bonded") {
            return Err(input("[system] molecule needs [forcefield] kind = \"bonded\""));
        }
        ff.insert("params".into(), toml::Value::String(params.display().to_string()));
    }
    resolve_outputs(&mut table, base, output_dir)?;
    let cfg_text = toml::to_string(&table).map_err(|e| input(format!("config '{shown}': {e}")))?;
    let cfg =
        SimConfig::from_toml_str(&cfg_text).map_err(|e| input(format!("config '{shown}': {e}")))?;
    Ok(Loaded { spec, cfg })
}

/// Removes and parses the `[system]` table.
fn take_system(table: &mut toml::Table, base: &Path) -> Result<SystemSpec, CliError> {
    let value = table
        .remove("system")
        .ok_or_else(|| input("config has no [system] table (set `structure` or `cubic_cells`)"))?;
    let sys = value.as_table().ok_or_else(|| input("[system] must be a table"))?;
    if let Some(key) = sys.keys().find(|k| !SYSTEM_KEYS.contains(&k.as_str())) {
        let expected: Vec<String> = SYSTEM_KEYS.iter().map(|k| format!("`{k}`")).collect();
        return Err(input(format!(
            "unknown key '{key}' in [system] (expected one of {})",
            expected.join(", ")
        )));
    }
    if sys.contains_key("molecule") {
        return take_molecule_grid(sys, base);
    }
    match (sys.get("structure"), sys.get("cubic_cells")) {
        (Some(path), None) => {
            if sys.contains_key("spacing") || sys.contains_key("mass") {
                return Err(input("[system] spacing and mass apply only to cubic_cells"));
            }
            let rel = path.as_str().ok_or_else(|| input("[system] structure must be a string"))?;
            Ok(SystemSpec::Structure(base.join(rel)))
        }
        (None, Some(cells)) => {
            let cells = cells
                .as_integer()
                .and_then(|n| usize::try_from(n).ok())
                .filter(|&n| n >= 1)
                .ok_or_else(|| input("[system] cubic_cells must be an integer >= 1"))?;
            let spacing = sys
                .get("spacing")
                .and_then(number)
                .ok_or_else(|| input("[system] spacing (nm) is required with cubic_cells"))?;
            let mass = sys
                .get("mass")
                .map_or(Some(1.0), number)
                .ok_or_else(|| input("[system] mass must be a number"))?;
            Ok(SystemSpec::Cubic { cells, spacing, mass })
        }
        _ => Err(input("[system] needs exactly one of `structure`, `cubic_cells` or `molecule`")),
    }
}

/// Reads an integer or float TOML value as `f64`.
fn number(v: &toml::Value) -> Option<f64> {
    v.as_float().or_else(|| v.as_integer().map(|n| n as f64))
}

/// Rewrites the output file names: under `output_dir` when given, otherwise
/// relative to the config file's directory.
fn resolve_outputs(
    table: &mut toml::Table,
    base: &Path,
    output_dir: Option<&Path>,
) -> Result<(), CliError> {
    if let Some(dir) = output_dir {
        fs::create_dir_all(dir).map_err(|e| {
            crate::error::runtime(format!("cannot create '{}': {e}", dir.display()))
        })?;
    }
    let Some(out) = table.get_mut("output").and_then(toml::Value::as_table_mut) else {
        return Ok(());
    };
    for key in OUTPUT_FILE_KEYS {
        let Some(value) = out.get_mut(*key) else { continue };
        let rel = value
            .as_str()
            .ok_or_else(|| input(format!("output.{key} must be a string")))?
            .to_string();
        let path = match output_dir {
            Some(dir) => {
                let name = Path::new(&rel)
                    .file_name()
                    .ok_or_else(|| usage(format!("output.{key} = '{rel}' has no file name")))?;
                dir.join(name)
            }
            None => base.join(&rel),
        };
        *value = toml::Value::String(path.to_string_lossy().into_owned());
    }
    Ok(())
}

/// Builds the starting system from a structure file or a cubic lattice.
pub(crate) fn build_system(spec: &SystemSpec) -> Result<System, CliError> {
    match spec {
        SystemSpec::Structure(path) => read_structure(path),
        SystemSpec::Molecules { params, molecule, grid, spacing } => {
            build_molecule_grid(params, molecule, *grid, *spacing)
        }
        SystemSpec::Cubic { cells, spacing, mass } => {
            if !(spacing.is_finite() && *spacing > 0.0) {
                return Err(input("[system] spacing must be a finite number > 0"));
            }
            if !(mass.is_finite() && *mass > 0.0) {
                return Err(input("[system] mass must be a finite number > 0"));
            }
            Ok(cubic_lattice(*cells, *spacing, *mass)?)
        }
    }
}

/// Parses the `molecule` form of `[system]`: a grid of copies of one template.
fn take_molecule_grid(sys: &toml::Table, base: &Path) -> Result<SystemSpec, CliError> {
    if ["structure", "cubic_cells", "mass"].iter().any(|k| sys.contains_key(*k)) {
        return Err(input(
            "[system] molecule cannot be combined with structure, cubic_cells or mass",
        ));
    }
    let molecule = sys
        .get("molecule")
        .and_then(|v| v.as_str())
        .filter(|s| !s.is_empty())
        .ok_or_else(|| input("[system] molecule must be a template name"))?;
    let params = sys
        .get("params")
        .and_then(|v| v.as_str())
        .ok_or_else(|| input("[system] params (an ffparams file) is required with molecule"))?;
    let counts = sys
        .get("grid")
        .and_then(|v| v.as_array())
        .filter(|a| a.len() == 3)
        .ok_or_else(|| input("[system] grid must be [nx, ny, nz]"))?;
    let mut grid = [0usize; 3];
    for (slot, v) in grid.iter_mut().zip(counts) {
        *slot = v
            .as_integer()
            .and_then(|n| usize::try_from(n).ok())
            .filter(|&n| n >= 1)
            .ok_or_else(|| input("[system] grid counts must be integers >= 1"))?;
    }
    let spacing = sys
        .get("spacing")
        .and_then(number)
        .filter(|s| s.is_finite() && *s > 0.0)
        .ok_or_else(|| input("[system] spacing (nm, > 0) is required with molecule"))?;
    Ok(SystemSpec::Molecules {
        params: base.join(params),
        molecule: molecule.to_owned(),
        grid,
        spacing,
    })
}

/// Builds the grid of molecule copies from an ffparams file with its topology.
fn build_molecule_grid(
    params: &Path,
    molecule: &str,
    grid: [usize; 3],
    spacing: f64,
) -> Result<System, CliError> {
    let text = fs::read_to_string(params)
        .map_err(|e| input(format!("cannot read [system] params '{}': {e}", params.display())))?;
    let set = ParamSet::from_toml_str(&text).map_err(sim_err)?;
    let [nx, ny, nz] = grid;
    let sim_box = SimBox::new(nx as f64 * spacing, ny as f64 * spacing, nz as f64 * spacing)
        .map_err(sim_err)?;
    let mut placements: Vec<(&str, Vec3)> = Vec::with_capacity(nx * ny * nz);
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let at = Vec3::new(i as f64, j as f64, k as f64) * spacing;
                placements.push((molecule, at));
            }
        }
    }
    Ok(set.assemble(&placements, sim_box).map_err(sim_err)?.system)
}

/// Reads a `.gro` or `.pdb` structure file, chosen by extension.
fn read_structure(path: &Path) -> Result<System, CliError> {
    let ext =
        path.extension().and_then(|e| e.to_str()).map(str::to_ascii_lowercase).unwrap_or_default();
    let file = File::open(path)
        .map_err(|e| input(format!("cannot open structure '{}': {e}", path.display())))?;
    let reader = BufReader::new(file);
    match ext.as_str() {
        "gro" => read_gro(reader).map_err(sim_err),
        "pdb" => read_pdb(reader).map_err(sim_err),
        _ => Err(input(format!(
            "structure '{}': unsupported extension (use .gro or .pdb)",
            path.display()
        ))),
    }
}

/// A cubic box of `cells^3` identical atoms on a simple cubic lattice.
fn cubic_lattice(cells: usize, spacing: f64, mass: f64) -> Result<System, CliError> {
    let edge = spacing * cells as f64;
    let sim_box = SimBox::new(edge, edge, edge).map_err(sim_err)?;
    let mut system = System::new(sim_box);
    for i in 0..cells {
        for j in 0..cells {
            for k in 0..cells {
                let p = Vec3::new(
                    (i as f64 + 0.5) * spacing,
                    (j as f64 + 0.5) * spacing,
                    (k as f64 + 0.5) * spacing,
                );
                system.push_atom(NewAtom {
                    name: "LJ".into(),
                    residue_name: "LJ".into(),
                    residue_id: 1,
                    element: "X".into(),
                    position: p,
                    velocity: Vec3::ZERO,
                    mass,
                    charge: 0.0,
                    atom_type: 0,
                });
            }
        }
    }
    Ok(system)
}

/// Runs `nsteps` steps, calling each observer at its own interval. With `show_steps`,
/// prints a line about every tenth of the run.
fn integrate(
    sim: &mut Simulation,
    observers: &mut [Box<dyn Observer>],
    nsteps: u64,
    show_steps: bool,
) -> Result<RunSummary, CliError> {
    let started = Instant::now();
    let initial = sim.energies();
    let builds_before = sim.eval.neighbours().builds();
    emit(sim, observers)?;
    let every = (nsteps / 10).max(1);
    for _ in 0..nsteps {
        sim.advance().map_err(sim_err)?;
        emit(sim, observers)?;
        if show_steps && (sim.step.is_multiple_of(every) || sim.step == nsteps) {
            let pct = 100 * sim.step / nsteps.max(1);
            eprintln!("step {}/{nsteps} ({pct}%)", sim.step);
        }
    }
    for o in observers.iter_mut() {
        o.finish().map_err(sim_err)?;
    }
    Ok(RunSummary {
        steps: nsteps,
        wall_seconds: started.elapsed().as_secs_f64(),
        neighbour_rebuilds: sim.eval.neighbours().builds() - builds_before,
        initial,
        last: sim.energies(),
    })
}

/// Sends the current state to every observer that is due at this step.
fn emit(sim: &Simulation, observers: &mut [Box<dyn Observer>]) -> Result<(), CliError> {
    let obs = Observation { system: &sim.system, energies: sim.energies() };
    for o in observers.iter_mut() {
        if sim.step.is_multiple_of(o.interval().max(1)) {
            o.observe(&obs).map_err(sim_err)?;
        }
    }
    Ok(())
}
