//! Simulation configuration in TOML.
//!
//! Example:
//!
//! ```toml
//! seed = 42
//! units = "md"          # or "reduced" (k_B = 1, for Lennard-Jones reduced units)
//! dt = 0.004            # ps
//! nsteps = 1000
//! temperature = 0.72    # K (or reduced units)
//! integrator = "velocity-verlet"   # or "leap-frog"
//!
//! [neighbour]
//! cutoff = 2.5          # nm
//! skin = 0.3            # nm
//!
//! [lennard_jones]
//! epsilon = 1.0
//! sigma = 1.0
//! modifier = "shift"    # "none", "shift" or "switch"
//! # switch_start = 2.0  # required for "switch"
//!
//! [thermostat]
//! kind = "berendsen"
//! tau = 0.5             # ps
//!
//! [output]
//! energy_interval = 100
//! energy_file = "energy.csv"
//! trajectory_interval = 100
//! trajectory_file = "traj.bin"
//! ```
//!
//! Unknown keys are rejected so that typos fail loudly.

use std::fs::{self, File};
use std::io::BufWriter;
use std::path::PathBuf;

use serde::{Deserialize, Serialize};

use crate::barostat::{Barostat, BAR_PER_MD_PRESSURE};
use crate::constraints::{Constraints, Lincs, Shake};
use crate::error::SimError;
use crate::ffparams::{FfOptions, ParamSet};
use crate::forcefield::{ForceField, LennardJones, LjType, Modifier};
use crate::integrator::{ForceEvaluator, Integrator, LeapFrog, VelocityVerlet};
use crate::minimize::{self, MinimizeSettings};
use crate::nonbonded::Electrostatics;
use crate::output::sink::FrameSinkObserver;
use crate::output::{EnergyLogger, Observer, TrajectoryWriter};
use crate::simulation::Simulation;
use crate::system::{System, UnitSystem};
use crate::thermostat::{Berendsen, NoseHoover, Thermostat, VRescale};
use crate::water::WaterModel;

/// Mixed into the run seed so that the thermostat noise differs from the
/// velocity-initialisation stream.
const THERMOSTAT_SEED_SALT: u64 = 0x5EED_0001;
/// Mixed into the run seed for the barostat noise.
const BAROSTAT_SEED_SALT: u64 = 0x5EED_0002;

/// Integrator selection.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum IntegratorKind {
    /// Velocity Verlet (default).
    #[default]
    VelocityVerlet,
    /// Leap-frog.
    LeapFrog,
}

/// Thermostat selection.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case", deny_unknown_fields)]
pub enum ThermostatConfig {
    /// No temperature control (NVE).
    #[default]
    None,
    /// Berendsen weak coupling to `temperature`.
    Berendsen {
        /// Coupling time in ps.
        tau: f64,
    },
    /// Stochastic velocity rescaling (Bussi et al. 2007) to `temperature`.
    VRescale {
        /// Coupling time in ps.
        tau: f64,
    },
    /// Nosé–Hoover (chain length 1) to `temperature`.
    NoseHoover {
        /// Period in ps.
        tau: f64,
    },
}

/// Pressure coupling selection. Pressures are in bar and compressibility in
/// bar⁻¹ for MD units; reduced units use the engine's values directly.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case", deny_unknown_fields)]
pub enum BarostatConfig {
    /// No pressure coupling (default).
    #[default]
    None,
    /// Isotropic Berendsen coupling.
    Berendsen {
        /// Target pressure.
        pressure: f64,
        /// Isothermal compressibility.
        compressibility: f64,
        /// Coupling time in ps.
        tau: f64,
    },
    /// Isotropic C-rescale coupling.
    CRescale {
        /// Target pressure.
        pressure: f64,
        /// Isothermal compressibility.
        compressibility: f64,
        /// Coupling time in ps.
        tau: f64,
    },
    /// Isotropic (scalar) Parrinello–Rahman coupling.
    ParrinelloRahman {
        /// Target pressure.
        pressure: f64,
        /// Isothermal compressibility.
        compressibility: f64,
        /// Box oscillation period in ps.
        tau: f64,
    },
}

/// Lennard-Jones modifier selection as written in TOML.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ModifierKind {
    /// Plain truncation.
    None,
    /// Energy shifted to zero at the cutoff (default).
    #[default]
    Shift,
    /// Polynomial switch between `switch_start` and the cutoff.
    Switch,
}

/// Lennard-Jones parameters and cutoff treatment.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct LjConfig {
    /// Well depth ε (kJ/mol, or reduced).
    pub epsilon: f64,
    /// Zero-crossing distance σ (nm, or reduced).
    pub sigma: f64,
    /// Cutoff treatment.
    pub modifier: ModifierKind,
    /// Switch start in nm (only for `modifier = "switch"`).
    pub switch_start: Option<f64>,
}

impl Default for LjConfig {
    fn default() -> Self {
        Self { epsilon: 1.0, sigma: 1.0, modifier: ModifierKind::Shift, switch_start: None }
    }
}

/// Neighbour-search parameters.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct NeighbourConfig {
    /// Interaction cutoff in nm.
    pub cutoff: f64,
    /// Verlet buffer width in nm.
    pub skin: f64,
}

impl Default for NeighbourConfig {
    fn default() -> Self {
        Self { cutoff: 2.5, skin: 0.3 }
    }
}

/// Output file and interval settings. Each file is optional.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct OutputConfig {
    /// Energy log interval in steps.
    pub energy_interval: Option<u64>,
    /// Energy log file (CSV).
    pub energy_file: Option<PathBuf>,
    /// Trajectory interval in steps.
    pub trajectory_interval: Option<u64>,
    /// Trajectory file (binary, see [`crate::output`]).
    pub trajectory_file: Option<PathBuf>,
    /// On-disk trajectory format: `"plain"` (default, see [`crate::output`]) or
    /// `"chunked"` (compressed, single-precision, from `tpt-dynamics-data`).
    pub trajectory_format: TrajectoryFormat,
}

/// Trajectory file format selection.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum TrajectoryFormat {
    /// The simple `TPTTRJ01` format written by [`crate::output::TrajectoryWriter`].
    #[default]
    Plain,
    /// The chunked `TPTCTRJ1` format from `tpt-dynamics-data`, streamed during the run.
    Chunked,
}

/// Frames per chunk for the chunked trajectory format.
const CHUNK_FRAMES: u32 = 64;

/// Force-field selection.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ForceFieldKind {
    /// Lennard-Jones particles with no bonded terms and no charges (default).
    #[default]
    LennardJones,
    /// Bonded and non-bonded terms from the ffparams file named by
    /// `forcefield.params`. The system must carry the topology and atom types
    /// that [`crate::ffparams::ParamSet::assemble`] produces.
    Bonded,
}

/// Electrostatics treatment for the bonded force field.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum ElectrostaticsKind {
    /// No electrostatics (default).
    #[default]
    None,
    /// Ewald real space plus smooth PME reciprocal space.
    Pme,
}

/// `[forcefield]` table.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct ForceFieldConfig {
    /// `"lennard-jones"` (default) or `"bonded"`.
    pub kind: ForceFieldKind,
    /// Path to the ffparams TOML file. Required for `kind = "bonded"`, and
    /// resolved against the process working directory when the simulation is built.
    pub params: Option<PathBuf>,
    /// `"none"` (default) or `"pme"`. PME needs `kind = "bonded"`, since the
    /// Lennard-Jones model carries no charges.
    pub electrostatics: ElectrostaticsKind,
    /// Real-space accuracy `erfc(α r_c) = tolerance` for PME (default 1e-5).
    pub pme_tolerance: Option<f64>,
    /// B-spline order for PME, 2 to 12 (default 4).
    pub pme_order: Option<usize>,
}

/// Water models that can be used with SETTLE.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum WaterName {
    /// TIP3P.
    Tip3p,
    /// SPC/E.
    SpcE,
}

impl WaterName {
    fn model(self) -> WaterModel {
        match self {
            Self::Tip3p => WaterModel::Tip3p,
            Self::SpcE => WaterModel::SpcE,
        }
    }
}

fn default_shake_tolerance() -> f64 {
    1e-6
}
fn default_shake_max_iter() -> usize {
    1000
}
fn default_lincs_order() -> usize {
    4
}
fn default_lincs_iterations() -> usize {
    1
}
fn default_settle_tolerance() -> f64 {
    1e-8
}

/// Holonomic constraint selection. SHAKE and LINCS constrain every topology
/// bond at its length in the starting structure. SETTLE constrains the listed
/// rigid water molecules.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "kebab-case", deny_unknown_fields)]
pub enum ConstraintConfig {
    /// No constraints (default).
    #[default]
    None,
    /// SHAKE iterative solver.
    Shake {
        /// Relative length tolerance.
        #[serde(default = "default_shake_tolerance")]
        tolerance: f64,
        /// Maximum solver sweeps.
        #[serde(default = "default_shake_max_iter")]
        max_iter: usize,
    },
    /// LINCS matrix-expansion solver.
    Lincs {
        /// Expansion order (4 is typical).
        #[serde(default = "default_lincs_order")]
        order: usize,
        /// Correction passes (1 is typical).
        #[serde(default = "default_lincs_iterations")]
        iterations: usize,
    },
    /// Analytic SETTLE for rigid three-site water.
    Settle {
        /// Water model whose geometry is enforced.
        water: WaterName,
        /// Molecules as `[oxygen, hydrogen1, hydrogen2]` atom indices.
        molecules: Vec<[usize; 3]>,
        /// Position tolerance in nm.
        #[serde(default = "default_settle_tolerance")]
        tolerance: f64,
    },
}

/// Minimiser choice for `[minimize]`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum MinimizeMethod {
    /// Adaptive steepest descent (default).
    #[default]
    SteepestDescent,
    /// L-BFGS.
    Lbfgs,
}

/// `[minimize]` table: energy minimisation before the run starts.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct MinimizeConfig {
    /// Minimise before the first step when true (default false).
    pub enabled: bool,
    /// Minimisation method.
    pub method: MinimizeMethod,
    /// Maximum iterations.
    pub max_iterations: usize,
    /// Convergence threshold on the largest force component (kJ mol⁻¹ nm⁻¹ or reduced).
    pub force_tolerance: f64,
}

impl Default for MinimizeConfig {
    fn default() -> Self {
        Self {
            enabled: false,
            method: MinimizeMethod::default(),
            max_iterations: 1000,
            force_tolerance: 10.0,
        }
    }
}

/// `[checkpoint]` table.
#[derive(Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct CheckpointConfig {
    /// File that [`SimConfig::write_checkpoint`] writes. Nothing is written
    /// unless the caller invokes that method after the run.
    pub write: Option<PathBuf>,
    /// Checkpoint file to restore after the simulation is built. Its state
    /// replaces the initial velocities and positions.
    pub restart: Option<PathBuf>,
}

/// Complete run configuration.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields, default)]
pub struct SimConfig {
    /// Random seed for velocity initialisation.
    pub seed: u64,
    /// Time step in ps.
    pub dt: f64,
    /// Number of steps.
    pub nsteps: u64,
    /// Unit system: `"md"` (K, ps, kJ/mol) or `"reduced"` (`k_B` = 1).
    pub units: UnitSystem,
    /// Target/initial temperature in K (or reduced temperature).
    pub temperature: f64,
    /// Integrator.
    pub integrator: IntegratorKind,
    /// Thermostat.
    pub thermostat: ThermostatConfig,
    /// Pressure coupling.
    pub barostat: BarostatConfig,
    /// Neighbour search.
    pub neighbour: NeighbourConfig,
    /// Lennard-Jones parameters.
    pub lennard_jones: LjConfig,
    /// Output settings.
    pub output: OutputConfig,
    /// Force field: Lennard-Jones or bonded terms from an ffparams file.
    pub forcefield: ForceFieldConfig,
    /// Holonomic constraints.
    pub constraints: ConstraintConfig,
    /// Energy minimisation before the run.
    pub minimize: MinimizeConfig,
    /// Checkpoint writing and restart.
    pub checkpoint: CheckpointConfig,
    /// Reproducibility mode. When true, the neighbour-list build runs serially
    /// in this crate and the force evaluation runs inside a one-thread rayon
    /// pool, so repeated runs with the same seed are bitwise identical
    /// regardless of the number of cores. When false, the neighbour list is
    /// built in parallel (results are still collected in atom order, so they
    /// are identical in practice, but this is not guaranteed by the code).
    pub deterministic: bool,
}

impl Default for SimConfig {
    fn default() -> Self {
        Self {
            seed: 1,
            dt: 0.002,
            nsteps: 1000,
            units: UnitSystem::default(),
            temperature: 1.0,
            integrator: IntegratorKind::default(),
            thermostat: ThermostatConfig::default(),
            barostat: BarostatConfig::default(),
            neighbour: NeighbourConfig::default(),
            lennard_jones: LjConfig::default(),
            output: OutputConfig::default(),
            forcefield: ForceFieldConfig::default(),
            constraints: ConstraintConfig::default(),
            minimize: MinimizeConfig::default(),
            checkpoint: CheckpointConfig::default(),
            deterministic: false,
        }
    }
}

impl SimConfig {
    /// Parses and validates a configuration from TOML text.
    ///
    /// # Errors
    /// Returns [`SimError::Config`] for syntax errors, unknown keys or invalid values.
    pub fn from_toml_str(text: &str) -> Result<Self, SimError> {
        let cfg: Self = toml::from_str(text).map_err(|e| SimError::Config(e.to_string()))?;
        cfg.validate()?;
        Ok(cfg)
    }

    /// Serialises the configuration to TOML text.
    ///
    /// # Errors
    /// Returns [`SimError::Config`] if serialisation fails.
    pub fn to_toml_string(&self) -> Result<String, SimError> {
        toml::to_string(self).map_err(|e| SimError::Config(e.to_string()))
    }

    /// Checks all values for physical and structural validity.
    ///
    /// # Errors
    /// Returns [`SimError::Config`] describing the first invalid value.
    pub fn validate(&self) -> Result<(), SimError> {
        let bad = |msg: &str| Err(SimError::Config(msg.to_string()));
        if !(self.dt.is_finite() && self.dt > 0.0) {
            return bad("dt must be a finite number > 0");
        }
        if !(self.temperature.is_finite() && self.temperature >= 0.0) {
            return bad("temperature must be >= 0");
        }
        if !(self.neighbour.cutoff.is_finite() && self.neighbour.cutoff > 0.0) {
            return bad("neighbour.cutoff must be > 0");
        }
        if !(self.neighbour.skin.is_finite() && self.neighbour.skin >= 0.0) {
            return bad("neighbour.skin must be >= 0");
        }
        self.lennard_jones()?;
        match self.thermostat {
            ThermostatConfig::None => {}
            ThermostatConfig::Berendsen { tau } => {
                Berendsen::new(self.temperature, tau, self.units)?;
            }
            ThermostatConfig::VRescale { tau } => {
                VRescale::new(self.temperature, tau, self.units, self.seed)?;
            }
            ThermostatConfig::NoseHoover { tau } => {
                NoseHoover::new(self.temperature, tau, self.units)?;
            }
        }
        match self.barostat {
            BarostatConfig::None => {}
            BarostatConfig::Berendsen { pressure, compressibility, tau }
            | BarostatConfig::CRescale { pressure, compressibility, tau }
            | BarostatConfig::ParrinelloRahman { pressure, compressibility, tau } => {
                if !pressure.is_finite() {
                    return bad("barostat pressure must be finite");
                }
                crate::barostat::Berendsen::new(pressure, tau, compressibility)?;
            }
        }
        if self.output.energy_interval == Some(0) || self.output.trajectory_interval == Some(0) {
            return bad("output intervals must be >= 1");
        }
        if self.output.energy_interval.is_some() != self.output.energy_file.is_some() {
            return bad("output.energy_interval and output.energy_file must be given together");
        }
        if self.output.trajectory_interval.is_some() != self.output.trajectory_file.is_some() {
            return bad(
                "output.trajectory_interval and output.trajectory_file must be given together",
            );
        }
        self.validate_forcefield()?;
        self.validate_constraints()?;
        self.validate_minimize()?;
        self.validate_checkpoint()
    }

    fn validate_forcefield(&self) -> Result<(), SimError> {
        let ff = &self.forcefield;
        let cfg = |msg: &str| Err(SimError::Config(msg.to_string()));
        match ff.kind {
            ForceFieldKind::LennardJones => {
                if ff.params.is_some() {
                    return cfg("forcefield.params is only used with forcefield.kind = \"bonded\"");
                }
                if ff.electrostatics == ElectrostaticsKind::Pme {
                    return cfg(
                        "forcefield.electrostatics = \"pme\" needs forcefield.kind = \"bonded\" \
                         (the Lennard-Jones model has no charges)",
                    );
                }
            }
            ForceFieldKind::Bonded => {
                match &ff.params {
                    None => return cfg(
                        "forcefield.kind = \"bonded\" needs forcefield.params (an ffparams file)",
                    ),
                    Some(p) if p.as_os_str().is_empty() => {
                        return cfg("forcefield.params must not be an empty path")
                    }
                    Some(_) => {}
                }
            }
        }
        if ff.electrostatics == ElectrostaticsKind::None
            && (ff.pme_tolerance.is_some() || ff.pme_order.is_some())
        {
            return cfg(
                "forcefield.pme_tolerance and pme_order apply only to electrostatics = \"pme\"",
            );
        }
        if let Some(t) = ff.pme_tolerance {
            if !(t.is_finite() && t > 0.0 && t < 1.0) {
                return cfg("forcefield.pme_tolerance must be in (0, 1)");
            }
        }
        if let Some(o) = ff.pme_order {
            if !(2..=12).contains(&o) {
                return cfg("forcefield.pme_order must be between 2 and 12");
            }
        }
        Ok(())
    }

    fn validate_constraints(&self) -> Result<(), SimError> {
        let cfg = |msg: &str| Err(SimError::Config(msg.to_string()));
        match &self.constraints {
            ConstraintConfig::None => {}
            ConstraintConfig::Shake { tolerance, max_iter } => {
                if !(tolerance.is_finite() && *tolerance > 0.0) {
                    return cfg("constraints.tolerance must be > 0 for shake");
                }
                if *max_iter == 0 {
                    return cfg("constraints.max_iter must be >= 1 for shake");
                }
            }
            ConstraintConfig::Lincs { order, iterations } => {
                if *order == 0 || *iterations == 0 {
                    return cfg("constraints.order and iterations must be >= 1 for lincs");
                }
            }
            ConstraintConfig::Settle { molecules, tolerance, .. } => {
                if molecules.is_empty() {
                    return cfg("constraints.molecules must list at least one water molecule");
                }
                if !(tolerance.is_finite() && *tolerance > 0.0) {
                    return cfg("constraints.tolerance must be > 0 for settle");
                }
            }
        }
        Ok(())
    }

    fn validate_minimize(&self) -> Result<(), SimError> {
        let m = &self.minimize;
        if !m.enabled {
            return Ok(());
        }
        if m.max_iterations == 0 {
            return Err(SimError::Config("minimize.max_iterations must be >= 1".into()));
        }
        if !(m.force_tolerance.is_finite() && m.force_tolerance > 0.0) {
            return Err(SimError::Config("minimize.force_tolerance must be > 0".into()));
        }
        if self.checkpoint.restart.is_some() {
            return Err(SimError::Config(
                "minimize.enabled cannot be combined with checkpoint.restart: \
                 the restored state would be minimised away"
                    .into(),
            ));
        }
        Ok(())
    }

    fn validate_checkpoint(&self) -> Result<(), SimError> {
        for (key, path) in
            [("write", &self.checkpoint.write), ("restart", &self.checkpoint.restart)]
        {
            if path.as_ref().is_some_and(|p| p.as_os_str().is_empty()) {
                return Err(SimError::Config(format!(
                    "checkpoint.{key} must not be an empty path"
                )));
            }
        }
        Ok(())
    }

    /// The Lennard-Jones cutoff treatment described by `[lennard_jones]`.
    fn modifier(&self) -> Result<Modifier, SimError> {
        let lj = &self.lennard_jones;
        Ok(match lj.modifier {
            ModifierKind::None => Modifier::None,
            ModifierKind::Shift => Modifier::Shift,
            ModifierKind::Switch => Modifier::Switch {
                r_switch: lj.switch_start.ok_or_else(|| {
                    SimError::Config("switch_start is required for the switch modifier".into())
                })?,
            },
        })
    }

    /// Builds the Lennard-Jones force field described by this configuration.
    ///
    /// # Errors
    /// Returns [`SimError::Config`] for invalid parameters.
    pub fn lennard_jones(&self) -> Result<LennardJones, SimError> {
        let lj = &self.lennard_jones;
        LennardJones::new(
            &[LjType { epsilon: lj.epsilon, sigma: lj.sigma }],
            self.neighbour.cutoff,
            self.modifier()?,
        )
    }

    /// The force field selected by `[forcefield]`. For the bonded kind this
    /// reads the ffparams file and builds the terms for `system`.
    fn build_force_field(&self, system: &System) -> Result<Box<dyn ForceField>, SimError> {
        match self.forcefield.kind {
            ForceFieldKind::LennardJones => Ok(Box::new(self.lennard_jones()?)),
            ForceFieldKind::Bonded => {
                let path = self.forcefield.params.as_ref().ok_or_else(|| {
                    SimError::Config("forcefield.params is required for the bonded kind".into())
                })?;
                let text = fs::read_to_string(path).map_err(|e| {
                    SimError::Config(format!(
                        "cannot read forcefield.params '{}': {e}",
                        path.display()
                    ))
                })?;
                let params = ParamSet::from_toml_str(&text)?;
                let electrostatics = match self.forcefield.electrostatics {
                    ElectrostaticsKind::None => Electrostatics::None,
                    ElectrostaticsKind::Pme => Electrostatics::Pme {
                        tolerance: self.forcefield.pme_tolerance.unwrap_or(1e-5),
                        order: self.forcefield.pme_order.unwrap_or(4),
                    },
                };
                let opts = FfOptions {
                    units: self.units,
                    cutoff: self.neighbour.cutoff,
                    modifier: self.modifier()?,
                    electrostatics,
                };
                Ok(Box::new(params.force_field(system, &opts)?))
            }
        }
    }

    /// The constraint solver selected by `[constraints]`, built from the
    /// simulation's current geometry. SHAKE and LINCS take every topology bond
    /// at its present length.
    fn build_constraints(
        &self,
        sim: &Simulation,
    ) -> Result<Option<Box<dyn Constraints>>, SimError> {
        let n = sim.system.len();
        let bond_pairs = || -> Result<Vec<(usize, usize, f64)>, SimError> {
            let s = &sim.system;
            if s.topology.bonds.is_empty() {
                return Err(SimError::Config(
                    "shake and lincs constrain topology bonds, but the system has none; \
                     use kind = \"bonded\" with an assembled system, or choose another constraint"
                        .into(),
                ));
            }
            Ok(s.topology
                .bonds
                .iter()
                .map(|&(i, j)| {
                    let d = s.sim_box.min_image(s.positions[i] - s.positions[j]).norm();
                    (i, j, d)
                })
                .collect())
        };
        Ok(match &self.constraints {
            ConstraintConfig::None => None,
            ConstraintConfig::Shake { tolerance, max_iter } => {
                Some(Box::new(Shake::new(bond_pairs()?, n, *tolerance, *max_iter)?))
            }
            ConstraintConfig::Lincs { order, iterations } => {
                Some(Box::new(Lincs::new(bond_pairs()?, n, *order, *iterations)?))
            }
            ConstraintConfig::Settle { water, molecules, tolerance } => {
                Some(Box::new(water.model().settle(molecules.clone(), n, *tolerance)?))
            }
        })
    }

    /// Writes a checkpoint of `sim` to `[checkpoint] write`, if that is set.
    /// The configuration does not run the loop itself, so the caller invokes
    /// this after the run.
    ///
    /// # Errors
    /// Propagates force-evaluation errors and returns [`SimError::Io`] if the
    /// file cannot be written.
    pub fn write_checkpoint(&self, sim: &mut Simulation) -> Result<(), SimError> {
        let Some(path) = &self.checkpoint.write else {
            return Ok(());
        };
        let bytes = sim.checkpoint_bytes()?;
        fs::write(path, bytes)?;
        Ok(())
    }

    /// Creates the output observers named in the configuration, opening files
    /// for writing (truncating existing files).
    ///
    /// # Errors
    /// Returns [`SimError::Io`] if a file cannot be created.
    pub fn open_observers(&self, n_atoms: usize) -> Result<Vec<Box<dyn Observer>>, SimError> {
        let mut out: Vec<Box<dyn Observer>> = Vec::new();
        if let (Some(every), Some(path)) = (self.output.energy_interval, &self.output.energy_file) {
            let w = BufWriter::new(File::create(path)?);
            out.push(Box::new(EnergyLogger::new(w, every)?));
        }
        if let (Some(every), Some(path)) =
            (self.output.trajectory_interval, &self.output.trajectory_file)
        {
            let w = BufWriter::new(File::create(path)?);
            match self.output.trajectory_format {
                TrajectoryFormat::Plain => {
                    out.push(Box::new(TrajectoryWriter::new(w, n_atoms, every)?));
                }
                TrajectoryFormat::Chunked => {
                    let sink = tpt_dynamics_data::trajectory::TrajectoryWriter::new(
                        w,
                        n_atoms,
                        tpt_dynamics_data::trajectory::Precision::F32,
                        CHUNK_FRAMES,
                    )
                    .map_err(|e| SimError::Invalid(format!("trajectory sink: {e}")))?;
                    out.push(Box::new(FrameSinkObserver::new(sink, n_atoms, every)));
                }
            }
        }
        Ok(out)
    }

    /// Builds a ready-to-run [`Simulation`] for `system` from this configuration.
    ///
    /// Velocities are drawn from a Maxwell-Boltzmann distribution at
    /// `temperature` using `seed`, with the centre-of-mass motion removed.
    ///
    /// # Errors
    /// Returns an error if the configuration or system is invalid.
    pub fn build_simulation(&self, mut system: System) -> Result<Simulation, SimError> {
        self.validate()?;
        let ff = self.build_force_field(&system)?;
        let mut eval = ForceEvaluator::new(ff, self.neighbour.skin)?;
        eval.set_deterministic(self.deterministic)?;
        if self.minimize.enabled {
            let settings = MinimizeSettings {
                max_iterations: self.minimize.max_iterations,
                force_tolerance: self.minimize.force_tolerance,
                ..MinimizeSettings::default()
            };
            match self.minimize.method {
                MinimizeMethod::SteepestDescent => {
                    minimize::steepest_descent(&mut system, &mut eval, &settings)?;
                }
                MinimizeMethod::Lbfgs => {
                    minimize::lbfgs(&mut system, &mut eval, &settings)?;
                }
            }
        }
        let integrator: Box<dyn Integrator> = match self.integrator {
            IntegratorKind::VelocityVerlet => Box::new(VelocityVerlet),
            IntegratorKind::LeapFrog => Box::new(LeapFrog),
        };
        let thermostat: Option<Box<dyn Thermostat>> = match self.thermostat {
            ThermostatConfig::None => None,
            ThermostatConfig::Berendsen { tau } => {
                Some(Box::new(Berendsen::new(self.temperature, tau, self.units)?))
            }
            ThermostatConfig::VRescale { tau } => Some(Box::new(VRescale::new(
                self.temperature,
                tau,
                self.units,
                self.seed ^ THERMOSTAT_SEED_SALT,
            )?)),
            ThermostatConfig::NoseHoover { tau } => {
                Some(Box::new(NoseHoover::new(self.temperature, tau, self.units)?))
            }
        };
        crate::simulation::maxwell_boltzmann_velocities(
            &mut system,
            self.temperature,
            self.seed,
            self.units,
        );
        let mut sim = Simulation::new(system, eval, integrator, thermostat, self.dt, self.units)?;
        sim.barostat = self.build_barostat()?;
        if let Some(constraints) = self.build_constraints(&sim)? {
            sim = sim.with_constraints(constraints)?;
        }
        if let Some(path) = &self.checkpoint.restart {
            let bytes = fs::read(path).map_err(|e| {
                SimError::Config(format!(
                    "cannot read checkpoint.restart '{}': {e}",
                    path.display()
                ))
            })?;
            sim.restore_checkpoint(&bytes)?;
        }
        Ok(sim)
    }

    /// The barostat described by the configuration, with pressures converted to
    /// the engine's units.
    fn build_barostat(&self) -> Result<Option<Box<dyn Barostat>>, SimError> {
        // 1 kJ mol⁻¹ nm⁻³ = BAR_PER_MD_PRESSURE bar, so bar → engine divides by it.
        let (scale_p, scale_beta) = match self.units {
            UnitSystem::Md => (1.0 / BAR_PER_MD_PRESSURE, BAR_PER_MD_PRESSURE),
            UnitSystem::Reduced => (1.0, 1.0),
        };
        Ok(match self.barostat {
            BarostatConfig::None => None,
            BarostatConfig::Berendsen { pressure, compressibility, tau } => {
                Some(Box::new(crate::barostat::Berendsen::new(
                    pressure * scale_p,
                    tau,
                    compressibility * scale_beta,
                )?))
            }
            BarostatConfig::CRescale { pressure, compressibility, tau } => {
                Some(Box::new(crate::barostat::CRescale::new(
                    pressure * scale_p,
                    tau,
                    compressibility * scale_beta,
                    self.seed ^ BAROSTAT_SEED_SALT,
                )?))
            }
            BarostatConfig::ParrinelloRahman { pressure, compressibility, tau } => {
                Some(Box::new(crate::barostat::ParrinelloRahman::new(
                    pressure * scale_p,
                    tau,
                    compressibility * scale_beta,
                )?))
            }
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_full_example() {
        let text = r#"
seed = 42
dt = 0.004
nsteps = 1000
temperature = 0.72
integrator = "leap-frog"

[neighbour]
cutoff = 2.5
skin = 0.3

[lennard_jones]
epsilon = 1.0
sigma = 1.0
modifier = "switch"
switch_start = 2.0

[thermostat]
kind = "berendsen"
tau = 0.5

[output]
energy_interval = 100
energy_file = "energy.csv"
"#;
        let cfg = SimConfig::from_toml_str(text).unwrap();
        assert_eq!(cfg.seed, 42);
        assert_eq!(cfg.integrator, IntegratorKind::LeapFrog);
        assert_eq!(cfg.thermostat, ThermostatConfig::Berendsen { tau: 0.5 });
        assert_eq!(cfg.lennard_jones.switch_start, Some(2.0));
        assert_eq!(cfg.output.energy_interval, Some(100));
    }

    #[test]
    fn defaults_are_valid_and_round_trip() {
        let cfg = SimConfig::default();
        cfg.validate().unwrap();
        let text = cfg.to_toml_string().unwrap();
        let back = SimConfig::from_toml_str(&text).unwrap();
        assert_eq!(back, cfg);
    }

    #[test]
    fn rejects_unknown_keys_and_bad_values() {
        let cases = [
            "dtt = 0.1\n",
            "dt = -0.1\n",
            "dt = 0.0\n",
            "temperature = -1.0\n",
            "[neighbour]\ncutoff = 0.0\n",
            "[neighbour]\nskin = -0.1\n",
            "[lennard_jones]\nmodifier = \"switch\"\n",
            "[lennard_jones]\nmodifier = \"switch\"\nswitch_start = 3.0\n",
            "[lennard_jones]\nsigma = 0.0\n",
            "integrator = \"runge-kutta\"\n",
            "[thermostat]\nkind = \"berendsen\"\ntau = 0.0\n",
            "[thermostat]\nkind = \"nose-hoover\"\n",
            "[output]\nenergy_interval = 0\nenergy_file = \"e.csv\"\n",
            "[output]\nenergy_interval = 10\n",
            "dt = \"fast\"\n",
        ];
        for case in cases {
            assert!(SimConfig::from_toml_str(case).is_err(), "expected error for {case:?}");
        }
    }

    #[test]
    fn builds_a_runnable_simulation() {
        let cfg = SimConfig { nsteps: 5, ..SimConfig::default() };
        let mut sys = crate::integrator::tests::lj_gas(6, 0.5, 9);
        sys.remove_com_velocity();
        let mut sim = cfg.build_simulation(sys).unwrap();
        assert_eq!(sim.integrator_name(), "velocity-verlet");
        sim.advance().unwrap();
        assert_eq!(sim.step, 1);
    }
}

#[cfg(test)]
mod options_tests {
    use std::path::PathBuf;

    use tpt_dynamics_core::Vec3;

    use super::*;
    use crate::ffparams::ParamSet;
    use crate::pbc::SimBox;
    use crate::system::{element_mass, NewAtom};
    use crate::topology::Topology;

    /// Synthetic parameter file shipped with the crate (butane, reduced units).
    fn synthetic_params() -> PathBuf {
        PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/data/synthetic-test.toml"))
    }

    /// `n` rigid TIP3P-geometry waters in a 6 nm box, atoms ordered O, H, H per molecule.
    fn waters(n: usize) -> System {
        let mut s = System::new(SimBox::new(6.0, 6.0, 6.0).unwrap());
        let r = 0.09572;
        let theta = 104.52f64.to_radians();
        for m in 0..n {
            let o = Vec3::new(1.0 + 0.7 * m as f64, 1.0, 1.0);
            let hs =
                [o + Vec3::new(r, 0.0, 0.0), o + Vec3::new(r * theta.cos(), r * theta.sin(), 0.0)];
            for (name, el, p) in [("OW", "O", o), ("HW1", "H", hs[0]), ("HW2", "H", hs[1])] {
                s.push_atom(NewAtom {
                    name: name.into(),
                    residue_name: "SOL".into(),
                    residue_id: i32::try_from(m + 1).unwrap(),
                    element: el.into(),
                    position: p,
                    velocity: Vec3::ZERO,
                    mass: element_mass(el).unwrap(),
                    charge: 0.0,
                    atom_type: 0,
                });
            }
        }
        s
    }

    /// Short, gentle Lennard-Jones settings that fit a 6 nm box.
    const BASE: &str = "dt = 0.001\nnsteps = 3\n[neighbour]\ncutoff = 1.0\nskin = 0.1\n\
                        [lennard_jones]\nepsilon = 0.0\nsigma = 0.3\n";

    fn cfg(extra: &str) -> SimConfig {
        SimConfig::from_toml_str(&format!("{BASE}{extra}")).unwrap()
    }

    fn min_dist(s: &System, i: usize, j: usize) -> f64 {
        s.sim_box.min_image(s.positions[i] - s.positions[j]).norm()
    }

    #[test]
    fn full_option_set_parses() {
        let text = "seed = 7\n\
            [forcefield]\nkind = \"bonded\"\nparams = \"ff.toml\"\nelectrostatics = \"pme\"\n\
            pme_tolerance = 1e-5\npme_order = 4\n\
            [constraints]\nkind = \"lincs\"\norder = 4\niterations = 1\n\
            [minimize]\nenabled = true\nmethod = \"lbfgs\"\nmax_iterations = 200\n\
            force_tolerance = 5.0\n\
            [checkpoint]\nwrite = \"out.ckpt\"\n";
        let c = SimConfig::from_toml_str(text).unwrap();
        assert_eq!(c.forcefield.kind, ForceFieldKind::Bonded);
        assert_eq!(c.forcefield.electrostatics, ElectrostaticsKind::Pme);
        assert_eq!(c.constraints, ConstraintConfig::Lincs { order: 4, iterations: 1 });
        assert_eq!(c.minimize.method, MinimizeMethod::Lbfgs);
        assert_eq!(c.checkpoint.write, Some(PathBuf::from("out.ckpt")));
    }

    #[test]
    fn minimise_and_settle_run_and_keep_water_rigid() {
        let c = cfg("[minimize]\nenabled = true\nmethod = \"lbfgs\"\nmax_iterations = 50\n\
             [constraints]\nkind = \"settle\"\nwater = \"tip3p\"\nmolecules = [[0,1,2],[3,4,5]]\n\
             [thermostat]\nkind = \"berendsen\"\ntau = 0.5\n");
        let mut sim = c.build_simulation(waters(2)).unwrap();
        sim.run(c.nsteps, &mut []).unwrap();
        for m in 0..2 {
            for h in [3 * m + 1, 3 * m + 2] {
                let d = min_dist(&sim.system, h, 3 * m);
                assert!((d - 0.09572).abs() < 1e-6, "O-H {d}");
            }
        }
    }

    #[test]
    fn shake_keeps_topology_bonds_at_their_length() {
        let bonds = [(0, 1), (0, 2), (3, 4), (3, 5)];
        let mut sys = waters(2);
        sys.topology = Topology::from_bonds(6, &bonds).unwrap();
        let c = cfg("[constraints]\nkind = \"shake\"\n");
        let mut sim = c.build_simulation(sys).unwrap();
        let before: Vec<f64> = bonds.iter().map(|&(i, j)| min_dist(&sim.system, i, j)).collect();
        sim.run(c.nsteps, &mut []).unwrap();
        for (k, &(i, j)) in bonds.iter().enumerate() {
            let d = min_dist(&sim.system, i, j);
            assert!((d - before[k]).abs() < 1e-5, "bond {i}-{j}: {d} vs {}", before[k]);
        }
    }

    #[test]
    fn lincs_runs_on_topology_bonds() {
        let mut sys = waters(1);
        sys.topology = Topology::from_bonds(3, &[(0, 1), (0, 2)]).unwrap();
        let c = cfg("[constraints]\nkind = \"lincs\"\n");
        let mut sim = c.build_simulation(sys).unwrap();
        sim.run(c.nsteps, &mut []).unwrap();
        assert!(sim.energies().total.is_finite());
    }

    #[test]
    fn bonded_forcefield_with_pme_runs() {
        let text = std::fs::read_to_string(synthetic_params()).unwrap();
        let sys = ParamSet::from_toml_str(&text)
            .unwrap()
            .assemble(&[("butane", Vec3::new(2.0, 2.0, 2.0))], SimBox::new(8.0, 8.0, 8.0).unwrap())
            .unwrap()
            .system;
        let toml = format!(
            "units = \"reduced\"\ntemperature = 0.5\ndt = 0.001\nnsteps = 2\n\
             [neighbour]\ncutoff = 2.0\nskin = 0.2\n\
             [forcefield]\nkind = \"bonded\"\nparams = '{}'\n\
             electrostatics = \"pme\"\npme_tolerance = 1e-4\n",
            synthetic_params().display()
        );
        let c = SimConfig::from_toml_str(&toml).unwrap();
        let mut sim = c.build_simulation(sys).unwrap();
        sim.run(c.nsteps, &mut []).unwrap();
        assert!(sim.energies().total.is_finite());
    }

    #[test]
    fn checkpoint_written_then_restored_reproduces_state() {
        let path = std::env::temp_dir().join(format!("tpt-sim-config-{}.ckpt", std::process::id()));
        let mut writer = cfg("");
        writer.checkpoint.write = Some(path.clone());
        let mut sim = writer.build_simulation(waters(2)).unwrap();
        sim.run(3, &mut []).unwrap();
        writer.write_checkpoint(&mut sim).unwrap();
        let expected = sim.energies().total;

        let mut reader = cfg("");
        reader.checkpoint.restart = Some(path.clone());
        let restored = reader.build_simulation(waters(2));
        std::fs::remove_file(&path).ok();
        let restored = restored.unwrap();
        assert_eq!(restored.step, 3);
        assert!((restored.energies().total - expected).abs() <= 1e-9 * expected.abs().max(1.0));
    }

    #[test]
    fn unknown_keys_in_new_tables_are_rejected() {
        let cases = [
            "[forcefield]\nkind = \"lennard-jones\"\nbogus = 1\n",
            "[constraints]\nkind = \"shake\"\nbogus = 1\n",
            "[minimize]\nenabled = true\nsteps = 5\n",
            "[checkpoint]\nfile = \"x.ckpt\"\n",
        ];
        for case in cases {
            assert!(SimConfig::from_toml_str(&format!("{BASE}{case}")).is_err(), "{case:?}");
        }
    }

    #[test]
    fn bad_combinations_give_config_errors() {
        let cases = [
            "[forcefield]\nelectrostatics = \"pme\"\n",
            "[forcefield]\nkind = \"bonded\"\n",
            "[forcefield]\npme_order = 4\n",
            "[forcefield]\nparams = \"ff.toml\"\n",
            "[forcefield]\nkind = \"bonded\"\nparams = \"ff.toml\"\nelectrostatics = \"pme\"\npme_order = 1\n",
            "[minimize]\nenabled = true\n[checkpoint]\nrestart = \"x.ckpt\"\n",
            "[constraints]\nkind = \"settle\"\nwater = \"tip3p\"\nmolecules = []\n",
            "[constraints]\nkind = \"settle\"\nwater = \"tip4p\"\nmolecules = [[0, 1, 2]]\n",
            "[constraints]\nkind = \"shake\"\ntolerance = 0.0\n",
            "[checkpoint]\nwrite = \"\"\n",
        ];
        for case in cases {
            match SimConfig::from_toml_str(&format!("{BASE}{case}")) {
                Err(SimError::Config(_)) => {}
                other => panic!("expected Config error for {case:?}, got {other:?}"),
            }
        }
    }

    #[test]
    fn constraints_without_bonds_and_missing_params_fail_at_build() {
        let c = cfg("[constraints]\nkind = \"shake\"\n");
        assert!(matches!(c.build_simulation(waters(1)), Err(SimError::Config(_))));
        let mut b = cfg("[forcefield]\nkind = \"bonded\"\nparams = \"x\"\n");
        b.forcefield.params = Some(PathBuf::from("definitely-missing-ffparams.toml"));
        assert!(matches!(b.build_simulation(waters(1)), Err(SimError::Config(_))));
    }
}
