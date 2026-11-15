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

use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

use serde::{Deserialize, Serialize};

use crate::barostat::{Barostat, BAR_PER_MD_PRESSURE};
use crate::error::SimError;
use crate::forcefield::{LennardJones, LjType, Modifier};
use crate::integrator::{ForceEvaluator, Integrator, LeapFrog, VelocityVerlet};
use crate::output::{EnergyLogger, Observer, TrajectoryWriter};
use crate::simulation::Simulation;
use crate::system::{System, UnitSystem};
use crate::thermostat::{Berendsen, NoseHoover, Thermostat, VRescale};

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
        Ok(())
    }

    /// Builds the Lennard-Jones force field described by this configuration.
    ///
    /// # Errors
    /// Returns [`SimError::Config`] for invalid parameters.
    pub fn lennard_jones(&self) -> Result<LennardJones, SimError> {
        let lj = &self.lennard_jones;
        let modifier = match lj.modifier {
            ModifierKind::None => Modifier::None,
            ModifierKind::Shift => Modifier::Shift,
            ModifierKind::Switch => Modifier::Switch {
                r_switch: lj.switch_start.ok_or_else(|| {
                    SimError::Config("switch_start is required for the switch modifier".into())
                })?,
            },
        };
        LennardJones::new(
            &[LjType { epsilon: lj.epsilon, sigma: lj.sigma }],
            self.neighbour.cutoff,
            modifier,
        )
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
            out.push(Box::new(TrajectoryWriter::new(w, n_atoms, every)?));
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
        let ff = self.lennard_jones()?;
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
        let mut eval = ForceEvaluator::new(Box::new(ff), self.neighbour.skin)?;
        eval.set_deterministic(self.deterministic)?;
        crate::simulation::maxwell_boltzmann_velocities(
            &mut system,
            self.temperature,
            self.seed,
            self.units,
        );
        let mut sim = Simulation::new(system, eval, integrator, thermostat, self.dt, self.units)?;
        sim.barostat = self.build_barostat()?;
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
