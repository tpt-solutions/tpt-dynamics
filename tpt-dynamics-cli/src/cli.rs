//! Command-line grammar for `tptdyn`, defined with the clap derive API.

use std::path::PathBuf;

use clap::{Parser, Subcommand, ValueEnum};
use clap_complete::Shell;
use serde::Deserialize;

/// Documented exit codes, shown at the end of `tptdyn --help`.
pub(crate) const EXIT_CODES: &str = "\
EXIT CODES:
  0  success
  1  runtime error (simulation failed, output could not be written)
  2  usage error (unknown option, missing or malformed argument)
  3  input or configuration error (missing or malformed file, invalid config)";

/// Unified command-line tool for tpt-dynamics.
#[derive(Debug, Parser)]
#[command(
    name = "tptdyn",
    version,
    about = "tpt-dynamics command-line tool",
    propagate_version = true,
    after_help = EXIT_CODES
)]
pub(crate) struct Cli {
    /// Diagnostic verbosity on stderr [default: warn].
    #[arg(long, global = true, value_name = "LEVEL", value_enum)]
    pub(crate) log_level: Option<LogLevel>,

    /// Suppress progress output and informational messages.
    #[arg(short, long, global = true)]
    pub(crate) quiet: bool,

    /// CLI defaults file (TOML) with optional keys `log_level`, `quiet`,
    /// `threads` and `output_dir`. Command-line flags override it.
    #[arg(long, global = true, value_name = "FILE")]
    pub(crate) config: Option<PathBuf>,

    /// Worker threads for parallel work [default: all logical cores].
    #[arg(long, global = true, value_name = "N", value_parser = parse_positive)]
    pub(crate) threads: Option<usize>,

    /// The subcommand to run.
    #[command(subcommand)]
    pub(crate) command: Command,
}

/// The `tptdyn` subcommands.
#[derive(Debug, Subcommand)]
pub(crate) enum Command {
    /// Run a molecular dynamics simulation from a TOML configuration.
    Simulate(SimulateArgs),
    /// Check a simulation configuration without running it.
    Validate(ValidateArgs),
    /// Scaffold a runnable project from a built-in template.
    Init(InitArgs),
    /// Summarise, filter and histogram an event table (CSV, Parquet or tptdyn binary).
    Analyze(AnalyzeArgs),
    /// Fill a 1D histogram from one CSV column.
    Histogram(HistogramArgs),
    /// Write a self-contained HTML report (energy, temperature, pressure) for a run.
    Report(ReportArgs),
    /// Print a shell completion script to stdout.
    Completions {
        /// Target shell.
        #[arg(value_name = "SHELL")]
        shell: Shell,
    },
    /// Write man pages for `tptdyn` and its subcommands.
    Man {
        /// Directory to write `tptdyn*.1` files into. Without it, all pages go to stdout.
        #[arg(long, value_name = "DIR")]
        dir: Option<PathBuf>,
    },
}

/// Arguments of `tptdyn simulate`.
#[derive(Debug, clap::Args)]
pub(crate) struct SimulateArgs {
    /// Simulation configuration (TOML). Relative paths inside it resolve against
    /// the directory that contains this file.
    #[arg(value_name = "CONFIG.toml")]
    pub(crate) sim_config: PathBuf,

    /// Write the energy log and trajectory here, keeping their configured file
    /// names, instead of the paths in `[output]`.
    #[arg(long, value_name = "DIR")]
    pub(crate) output_dir: Option<PathBuf>,

    /// After the run, write report.html next to the energy log (see `tptdyn report`).
    #[arg(long)]
    pub(crate) report: bool,
}

/// Arguments of `tptdyn validate`.
#[derive(Debug, clap::Args)]
pub(crate) struct ValidateArgs {
    /// Simulation configuration (TOML) to check. Nothing is run or written.
    #[arg(value_name = "CONFIG.toml")]
    pub(crate) sim_config: PathBuf,
}

/// Arguments of `tptdyn init`.
#[derive(Debug, clap::Args)]
pub(crate) struct InitArgs {
    /// Template to scaffold.
    #[arg(value_name = "TEMPLATE", value_enum)]
    pub(crate) template: Template,

    /// Directory to write into (created if missing). Default: the current directory.
    #[arg(value_name = "DIR", default_value = ".")]
    pub(crate) dir: PathBuf,

    /// Overwrite files that already exist.
    #[arg(long)]
    pub(crate) force: bool,
}

/// Project templates for `tptdyn init`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum Template {
    /// Lennard-Jones fluid on a cubic lattice: simulate config with energy and trajectory output.
    LjFluid,
    /// Synthetic event CSV generator plus an analyze/histogram command sheet.
    HepAnalysis,
}

/// Arguments of `tptdyn report`.
#[derive(Debug, clap::Args)]
pub(crate) struct ReportArgs {
    /// Run directory containing energy.csv. The report is written here unless --out is given.
    #[arg(value_name = "RUN_DIR")]
    pub(crate) run_dir: PathBuf,

    /// Trajectory file to summarise (frame count and atoms), relative to `RUN_DIR` if not absolute.
    #[arg(long, value_name = "FILE")]
    pub(crate) trajectory: Option<PathBuf>,

    /// Output HTML file. Default: `RUN_DIR/report.html`.
    #[arg(long, value_name = "FILE")]
    pub(crate) out: Option<PathBuf>,
}

/// Arguments of `tptdyn analyze`.
#[derive(Debug, clap::Args)]
pub(crate) struct AnalyzeArgs {
    /// Event table: CSV with a header row (use --no-header otherwise), a Parquet file,
    /// or a tptdyn columnar/row-store file. The format is chosen by magic bytes; anything
    /// else is read as CSV. Not used with --trajectory.
    #[arg(
        value_name = "INPUT.csv",
        required_unless_present = "trajectory",
        conflicts_with = "trajectory"
    )]
    pub(crate) input: Option<PathBuf>,

    /// Analyse a trajectory file instead of a CSV: writes per-frame RMSD against the
    /// first frame and radius of gyration to --out. Example: --trajectory run.trj
    #[arg(long, value_name = "FILE", conflicts_with_all = ["filters", "wheres", "hists"])]
    pub(crate) trajectory: Option<PathBuf>,

    /// Structure file (.pdb or .gro) for the atom masses used by --trajectory.
    /// Without it, every atom has equal weight.
    #[arg(long, value_name = "FILE", requires = "trajectory")]
    pub(crate) structure: Option<PathBuf>,

    /// CSV file for --trajectory. Columns: `step`, `time_ps`, `rmsd_nm`, `rg_nm`.
    #[arg(long, value_name = "FILE", requires = "trajectory")]
    pub(crate) out: Option<PathBuf>,

    /// Field separator for CSV input (a single ASCII character).
    #[arg(long, value_name = "CHAR", default_value_t = ',')]
    pub(crate) delimiter: char,

    /// The first line is data, not a header; columns are named c0, c1, ... (CSV only).
    #[arg(long)]
    pub(crate) no_header: bool,

    /// Keep rows where `COLUMN OP VALUE` holds. OP is one of > >= < <= == !=.
    /// Repeat the option to combine filters with AND. Example: --filter 'pt > 20'
    /// Filters and --where may be combined; every condition must hold.
    #[arg(long = "filter", value_name = "'COLUMN OP VALUE'")]
    pub(crate) filters: Vec<String>,

    /// Keep rows where the EXPR holds. EXPR may use column names, numbers, + - * /,
    /// comparisons (< <= > >= == !=), and/or/not, and parentheses.
    /// Repeat the option to combine with AND. Example: --where '(pt > 20 or `n_jets` >= 2) and eta < 2.5'
    #[arg(long = "where", value_name = "'EXPR'")]
    pub(crate) wheres: Vec<String>,

    /// Write a uniform histogram of COLUMN: COLUMN,BINS,LO,HI. Repeatable.
    /// Example: --hist pt,50,0,100
    #[arg(long = "hist", value_name = "COLUMN,BINS,LO,HI")]
    pub(crate) hists: Vec<String>,

    /// Directory for histogram files, named COLUMN.json or COLUMN.bin.
    #[arg(long, value_name = "DIR", default_value = ".")]
    pub(crate) out_dir: PathBuf,

    /// Histogram file format.
    #[arg(long, value_name = "FORMAT", value_enum, default_value_t = Format::Json)]
    pub(crate) format: Format,
}

/// Arguments of `tptdyn histogram`.
#[derive(Debug, clap::Args)]
pub(crate) struct HistogramArgs {
    /// CSV file to read.
    #[arg(value_name = "INPUT.csv")]
    pub(crate) input: PathBuf,

    /// Column name (with a header) or zero-based index (with --no-header).
    #[arg(long, value_name = "NAME|INDEX")]
    pub(crate) column: String,

    /// Number of uniform bins.
    #[arg(long, value_name = "N", value_parser = parse_positive)]
    pub(crate) bins: usize,

    /// Histogram range as LO,HI.
    #[arg(long, value_name = "LO,HI", allow_hyphen_values = true)]
    pub(crate) range: String,

    /// Output file.
    #[arg(long, value_name = "FILE")]
    pub(crate) out: PathBuf,

    /// Output format. Default: `bin` for a .bin file name, otherwise `json`.
    #[arg(long, value_name = "FORMAT", value_enum)]
    pub(crate) format: Option<Format>,

    /// Field separator (a single ASCII character).
    #[arg(long, value_name = "CHAR", default_value_t = ',')]
    pub(crate) delimiter: char,

    /// The first line is data, not a header.
    #[arg(long)]
    pub(crate) no_header: bool,
}

/// Diagnostic verbosity levels.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, ValueEnum, Deserialize)]
#[serde(rename_all = "lowercase")]
pub(crate) enum LogLevel {
    /// Log nothing.
    Off,
    /// Errors only.
    Error,
    /// Warnings and errors (default).
    Warn,
    /// Informational messages such as written output files.
    Info,
    /// Debugging detail.
    Debug,
    /// Everything.
    Trace,
}

impl From<LogLevel> for log::LevelFilter {
    fn from(level: LogLevel) -> Self {
        match level {
            LogLevel::Off => Self::Off,
            LogLevel::Error => Self::Error,
            LogLevel::Warn => Self::Warn,
            LogLevel::Info => Self::Info,
            LogLevel::Debug => Self::Debug,
            LogLevel::Trace => Self::Trace,
        }
    }
}

/// Histogram file formats.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub(crate) enum Format {
    /// Self-describing JSON text.
    Json,
    /// Checksummed binary (tpt-dynamics-data serialise format).
    Bin,
}

/// Parses an integer of at least 1 (used for `--threads` and `--bins`).
fn parse_positive(text: &str) -> Result<usize, String> {
    match text.parse::<usize>() {
        Ok(n) if n >= 1 => Ok(n),
        _ => Err(format!("'{text}' must be an integer >= 1")),
    }
}
