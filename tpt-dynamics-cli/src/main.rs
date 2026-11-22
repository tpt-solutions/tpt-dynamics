//! `tptdyn`: unified command-line tool for tpt-dynamics.
//!
//! Subcommands: `init` (project templates), `validate` (check a simulation config),
//! `simulate` (MD engine), `report` (HTML report of a run), `analyze` (summaries,
//! filters and histograms of CSV, Parquet or tptdyn binary tables), `histogram`
//! (1D histogram from one CSV column), `completions` and `man`. Exit codes are
//! documented in `tptdyn --help` and in the README.

#![forbid(unsafe_code)]

mod analyze;
mod cli;
mod csvin;
mod error;
mod expr;
mod histogram;
mod init;
mod logger;
mod progress;
mod report;
mod simulate;
mod table;
mod traj_analysis;
mod validate;

use std::fs;
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use clap::{CommandFactory, Parser};
use serde::Deserialize;

use crate::cli::{Cli, Command, LogLevel};
use crate::error::{input, runtime, CliError};

/// Defaults read from the `--config` file. Unknown keys are rejected.
#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
struct CliDefaults {
    log_level: Option<LogLevel>,
    quiet: Option<bool>,
    threads: Option<usize>,
    output_dir: Option<PathBuf>,
}

fn main() -> ExitCode {
    let cli = match Cli::try_parse() {
        Ok(cli) => cli,
        Err(err) => {
            let code = if err.use_stderr() { 2 } else { 0 };
            let _ = err.print();
            return ExitCode::from(code);
        }
    };
    match run(cli) {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("tptdyn: error: {err}");
            err.exit_code()
        }
    }
}

/// Resolves global options (flags override the defaults file) and dispatches.
fn run(cli: Cli) -> Result<(), CliError> {
    let defaults = match &cli.config {
        Some(path) => load_defaults(path)?,
        None => CliDefaults::default(),
    };
    let level = cli.log_level.or(defaults.log_level).unwrap_or(LogLevel::Warn);
    let quiet = cli.quiet || defaults.quiet.unwrap_or(false);
    logger::init(level.into());
    if let Some(n) = cli.threads.or(defaults.threads) {
        // The global pool can be configured once per process; a second attempt is harmless.
        if rayon::ThreadPoolBuilder::new().num_threads(n).build_global().is_err() {
            log::debug!("thread pool already initialised; --threads {n} ignored");
        }
    }

    match cli.command {
        Command::Simulate(args) => {
            let out = args.output_dir.or(defaults.output_dir);
            simulate::run(&args.sim_config, out.as_deref(), quiet, args.report)
        }
        Command::Validate(args) => validate::run(&args.sim_config),
        Command::Init(args) => init::run(&args, quiet),
        Command::Report(args) => report::run(&args, quiet),
        Command::Analyze(args) => analyze::run(&args, quiet),
        Command::Histogram(args) => histogram::run(&args, quiet),
        Command::Completions { shell } => {
            let mut cmd = Cli::command();
            clap_complete::generate(shell, &mut cmd, "tptdyn", &mut io::stdout());
            Ok(())
        }
        Command::Man { dir } => write_man_pages(dir.as_deref()),
    }
}

/// Reads the CLI defaults file.
fn load_defaults(path: &Path) -> Result<CliDefaults, CliError> {
    let text = fs::read_to_string(path)
        .map_err(|e| input(format!("cannot read --config '{}': {e}", path.display())))?;
    toml::from_str(&text).map_err(|e| input(format!("--config '{}': {e}", path.display())))
}

/// Writes `tptdyn.1` and `tptdyn-<subcommand>.1` into `dir`, or every page
/// concatenated to stdout when no directory is given.
fn write_man_pages(dir: Option<&Path>) -> Result<(), CliError> {
    let root = Cli::command();
    let mut pages = vec![("tptdyn".to_string(), root.clone())];
    for sub in root.get_subcommands() {
        let name = sub.get_name();
        if name == "help" {
            continue;
        }
        let page_name = format!("tptdyn-{name}");
        pages.push((page_name.clone(), sub.clone().bin_name(page_name)));
    }

    if let Some(dir) = dir {
        fs::create_dir_all(dir)
            .map_err(|e| runtime(format!("cannot create '{}': {e}", dir.display())))?;
        for (name, cmd) in pages {
            let path = dir.join(format!("{name}.1"));
            let mut buf = Vec::new();
            clap_mangen::Man::new(cmd)
                .render(&mut buf)
                .map_err(|e| runtime(format!("cannot render man page: {e}")))?;
            fs::write(&path, buf)
                .map_err(|e| runtime(format!("cannot write '{}': {e}", path.display())))?;
            log::info!("wrote {}", path.display());
        }
    } else {
        let mut out = io::stdout().lock();
        for (_, cmd) in pages {
            clap_mangen::Man::new(cmd)
                .render(&mut out)
                .map_err(|e| runtime(format!("cannot render man page: {e}")))?;
        }
        out.flush().map_err(|e| runtime(format!("cannot write to stdout: {e}")))?;
    }
    Ok(())
}
