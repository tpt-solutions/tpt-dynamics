//! Error classes and their documented exit codes.

use std::fmt;
use std::process::ExitCode;

use tpt_dynamics_sim::SimError;

/// A failure reported by `tptdyn`. The variant selects the exit code.
#[derive(Debug)]
pub(crate) enum CliError {
    /// Invalid arguments or option values (exit 2).
    Usage(String),
    /// Missing, unreadable or malformed input or configuration (exit 3).
    Input(String),
    /// A failure while running or writing results (exit 1).
    Runtime(String),
}

impl CliError {
    /// The process exit code for this error class.
    pub(crate) fn exit_code(&self) -> ExitCode {
        ExitCode::from(match self {
            Self::Runtime(_) => 1,
            Self::Usage(_) => 2,
            Self::Input(_) => 3,
        })
    }
}

impl fmt::Display for CliError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Usage(m) | Self::Input(m) | Self::Runtime(m) => f.write_str(m),
        }
    }
}

/// Maps an engine error: I/O failures are runtime errors, everything else
/// (bad configuration, malformed files, invalid systems) is an input error.
pub(crate) fn sim_err(e: SimError) -> CliError {
    match e {
        SimError::Io(_) => CliError::Runtime(e.to_string()),
        other => CliError::Input(other.to_string()),
    }
}

/// Shorthand for an input error.
pub(crate) fn input(msg: impl Into<String>) -> CliError {
    CliError::Input(msg.into())
}

/// Shorthand for a usage error.
pub(crate) fn usage(msg: impl Into<String>) -> CliError {
    CliError::Usage(msg.into())
}

/// Shorthand for a runtime error.
pub(crate) fn runtime(msg: impl Into<String>) -> CliError {
    CliError::Runtime(msg.into())
}
