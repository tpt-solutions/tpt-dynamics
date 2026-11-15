//! Error type shared by every module of the simulation engine.

use std::fmt;

/// Errors produced by the simulation engine: I/O, malformed input files,
/// invalid configuration and invalid system state.
#[derive(Debug)]
pub enum SimError {
    /// An underlying I/O operation failed.
    Io(std::io::Error),
    /// A structure or trajectory file could not be parsed.
    Parse {
        /// Name of the file format (`"PDB"`, `"GRO"`, ...).
        format: &'static str,
        /// One-based line number where the problem was found (0 if not line based).
        line: usize,
        /// Human-readable description of the problem.
        message: String,
    },
    /// The simulation configuration is invalid.
    Config(String),
    /// The system or a requested operation is physically or structurally invalid.
    Invalid(String),
}

impl SimError {
    /// Builds a [`SimError::Parse`] for a line-based format.
    pub fn parse(format: &'static str, line: usize, message: impl Into<String>) -> Self {
        Self::Parse { format, line, message: message.into() }
    }
}

impl fmt::Display for SimError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(e) => write!(f, "I/O error: {e}"),
            Self::Parse { format, line, message } => {
                if *line == 0 {
                    write!(f, "{format} parse error: {message}")
                } else {
                    write!(f, "{format} parse error at line {line}: {message}")
                }
            }
            Self::Config(msg) => write!(f, "invalid configuration: {msg}"),
            Self::Invalid(msg) => write!(f, "invalid system: {msg}"),
        }
    }
}

impl std::error::Error for SimError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(e) => Some(e),
            _ => None,
        }
    }
}

impl From<std::io::Error> for SimError {
    fn from(e: std::io::Error) -> Self {
        Self::Io(e)
    }
}
