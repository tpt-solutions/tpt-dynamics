//! `tptdyn validate CONFIG.toml`: parses a configuration the way `simulate` does,
//! builds the starting system and engine (without stepping), and prints a summary.
//!
//! Unknown keys get a "did you mean" hint: the closest expected key within an edit
//! distance of 2, taken from the list the parser reports.

use std::path::Path;

use crate::error::{input, sim_err, CliError};
use crate::simulate::{build_system, load_config};

/// Largest edit distance at which a key is suggested.
const MAX_SUGGESTION_DISTANCE: usize = 2;

/// Bytes per atom in one plain trajectory frame (single-precision position).
const BYTES_PER_ATOM: u64 = 12;
/// Bytes per plain trajectory frame before the positions: step, time and box.
const FRAME_HEADER_BYTES: u64 = 8 + 8 + 3 * 8;
/// Bytes of the plain trajectory file header: magic, version and atom count.
const FILE_HEADER_BYTES: u64 = 8 + 4 + 8;

/// Runs the `validate` subcommand.
pub(crate) fn run(config_path: &Path) -> Result<(), CliError> {
    let shown = config_path.display().to_string();
    let loaded = load_config(config_path, None).map_err(with_hint)?;
    let system = build_system(&loaded.spec)?;
    let atoms = system.len();
    let l = system.sim_box.lengths;
    // Building the simulation checks the system against the force field settings.
    loaded.cfg.build_simulation(system).map_err(sim_err)?;

    let cfg = &loaded.cfg;
    println!("config ok: {shown}");
    println!("  atoms         {atoms}");
    println!("  box (nm)      {:.3} x {:.3} x {:.3}", l.x, l.y, l.z);
    println!(
        "  steps         {} (dt {} ps, {:.3} ps in total)",
        cfg.nsteps,
        cfg.dt,
        cfg.dt * cfg.nsteps as f64
    );
    if let (Some(every), Some(file)) = (cfg.output.energy_interval, &cfg.output.energy_file) {
        let rows = cfg.nsteps / every.max(1) + 1;
        println!("  energy log    {rows} rows -> {}", file.display());
    }
    if let (Some(every), Some(file)) = (cfg.output.trajectory_interval, &cfg.output.trajectory_file)
    {
        let frames = cfg.nsteps / every.max(1) + 1;
        let mb = trajectory_bytes(frames, atoms as u64) as f64 / 1.0e6;
        let note = match cfg.output.trajectory_format {
            tpt_dynamics_sim::config::TrajectoryFormat::Plain => "",
            tpt_dynamics_sim::config::TrajectoryFormat::Chunked => {
                ", upper bound: chunked output is compressed"
            }
        };
        println!("  trajectory    {frames} frames, about {mb:.3} MB{note} -> {}", file.display());
    }
    Ok(())
}

/// Size in bytes of a plain trajectory with `frames` frames of `atoms` atoms.
fn trajectory_bytes(frames: u64, atoms: u64) -> u64 {
    FILE_HEADER_BYTES + frames * (FRAME_HEADER_BYTES + BYTES_PER_ATOM * atoms)
}

/// Adds a "did you mean" line to an input error that reports an unknown key.
fn with_hint(err: CliError) -> CliError {
    match err {
        CliError::Input(msg) => match unknown_key_hint(&msg) {
            Some(hint) => input(format!("{msg}\n  hint: {hint}")),
            None => input(msg),
        },
        other => other,
    }
}

/// Finds the unknown key named in a parser message and suggests the closest
/// expected key. The expected keys are the backquoted (or single-quoted) names that
/// follow the word "expected" in the same message.
fn unknown_key_hint(message: &str) -> Option<String> {
    let (key, from) = if let Some(i) = message.find("unknown field `") {
        let start = i + "unknown field `".len();
        let end = message[start..].find('`')? + start;
        (&message[start..end], end)
    } else {
        let i = message.find("unknown key '")?;
        let start = i + "unknown key '".len();
        let end = message[start..].find('\'')? + start;
        (&message[start..end], end)
    };
    let expected_at = message[from..].find("expected")? + from;
    let candidates = quoted_names(&message[expected_at..]);
    let best = candidates
        .iter()
        .map(|c| (levenshtein(key, c), c))
        .filter(|(d, _)| *d <= MAX_SUGGESTION_DISTANCE)
        .min_by_key(|(d, _)| *d)?;
    Some(format!("did you mean `{}`?", best.1))
}

/// Names in backquotes or single quotes, in order of appearance.
fn quoted_names(text: &str) -> Vec<&str> {
    let mut names = Vec::new();
    let mut rest = text;
    while let Some(open) = rest.find(['`', '\'']) {
        let quote = rest.as_bytes()[open] as char;
        let after = &rest[open + 1..];
        let Some(close) = after.find(quote) else { break };
        let name = &after[..close];
        if !name.is_empty() && !name.contains(' ') {
            names.push(name);
        }
        rest = &after[close + 1..];
    }
    names
}

/// Classic edit distance (insertions, deletions and substitutions).
pub(crate) fn levenshtein(a: &str, b: &str) -> usize {
    let b: Vec<char> = b.chars().collect();
    let mut prev: Vec<usize> = (0..=b.len()).collect();
    for (i, ca) in a.chars().enumerate() {
        let mut cur = vec![i + 1; b.len() + 1];
        for (j, cb) in b.iter().enumerate() {
            let subst = prev[j] + usize::from(ca != *cb);
            cur[j + 1] = subst.min(prev[j + 1] + 1).min(cur[j] + 1);
        }
        prev = cur;
    }
    prev[b.len()]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn levenshtein_basics() {
        assert_eq!(levenshtein("", ""), 0);
        assert_eq!(levenshtein("abc", ""), 3);
        assert_eq!(levenshtein("energ_interval", "energy_interval"), 1);
        assert_eq!(levenshtein("kitten", "sitting"), 3);
    }

    #[test]
    fn suggests_closest_expected_key() {
        let msg =
            "unknown field `energ_interval`, expected one of `energy_interval`, `energy_file`";
        assert_eq!(unknown_key_hint(msg).as_deref(), Some("did you mean `energy_interval`?"));
    }

    #[test]
    fn suggests_for_single_expected_key_and_system_table() {
        let msg = "unknown field `cutof`, expected `cutoff` or `skin`";
        assert_eq!(unknown_key_hint(msg).as_deref(), Some("did you mean `cutoff`?"));
        let sys = "unknown key 'cubic_cell' in [system] (expected one of `structure`, `cubic_cells`, `spacing`, `mass`)";
        assert_eq!(unknown_key_hint(sys).as_deref(), Some("did you mean `cubic_cells`?"));
    }

    #[test]
    fn no_suggestion_when_far_from_every_key() {
        let msg = "unknown field `zzzzzz`, expected `cutoff` or `skin`";
        assert_eq!(unknown_key_hint(msg), None);
        assert_eq!(unknown_key_hint("invalid value for dt"), None);
    }

    #[test]
    fn trajectory_size_for_plain_format() {
        // 3 frames of 2 atoms: 20-byte header, then 3 * (40-byte frame header + 2 * 12).
        assert_eq!(trajectory_bytes(3, 2), 20 + 3 * 64);
    }
}
