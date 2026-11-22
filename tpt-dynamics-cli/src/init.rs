//! `tptdyn init TEMPLATE [DIR]`: scaffolds a runnable project from templates that
//! are embedded in the binary.
//!
//! Every target file is checked before anything is written, so a refused init
//! leaves the directory untouched. `--force` replaces existing files.

use std::fs;
use std::path::Path;

use crate::cli::{InitArgs, Template};
use crate::error::{input, runtime, CliError};

/// One file of a template: its path inside the project and its text.
type TemplateFile = (&'static str, String);

/// Molecules of the water-box template: the grid's 8 x 8 x 8 copies of SOL.
const WATER_BOX_MOLECULES: usize = 8 * 8 * 8;

/// The files of a template.
fn files(template: Template) -> Vec<TemplateFile> {
    match template {
        Template::LjFluid => vec![
            ("simulate.toml", include_str!("../templates/lj-fluid/simulate.toml").to_owned()),
            ("README.txt", include_str!("../templates/lj-fluid/README.txt").to_owned()),
        ],
        Template::HepAnalysis => vec![
            ("make_events.py", include_str!("../templates/hep-analysis/make_events.py").to_owned()),
            ("analyze.txt", include_str!("../templates/hep-analysis/analyze.txt").to_owned()),
            ("README.txt", include_str!("../templates/hep-analysis/README.txt").to_owned()),
        ],
        Template::WaterBox => vec![
            ("simulate.toml", water_box_config()),
            ("README.txt", include_str!("../templates/water-box/README.txt").to_owned()),
        ],
    }
}

/// The water-box `simulate.toml`, with the SETTLE molecule list written out: one
/// `[O, H1, H2]` triple of atom indices per molecule, eight triples per line.
fn water_box_config() -> String {
    let triples: Vec<String> = (0..WATER_BOX_MOLECULES)
        .map(|m| format!("[{}, {}, {}]", 3 * m, 3 * m + 1, 3 * m + 2))
        .collect();
    let lines: Vec<String> = triples.chunks(8).map(|c| format!("  {},", c.join(", "))).collect();
    let list = format!("\n{}\n", lines.join("\n"));
    include_str!("../templates/water-box/simulate.toml").replace("@MOLECULES@", &list)
}

/// Runs the `init` subcommand.
pub(crate) fn run(args: &InitArgs, quiet: bool) -> Result<(), CliError> {
    write_template(args.template, &args.dir, args.force)?;
    if !quiet {
        let targets = files(args.template);
        let names: Vec<&str> = targets.iter().map(|(name, _)| *name).collect();
        println!("created {} in {}: {}", names.len(), args.dir.display(), names.join(", "));
    }
    Ok(())
}

/// Writes every file of `template` under `dir`. Refuses to overwrite an existing
/// file unless `force` is set; nothing is written when the check fails.
fn write_template(template: Template, dir: &Path, force: bool) -> Result<(), CliError> {
    let targets = files(template);
    if !force {
        if let Some((name, _)) = targets.iter().find(|(name, _)| dir.join(name).exists()) {
            let path = dir.join(name);
            return Err(input(format!(
                "'{}' already exists; use --force to overwrite it",
                path.display()
            )));
        }
    }
    fs::create_dir_all(dir)
        .map_err(|e| runtime(format!("cannot create '{}': {e}", dir.display())))?;
    for (name, text) in targets {
        let path = dir.join(name);
        fs::write(&path, text)
            .map_err(|e| runtime(format!("cannot write '{}': {e}", path.display())))?;
        log::info!("wrote {}", path.display());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn embedded_templates_are_not_empty() {
        for template in [Template::LjFluid, Template::HepAnalysis, Template::WaterBox] {
            for (name, text) in files(template) {
                assert!(!text.trim().is_empty(), "{name} is empty");
            }
        }
    }

    #[test]
    fn water_box_lists_every_molecule_and_no_placeholder() {
        let text = water_box_config();
        assert!(!text.contains("@MOLECULES@"));
        // First and last molecule of 512; no triple beyond the last molecule.
        assert!(text.contains("[0, 1, 2]"));
        assert!(text.contains("[1533, 1534, 1535]"));
        assert!(!text.contains("[1536"));
    }
}
