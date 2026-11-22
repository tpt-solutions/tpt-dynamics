//! `tptdyn init TEMPLATE [DIR]`: scaffolds a runnable project from templates that
//! are embedded in the binary.
//!
//! Every target file is checked before anything is written, so a refused init
//! leaves the directory untouched. `--force` replaces existing files.

use std::fs;
use std::path::Path;

use crate::cli::{InitArgs, Template};
use crate::error::{input, runtime, CliError};

/// One file of a template: its path inside the project and its embedded text.
type TemplateFile = (&'static str, &'static str);

const LJ_FLUID: &[TemplateFile] = &[
    ("simulate.toml", include_str!("../templates/lj-fluid/simulate.toml")),
    ("README.txt", include_str!("../templates/lj-fluid/README.txt")),
];

const HEP_ANALYSIS: &[TemplateFile] = &[
    ("make_events.py", include_str!("../templates/hep-analysis/make_events.py")),
    ("analyze.txt", include_str!("../templates/hep-analysis/analyze.txt")),
    ("README.txt", include_str!("../templates/hep-analysis/README.txt")),
];

/// The files of a template.
fn files(template: Template) -> &'static [TemplateFile] {
    match template {
        Template::LjFluid => LJ_FLUID,
        Template::HepAnalysis => HEP_ANALYSIS,
    }
}

/// Runs the `init` subcommand.
pub(crate) fn run(args: &InitArgs, quiet: bool) -> Result<(), CliError> {
    write_template(args.template, &args.dir, args.force)?;
    if !quiet {
        let names: Vec<&str> = files(args.template).iter().map(|(name, _)| *name).collect();
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
        for template in [Template::LjFluid, Template::HepAnalysis] {
            for (name, text) in files(template) {
                assert!(!text.trim().is_empty(), "{name} is empty");
            }
        }
    }
}
