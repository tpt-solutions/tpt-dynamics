//! End-to-end tests of the `tptdyn` binary: help and version, exit codes,
//! histogram and analyze output, and a short simulation.

use std::env;
use std::ffi::OsStr;
use std::fs::{self, File};
use std::path::{Path, PathBuf};
use std::process;

use tpt_dynamics_data::hist1d::Histogram1D;
use tpt_dynamics_data::lazy::{Column, DataFrame};
use tpt_dynamics_data::rowstore::write_rowstore;
use tpt_dynamics_data::trajectory::{Frame, FrameSink, Precision, TrajectoryWriter};
use tpt_dynamics_sim::output::read_trajectory;

/// A fresh scratch directory unique to this test process and name.
fn scratch(name: &str) -> PathBuf {
    let dir = env::temp_dir().join(format!("tptdyn-cli-test-{}-{name}", process::id()));
    let _ = fs::remove_dir_all(&dir);
    fs::create_dir_all(&dir).expect("create scratch dir");
    dir
}

/// Minimal stand-in for `assert_cmd`. The real crate enables `predicates/diff`,
/// which pulls in `normalize-line-endings` (Apache-2.0 only), so the test harness
/// is written locally to keep the dependency tree within the licence policy.
struct Cmd {
    inner: process::Command,
}

fn tptdyn() -> Cmd {
    Cmd { inner: process::Command::new(env!("CARGO_BIN_EXE_tptdyn")) }
}

impl Cmd {
    fn arg<S: AsRef<OsStr>>(&mut self, a: S) -> &mut Self {
        self.inner.arg(a);
        self
    }

    fn args<I, S>(&mut self, items: I) -> &mut Self
    where
        I: IntoIterator<Item = S>,
        S: AsRef<OsStr>,
    {
        self.inner.args(items);
        self
    }

    fn assert(&mut self) -> Assert {
        Assert { out: self.inner.output().expect("run tptdyn") }
    }
}

/// Output matcher: substring (from [`contains`]) or exact equality (from `&str`).
enum Matcher {
    Contains(String),
    Exact(String),
}

fn contains(s: &str) -> Matcher {
    Matcher::Contains(s.to_owned())
}

impl From<&str> for Matcher {
    fn from(s: &str) -> Self {
        Matcher::Exact(s.to_owned())
    }
}

impl Matcher {
    fn check(&self, text: &str, which: &str, out: &process::Output) {
        let ok = match self {
            Matcher::Contains(s) => text.contains(s.as_str()),
            Matcher::Exact(s) => text == s,
        };
        assert!(
            ok,
            "{which} did not match {}; got {text:?}; stderr: {}",
            match self {
                Matcher::Contains(s) => format!("contains {s:?}"),
                Matcher::Exact(s) => format!("== {s:?}"),
            },
            String::from_utf8_lossy(&out.stderr)
        );
    }
}

struct Assert {
    out: process::Output,
}

impl Assert {
    fn success(self) -> Self {
        assert!(
            self.out.status.success(),
            "expected success, got {:?}; stderr: {}",
            self.out.status,
            String::from_utf8_lossy(&self.out.stderr)
        );
        self
    }

    fn code(self, expected: i32) -> Self {
        assert_eq!(
            self.out.status.code(),
            Some(expected),
            "stderr: {}",
            String::from_utf8_lossy(&self.out.stderr)
        );
        self
    }

    fn stdout(self, m: impl Into<Matcher>) -> Self {
        let text = String::from_utf8_lossy(&self.out.stdout).into_owned();
        m.into().check(&text, "stdout", &self.out);
        self
    }

    fn stderr(self, m: impl Into<Matcher>) -> Self {
        let text = String::from_utf8_lossy(&self.out.stderr).into_owned();
        m.into().check(&text, "stderr", &self.out);
        self
    }
}

/// Small CSV: x values 0.5 1.5 1.7 2.2 -1 9 (one underflow, one overflow for [0,4)).
fn write_events(dir: &Path) -> PathBuf {
    let path = dir.join("events.csv");
    fs::write(&path, "x,y\n0.5,1\n1.5,2\n1.7,3\n2.2,4\n-1,5\n9,6\n").expect("write csv");
    path
}

#[test]
fn help_lists_subcommands_and_exit_codes() {
    tptdyn()
        .arg("--help")
        .assert()
        .success()
        .stdout(contains("simulate"))
        .stdout(contains("analyze"))
        .stdout(contains("histogram"))
        .stdout(contains("EXIT CODES"));
}

#[test]
fn version_prints_package_version() {
    tptdyn().arg("--version").assert().success().stdout(contains(env!("CARGO_PKG_VERSION")));
}

#[test]
fn unknown_option_is_usage_error() {
    tptdyn().arg("--bogus").assert().code(2);
}

#[test]
fn missing_subcommand_arguments_is_usage_error() {
    tptdyn().arg("histogram").assert().code(2);
}

#[test]
fn invalid_option_values_are_usage_errors() {
    let dir = scratch("usage");
    let csv = write_events(&dir);
    let out = dir.join("h.json");
    // Zero bins.
    tptdyn()
        .args(["histogram", csv.to_str().unwrap(), "--column", "x", "--bins", "0"])
        .args(["--range", "0,4", "--out", out.to_str().unwrap()])
        .assert()
        .code(2);
    // Unknown log level.
    tptdyn().args(["--log-level", "loud", "--version"]).assert().code(2);
    // Empty range.
    tptdyn()
        .args(["histogram", csv.to_str().unwrap(), "--column", "x", "--bins", "4"])
        .args(["--range", "4,0", "--out", out.to_str().unwrap()])
        .assert()
        .code(2);
    // Malformed filter.
    tptdyn().args(["analyze", csv.to_str().unwrap(), "--filter", "x ~ 1"]).assert().code(2);
}

#[test]
fn missing_files_are_input_errors() {
    let dir = scratch("missing");
    let nope = dir.join("nope.csv");
    tptdyn()
        .args(["histogram", nope.to_str().unwrap(), "--column", "x", "--bins", "4"])
        .args(["--range", "0,4", "--out", "h.json"])
        .assert()
        .code(3)
        .stderr(contains("cannot open"));
    tptdyn().args(["analyze", nope.to_str().unwrap()]).assert().code(3);
    tptdyn().args(["simulate", dir.join("nope.toml").to_str().unwrap()]).assert().code(3);
    tptdyn()
        .args(["--config", dir.join("nope.toml").to_str().unwrap(), "analyze", "x.csv"])
        .assert()
        .code(3);
}

#[test]
fn unknown_column_is_input_error() {
    let dir = scratch("badcol");
    let csv = write_events(&dir);
    tptdyn()
        .args(["histogram", csv.to_str().unwrap(), "--column", "zz", "--bins", "4"])
        .args(["--range", "0,4", "--out", dir.join("h.json").to_str().unwrap()])
        .assert()
        .code(3)
        .stderr(contains("no column named 'zz'"));
    tptdyn().args(["analyze", csv.to_str().unwrap(), "--filter", "zz > 1"]).assert().code(3);
}

#[test]
fn histogram_bin_counts_match_expected() {
    let dir = scratch("hist");
    let csv = write_events(&dir);
    let json = dir.join("h.json");
    tptdyn()
        .args(["--quiet", "histogram", csv.to_str().unwrap(), "--column", "x", "--bins", "4"])
        .args(["--range", "0,4", "--out", json.to_str().unwrap()])
        .assert()
        .success();
    let h = Histogram1D::from_json(&fs::read_to_string(&json).unwrap()).expect("valid json");
    assert_eq!(h.contents(), &[1.0, 2.0, 1.0, 0.0]);
    assert_eq!(h.underflow(), 1.0);
    assert_eq!(h.overflow(), 1.0);
    assert_eq!(h.entries(), 6);

    // Binary output, selected by the .bin extension, with --no-header and an index.
    let raw = dir.join("raw.csv");
    fs::write(&raw, "0.5,1\n1.5,2\n1.7,3\n2.2,4\n-1,5\n9,6\n").unwrap();
    let bin = dir.join("h.bin");
    tptdyn()
        .args(["--quiet", "histogram", raw.to_str().unwrap(), "--column", "0"])
        .args(["--bins", "4", "--range", "0,4", "--no-header", "--out", bin.to_str().unwrap()])
        .assert()
        .success();
    let hb = Histogram1D::from_bytes(&fs::read(&bin).unwrap()).expect("valid binary");
    assert_eq!(hb.contents(), &[1.0, 2.0, 1.0, 0.0]);
    assert_eq!(hb.underflow(), 1.0);
    assert_eq!(hb.overflow(), 1.0);
}

#[test]
fn histogram_reports_bad_numeric_field_with_line() {
    let dir = scratch("badnum");
    let csv = dir.join("bad.csv");
    fs::write(&csv, "x\n1.0\nabc\n").unwrap();
    tptdyn()
        .args(["histogram", csv.to_str().unwrap(), "--column", "x", "--bins", "4"])
        .args(["--range", "0,4", "--out", dir.join("h.json").to_str().unwrap()])
        .assert()
        .code(3)
        .stderr(contains("line 3"));
}

#[test]
fn analyze_filters_and_writes_histograms() {
    let dir = scratch("analyze");
    let csv = write_events(&dir);
    let out_dir = dir.join("out");
    tptdyn()
        .args(["analyze", csv.to_str().unwrap(), "--filter", "x > 1"])
        .args(["--hist", "x,2,0,4", "--out-dir", out_dir.to_str().unwrap()])
        .assert()
        .success()
        .stdout(contains("rows: 4 (of 6 read)"))
        .stdout(contains("column"));
    let h = Histogram1D::from_json(&fs::read_to_string(out_dir.join("x.json")).unwrap()).unwrap();
    // Rows with x > 1: 1.5, 1.7, 2.2, 9.
    assert_eq!(h.contents(), &[2.0, 1.0]);
    assert_eq!(h.overflow(), 1.0);
}

#[test]
fn analyze_summary_without_filters() {
    let dir = scratch("summary");
    let csv = write_events(&dir);
    tptdyn()
        .args(["analyze", csv.to_str().unwrap()])
        .assert()
        .success()
        .stdout(contains("rows: 6 (of 6 read)"));
}

#[test]
fn completions_and_man_pages() {
    tptdyn().args(["completions", "bash"]).assert().success().stdout(contains("tptdyn"));
    tptdyn().args(["completions", "nosuchshell"]).assert().code(2);

    let dir = scratch("man");
    tptdyn().args(["man", "--dir", dir.to_str().unwrap()]).assert().success();
    assert!(dir.join("tptdyn.1").is_file());
    assert!(dir.join("tptdyn-simulate.1").is_file());
}

/// A 27-atom lattice run for 50 steps with energy and trajectory output every 10 steps.
fn tiny_config(dir: &Path, dt: &str, extra: &str) -> PathBuf {
    let path = dir.join("tiny.toml");
    let text = format!(
        "seed = 7\n\
         dt = {dt}\n\
         nsteps = 50\n\
         temperature = 0.7\n\
         {extra}\n\
         [system]\n\
         cubic_cells = 3\n\
         spacing = 1.6\n\
         \n\
         [neighbour]\n\
         cutoff = 2.0\n\
         skin = 0.3\n\
         \n\
         [output]\n\
         energy_interval = 10\n\
         energy_file = \"energy.csv\"\n\
         trajectory_interval = 10\n\
         trajectory_file = \"traj.bin\"\n"
    );
    fs::write(&path, text).unwrap();
    path
}

#[test]
fn simulate_writes_energy_log_and_trajectory() {
    let dir = scratch("simulate");
    let cfg = tiny_config(&dir, "0.002", "");
    tptdyn().args(["simulate", cfg.to_str().unwrap(), "--quiet"]).assert().success();

    let energy = fs::read_to_string(dir.join("energy.csv")).unwrap();
    let lines: Vec<&str> = energy.lines().collect();
    assert_eq!(lines[0], "step,time,kinetic,potential,total,temperature,pressure");
    // Steps 0, 10, ..., 50.
    assert_eq!(lines.len(), 1 + 6, "energy log:\n{energy}");

    let traj = read_trajectory(File::open(dir.join("traj.bin")).unwrap()).unwrap();
    assert_eq!(traj.n_atoms, 27);
    assert_eq!(traj.frames.len(), 6);
}

/// `simulate` output must be readable by `analyze --trajectory` (both formats).
#[test]
fn simulate_output_round_trips_through_analyze() {
    for format in ["plain", "chunked"] {
        let dir = scratch(&format!("roundtrip-{format}"));
        let cfg = tiny_config(&dir, "0.002", "");
        if format == "chunked" {
            let text = fs::read_to_string(&cfg).unwrap();
            let text = text.replace(
                "trajectory_file = \"traj.bin\"",
                "trajectory_file = \"traj.bin\"\ntrajectory_format = \"chunked\"",
            );
            fs::write(&cfg, text).unwrap();
        }
        tptdyn().args(["simulate", cfg.to_str().unwrap(), "--quiet"]).assert().success();

        let out = dir.join("analysis.csv");
        tptdyn()
            .args(["analyze", "--trajectory", dir.join("traj.bin").to_str().unwrap()])
            .args(["--out", out.to_str().unwrap()])
            .assert()
            .success();
        let rows = read_traj_csv(&out);
        assert_eq!(rows.len(), 6, "{format}: one row per output frame");
        // Superposition rounding leaves a tiny residual, as in the rigid-motion test.
        assert!(rows[0].0.abs() < 1e-6, "{format}: first frame rmsd {}", rows[0].0);
    }
}

#[test]
fn simulate_output_dir_redirects_outputs() {
    let dir = scratch("simulate_outdir");
    let cfg = tiny_config(&dir, "0.002", "");
    let out = dir.join("results");
    tptdyn()
        .args(["simulate", cfg.to_str().unwrap(), "--output-dir", out.to_str().unwrap()])
        .args(["--quiet"])
        .assert()
        .success();
    assert!(out.join("energy.csv").is_file());
    assert!(out.join("traj.bin").is_file());
    assert!(!dir.join("energy.csv").exists());
}

#[test]
fn simulate_invalid_config_is_input_error() {
    let dir = scratch("simulate_bad");
    let cfg = tiny_config(&dir, "-1.0", "");
    tptdyn()
        .args(["simulate", cfg.to_str().unwrap(), "--quiet"])
        .assert()
        .code(3)
        .stderr(contains("dt must be"));

    let unknown = tiny_config(&dir, "0.002", "not_a_key = 1");
    tptdyn()
        .args(["simulate", unknown.to_str().unwrap(), "--quiet"])
        .assert()
        .code(3)
        .stderr(contains("not_a_key"));
}

#[test]
fn cli_defaults_file_is_applied_and_validated() {
    let dir = scratch("defaults");
    let defaults = dir.join("cli.toml");
    fs::write(&defaults, "quiet = true\nlog_level = \"error\"\n").unwrap();
    let csv = write_events(&dir);
    tptdyn()
        .args(["--config", defaults.to_str().unwrap(), "analyze", csv.to_str().unwrap()])
        .assert()
        .success()
        .stdout("");

    let bad = dir.join("bad.toml");
    fs::write(&bad, "colour = \"red\"\n").unwrap();
    tptdyn()
        .args(["--config", bad.to_str().unwrap(), "analyze", csv.to_str().unwrap()])
        .assert()
        .code(3);
}

// ---------------------------------------------------------------------------
// analyze --trajectory
// ---------------------------------------------------------------------------

/// Base shape of the four test atoms (nm).
const TRAJ_BASE: [[f64; 3]; 4] =
    [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 3.0]];

/// Writes a four-frame trajectory: the base shape, a translated copy, a copy rotated
/// 90 degrees about z, and a deformed copy (last atom moved by 1 nm in z).
fn write_test_trajectory(path: &Path) {
    let translated: Vec<[f64; 3]> =
        TRAJ_BASE.iter().map(|p| [p[0] + 5.0, p[1] - 2.0, p[2] + 1.0]).collect();
    let rotated: Vec<[f64; 3]> = TRAJ_BASE.iter().map(|p| [-p[1], p[0], p[2]]).collect();
    let mut deformed = TRAJ_BASE.to_vec();
    deformed[3][2] += 1.0;
    let frames = [TRAJ_BASE.to_vec(), translated, rotated, deformed];

    let file = File::create(path).expect("create trajectory");
    let mut writer = TrajectoryWriter::new(file, 4, Precision::F64, 2).expect("writer");
    for (i, positions) in frames.iter().enumerate() {
        let frame = Frame {
            step: 10 * i as u64,
            time: i as f64,
            box_lengths: [10.0, 10.0, 10.0],
            positions: positions.clone(),
        };
        writer.push_frame(&frame).expect("push frame");
    }
    writer.finish().expect("finish trajectory");
}

/// Writes a GRO structure with the given atom names at the base positions.
fn write_test_gro(path: &Path, names: &[&str]) {
    use std::fmt::Write as _;
    let mut text = String::from(
        "test
",
    );
    writeln!(text, "{:5}", names.len()).expect("format");
    for (i, name) in names.iter().enumerate() {
        let p = TRAJ_BASE[i];
        writeln!(
            text,
            "{:>5}{:<5}{:>5}{:>5}{:8.3}{:8.3}{:8.3}",
            1,
            "MOL",
            name,
            i + 1,
            p[0],
            p[1],
            p[2]
        )
        .expect("format");
    }
    writeln!(text, "   10.00000   10.00000   10.00000").expect("format");
    fs::write(path, text).expect("write gro");
}

/// Parses the CSV written by `analyze --trajectory` into (rmsd, rg) rows.
fn read_traj_csv(path: &Path) -> Vec<(f64, f64)> {
    let text = fs::read_to_string(path).expect("read csv");
    let mut lines = text.lines();
    assert_eq!(lines.next(), Some("step,time_ps,rmsd_nm,rg_nm"));
    lines
        .map(|line| {
            let cols: Vec<&str> = line.split(',').collect();
            assert_eq!(cols.len(), 4, "row: {line}");
            (cols[2].parse().expect("rmsd"), cols[3].parse().expect("rg"))
        })
        .collect()
}

/// Radius of gyration of `pos` with optional masses, computed directly.
fn expected_rg(pos: &[[f64; 3]], masses: Option<&[f64]>) -> f64 {
    let w = |i: usize| masses.map_or(1.0, |m| m[i]);
    let total: f64 = (0..pos.len()).map(w).sum();
    let mut c = [0.0; 3];
    for (i, p) in pos.iter().enumerate() {
        for k in 0..3 {
            c[k] += w(i) * p[k] / total;
        }
    }
    let s: f64 = pos
        .iter()
        .enumerate()
        .map(|(i, p)| w(i) * (0..3).map(|k| (p[k] - c[k]).powi(2)).sum::<f64>())
        .sum();
    (s / total).sqrt()
}

#[test]
fn analyze_trajectory_writes_rmsd_and_rg_per_frame() {
    let dir = scratch("traj-uniform");
    let traj = dir.join("run.trj");
    let out = dir.join("analysis.csv");
    write_test_trajectory(&traj);

    tptdyn()
        .args(["analyze", "--trajectory", traj.to_str().unwrap()])
        .args(["--out", out.to_str().unwrap()])
        .assert()
        .success();

    let rows = read_traj_csv(&out);
    assert_eq!(rows.len(), 4);
    let rg_base = expected_rg(&TRAJ_BASE, None);
    // Rigid translation and rotation: RMSD after superposition is zero.
    assert!(rows[0].0.abs() < 1e-6, "frame 0 rmsd {}", rows[0].0);
    assert!(rows[1].0.abs() < 1e-6, "translated rmsd {}", rows[1].0);
    assert!(rows[2].0.abs() < 1e-6, "rotated rmsd {}", rows[2].0);
    // The deformed frame moves one atom by 1 nm, so the RMSD is clearly positive.
    assert!(rows[3].0 > 0.1, "deformed rmsd {}", rows[3].0);
    // Rg is invariant under rigid motion; uniform masses give the geometric value.
    for (_, rg) in &rows[..3] {
        assert!((rg - rg_base).abs() < 1e-9, "rg {rg} vs {rg_base}");
    }
    let mut deformed = TRAJ_BASE.to_vec();
    deformed[3][2] += 1.0;
    assert!((rows[3].1 - expected_rg(&deformed, None)).abs() < 1e-9);
}

#[test]
fn analyze_trajectory_uses_structure_masses() {
    let dir = scratch("traj-masses");
    let traj = dir.join("run.trj");
    let gro = dir.join("system.gro");
    let out = dir.join("analysis.csv");
    write_test_trajectory(&traj);
    write_test_gro(&gro, &["C", "O", "N", "H"]);

    tptdyn()
        .args(["analyze", "--trajectory", traj.to_str().unwrap()])
        .args(["--structure", gro.to_str().unwrap()])
        .args(["--out", out.to_str().unwrap()])
        .args(["--quiet"])
        .assert()
        .success();

    let masses: Vec<f64> = ["C", "O", "N", "H"]
        .iter()
        .map(|e| tpt_dynamics_sim::system::element_mass(e).expect("element mass"))
        .collect();
    let rows = read_traj_csv(&out);
    let rg_mass = expected_rg(&TRAJ_BASE, Some(&masses));
    let rg_uniform = expected_rg(&TRAJ_BASE, None);
    assert!((rows[0].1 - rg_mass).abs() < 1e-9, "{} vs {rg_mass}", rows[0].1);
    assert!((rg_mass - rg_uniform).abs() > 1e-6, "masses should change Rg");
    // RMSD of a rigid copy stays zero with mass weights.
    assert!(rows[1].0.abs() < 1e-6);
}

#[test]
fn analyze_trajectory_usage_and_input_errors() {
    let dir = scratch("traj-errors");
    let traj = dir.join("run.trj");
    let csv = dir.join("events.csv");
    write_test_trajectory(&traj);
    fs::write(&csv, "x\n1\n").expect("write csv");
    let out = dir.join("analysis.csv");

    // Missing --out is a usage error.
    tptdyn().args(["analyze", "--trajectory", traj.to_str().unwrap()]).assert().code(2);
    // A CSV input cannot be combined with --trajectory.
    tptdyn()
        .args(["analyze", csv.to_str().unwrap(), "--trajectory", traj.to_str().unwrap()])
        .args(["--out", out.to_str().unwrap()])
        .assert()
        .code(2);
    // No input at all.
    tptdyn().args(["analyze"]).assert().code(2);
    // Structure file with an unsupported extension.
    let xyz = dir.join("system.xyz");
    fs::write(&xyz, "").expect("write xyz");
    tptdyn()
        .args(["analyze", "--trajectory", traj.to_str().unwrap()])
        .args(["--structure", xyz.to_str().unwrap(), "--out", out.to_str().unwrap()])
        .assert()
        .code(2);
    // Structure with the wrong atom count is an input error.
    let gro3 = dir.join("three.gro");
    write_test_gro(&gro3, &["C", "O", "N"]);
    tptdyn()
        .args(["analyze", "--trajectory", traj.to_str().unwrap()])
        .args(["--structure", gro3.to_str().unwrap(), "--out", out.to_str().unwrap()])
        .assert()
        .code(3);
    // A missing trajectory is an input error.
    tptdyn()
        .args(["analyze", "--trajectory", dir.join("nope.trj").to_str().unwrap()])
        .args(["--out", out.to_str().unwrap()])
        .assert()
        .code(3);
}

/// Path of a fixture file under `tests/fixtures`.
fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("tests").join("fixtures").join(name)
}

/// A small energy log with known statistics: T = 1, 2, 3 and E = 5, 4.9, 4.8.
fn write_energy_log(dir: &Path) -> PathBuf {
    let path = dir.join("energy.csv");
    fs::write(
        &path,
        "step,time,kinetic,potential,total,temperature,pressure\n\
0,0.0,10.0,-5.0,5.0,1.0,2.0\n\
10,0.02,10.0,-5.1,4.9,2.0,3.0\n\
20,0.04,10.0,-5.2,4.8,3.0,4.0\n",
    )
    .expect("write energy log");
    path
}

#[test]
fn init_scaffolds_templates_and_refuses_overwrite() {
    let dir = scratch("init");
    let lj = dir.join("lj");
    let lj_str = lj.to_str().unwrap();
    tptdyn().args(["init", "lj-fluid", lj_str, "--quiet"]).assert().success();
    let cfg = fs::read_to_string(lj.join("simulate.toml")).unwrap();
    assert!(cfg.contains("cubic_cells = 4"), "template content missing");
    assert!(lj.join("README.txt").exists());

    // A second init without --force is an input error and leaves files alone.
    tptdyn()
        .args(["init", "lj-fluid", lj_str, "--quiet"])
        .assert()
        .code(3)
        .stderr(contains("already exists; use --force"));
    tptdyn().args(["init", "lj-fluid", lj_str, "--force", "--quiet"]).assert().success();

    let hep = dir.join("hep");
    tptdyn().args(["init", "hep-analysis", hep.to_str().unwrap(), "--quiet"]).assert().success();
    let gen = fs::read_to_string(hep.join("make_events.py")).unwrap();
    assert!(gen.contains("def main"));
    assert!(fs::read_to_string(hep.join("analyze.txt")).unwrap().contains("--where"));
}

#[test]
fn init_refuses_when_any_target_exists() {
    let dir = scratch("init-partial");
    fs::write(dir.join("README.txt"), "mine").unwrap();
    tptdyn().args(["init", "hep-analysis", dir.to_str().unwrap(), "--quiet"]).assert().code(3);
    // Nothing else was written, and the existing file is untouched.
    assert!(!dir.join("make_events.py").exists());
    assert_eq!(fs::read_to_string(dir.join("README.txt")).unwrap(), "mine");
}

#[test]
fn validate_accepts_scaffolded_config_and_summarises() {
    let dir = scratch("validate-ok");
    tptdyn().args(["init", "lj-fluid", dir.to_str().unwrap(), "--quiet"]).assert().success();
    let cfg = dir.join("simulate.toml");
    tptdyn()
        .args(["validate", cfg.to_str().unwrap()])
        .assert()
        .success()
        .stdout(contains("config ok"))
        .stdout(contains("atoms         64"))
        .stdout(contains("2000"))
        .stdout(contains("MB"));
    assert!(!dir.join("energy.csv").exists(), "validate must not run the simulation");
}

#[test]
fn validate_suggests_closest_key_for_typo() {
    let dir = scratch("validate-typo");
    let cfg = tiny_config(&dir, "0.002", "");
    let text = fs::read_to_string(&cfg).unwrap().replace("energy_interval", "energy_intervall");
    fs::write(&cfg, text).unwrap();
    tptdyn()
        .args(["validate", cfg.to_str().unwrap()])
        .assert()
        .code(3)
        .stderr(contains("energy_intervall"))
        .stderr(contains("did you mean `energy_interval`?"));
}

#[test]
fn validate_missing_file_is_input_error() {
    let dir = scratch("validate-missing");
    let missing = dir.join("nope.toml");
    tptdyn()
        .args(["validate", missing.to_str().unwrap()])
        .assert()
        .code(3)
        .stderr(contains("cannot read config"));
}

#[test]
fn report_writes_self_contained_html_with_summary() {
    let dir = scratch("report");
    write_energy_log(&dir);
    tptdyn().args(["report", dir.to_str().unwrap(), "--quiet"]).assert().success();
    let html = fs::read_to_string(dir.join("report.html")).unwrap();
    assert!(html.contains("<svg"));
    assert!(html.contains("Mean temperature (K)</th><td>2.0000</td>"), "{html}");
    assert!(html.contains("Std dev of temperature (K)</th><td>1.0000</td>"));
    // Final drift: (4.8 - 5.0) / 5.0 = -4 %.
    assert!(html.contains("Energy drift, final (%)</th><td>-4.0000</td>"));
    assert!(!html.contains("<script"));
    assert!(!html.contains("<link"));
}

#[test]
fn report_reports_missing_energy_log_as_input_error() {
    let dir = scratch("report-missing");
    tptdyn()
        .args(["report", dir.to_str().unwrap()])
        .assert()
        .code(3)
        .stderr(contains("cannot read energy log"));
}

#[test]
fn simulate_report_flag_writes_report_next_to_energy_log() {
    let dir = scratch("sim-report");
    let cfg = tiny_config(&dir, "0.002", "");
    tptdyn().args(["simulate", cfg.to_str().unwrap(), "--quiet", "--report"]).assert().success();
    let html = fs::read_to_string(dir.join("report.html")).unwrap();
    assert!(html.contains("<svg"));
    assert!(html.contains("Energy samples</th><td>6</td>"), "{html}");
}

/// Runs `analyze` with the given options and checks that the printed row count appears.
fn assert_rows(csv: &Path, extra: &[&str], rows: &str) {
    tptdyn()
        .args(["analyze", csv.to_str().unwrap()])
        .args(extra)
        .assert()
        .success()
        .stdout(contains(rows));
}

#[test]
fn analyze_where_expressions_filter_rows() {
    let dir = scratch("where");
    let csv = write_events(&dir);
    // The data: (x, y) = (0.5,1) (1.5,2) (1.7,3) (2.2,4) (-1,5) (9,6).
    assert_rows(&csv, &["--where", "x > 1 and y < 5"], "rows: 3 (of 6 read)");
    assert_rows(&csv, &["--where", "x < 0 or y == 6"], "rows: 2 (of 6 read)");
    assert_rows(&csv, &["--where", "y * 0.5 > 1.2 and not (x < 0)"], "rows: 3 (of 6 read)");
    assert_rows(&csv, &["--where", "(x + 1) * 2 > 5 - y"], "rows: 4 (of 6 read)");
    // "and" binds tighter than "or": this is x < 0 or (x > 1 and y > 5).
    assert_rows(&csv, &["--where", "x < 0 or x > 1 and y > 5"], "rows: 2 (of 6 read)");
    // Division by zero is not an error: x / 0 is +inf for x > 0, which passes "> 1000".
    assert_rows(&csv, &["--where", "x / 0 > 1000"], "rows: 5 (of 6 read)");
    // 0 / 0 is NaN, and NaN compares false.
    assert_rows(&csv, &["--where", "0 / 0 > 1"], "rows: 0 (of 6 read)");
    // Filters and expressions are combined with AND.
    assert_rows(&csv, &["--filter", "y > 2", "--where", "x < 2"], "rows: 2 (of 6 read)");
}

#[test]
fn analyze_where_errors_have_documented_exit_codes() {
    let dir = scratch("where-errors");
    let csv = write_events(&dir);
    tptdyn()
        .args(["analyze", csv.to_str().unwrap(), "--where", "zz > 1"])
        .assert()
        .code(3)
        .stderr(contains("no column named 'zz'"));
    tptdyn()
        .args(["analyze", csv.to_str().unwrap(), "--where", "x > "])
        .assert()
        .code(2)
        .stderr(contains("--where 'x > '"));
    tptdyn()
        .args(["analyze", csv.to_str().unwrap(), "--where", "x = 1"])
        .assert()
        .code(2)
        .stderr(contains("use '==' for equality"));
}

#[test]
fn analyze_reads_row_store_and_parquet_inputs() {
    let dir = scratch("binary-inputs");
    let df = DataFrame::new([
        ("x", Column::F64(vec![0.5, 1.5, 1.7, 2.2, -1.0, 9.0])),
        ("y", Column::I64(vec![1, 2, 3, 4, 5, 6])),
    ])
    .expect("frame");
    let rows = dir.join("events.rows");
    write_rowstore(&rows, &df).expect("write row store");
    tptdyn()
        .args(["analyze", rows.to_str().unwrap(), "--where", "x > 1 and y < 5"])
        .assert()
        .success()
        .stdout(contains("rows: 3 (of 6 read)"));

    // A Parquet file from the data crate's fixtures, detected by its PAR1 magic.
    let parquet = fixture("flat_none_plain.parquet");
    tptdyn()
        .args(["analyze", parquet.to_str().unwrap()])
        .assert()
        .success()
        .stdout(contains("rows:"));
}

#[test]
fn analyze_rejects_trajectory_given_as_table() {
    let dir = scratch("analyze-traj");
    let traj = dir.join("traj.bin");
    write_test_trajectory(&traj);
    tptdyn()
        .args(["analyze", traj.to_str().unwrap()])
        .assert()
        .code(3)
        .stderr(contains("--trajectory"));
}
