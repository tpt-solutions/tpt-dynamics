//! `tptdyn report RUN_DIR` and the `--report` flag of `simulate`: a single,
//! self-contained `report.html` with inline SVG plots and a summary table.
//!
//! Input is the energy log written by `simulate` (header
//! `step,time,kinetic,potential,total,temperature,pressure`; columns are found by
//! name) and, optionally, a trajectory whose frame count is read from its header and
//! file size. The page has no scripts, external stylesheets or network references.

use std::fmt::Write as _;
use std::fs::{self, File};
use std::io::Read;
use std::path::Path;

use crate::cli::ReportArgs;
use crate::error::{input, runtime, CliError};

/// Pressure conversion from the engine's kJ mol⁻¹ nm⁻³ to bar.
const BAR_PER_ENGINE_PRESSURE: f64 = 16.6054;

/// Energy-log columns the report needs.
const REQUIRED_COLUMNS: [&str; 4] = ["time", "total", "temperature", "pressure"];

/// Columns of the energy log, one vector per column the report uses.
#[derive(Debug, Default, PartialEq)]
struct Energy {
    time: Vec<f64>,
    total: Vec<f64>,
    temperature: Vec<f64>,
    pressure: Vec<f64>,
}

/// Runs the `report` subcommand.
pub(crate) fn run(args: &ReportArgs, quiet: bool) -> Result<(), CliError> {
    let energy = args.run_dir.join("energy.csv");
    let trajectory = args.trajectory.as_ref().map(|t| {
        if t.is_absolute() {
            t.clone()
        } else {
            args.run_dir.join(t)
        }
    });
    let out = args.out.clone().unwrap_or_else(|| args.run_dir.join("report.html"));
    write_report(&energy, trajectory.as_deref(), &out)?;
    if !quiet {
        println!("report written to {}", out.display());
    }
    Ok(())
}

/// Reads the energy log (and the trajectory summary, when given) and writes the HTML
/// report to `out`.
///
/// # Errors
/// Input errors for a missing or malformed energy log; runtime errors when `out`
/// cannot be written.
pub(crate) fn write_report(
    energy_path: &Path,
    trajectory: Option<&Path>,
    out: &Path,
) -> Result<(), CliError> {
    let text = fs::read_to_string(energy_path)
        .map_err(|e| input(format!("cannot read energy log '{}': {e}", energy_path.display())))?;
    let energy = parse_energy(&text, energy_path)?;
    let traj_note = trajectory.map(trajectory_summary).transpose()?;
    let html = render_html(&energy, energy_path, traj_note.as_deref());
    fs::write(out, html).map_err(|e| runtime(format!("cannot write '{}': {e}", out.display())))
}

/// Parses the energy CSV, locating columns by name in the header row.
fn parse_energy(text: &str, path: &Path) -> Result<Energy, CliError> {
    let shown = path.display();
    let mut lines = text.lines().enumerate().filter(|(_, l)| !l.trim().is_empty());
    let (_, header) =
        lines.next().ok_or_else(|| input(format!("energy log '{shown}' is empty")))?;
    let names: Vec<&str> = header.split(',').map(str::trim).collect();
    let mut index = Vec::with_capacity(REQUIRED_COLUMNS.len());
    for want in REQUIRED_COLUMNS {
        let pos = names.iter().position(|n| *n == want).ok_or_else(|| {
            input(format!("energy log '{shown}' has no column '{want}' (header: {header})"))
        })?;
        index.push(pos);
    }

    let mut energy = Energy::default();
    for (line_no, line) in lines {
        let fields: Vec<&str> = line.split(',').map(str::trim).collect();
        let mut values = [0.0_f64; 4];
        for (slot, &col) in values.iter_mut().zip(&index) {
            let field = fields.get(col).copied().unwrap_or("");
            *slot = field.parse::<f64>().map_err(|_| {
                input(format!(
                    "energy log '{shown}' line {}: '{field}' is not a number",
                    line_no + 1
                ))
            })?;
        }
        energy.time.push(values[0]);
        energy.total.push(values[1]);
        energy.temperature.push(values[2]);
        energy.pressure.push(values[3]);
    }
    if energy.time.is_empty() {
        return Err(input(format!("energy log '{shown}' has no data rows")));
    }
    Ok(energy)
}

/// Mean and sample standard deviation (n - 1 denominator; 0 for a single value).
fn mean_and_std(values: &[f64]) -> (f64, f64) {
    let n = values.len() as f64;
    let mean = values.iter().sum::<f64>() / n;
    if values.len() < 2 {
        return (mean, 0.0);
    }
    let var = values.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0);
    (mean, var.sqrt())
}

/// Energy drift relative to the first sample, in percent of |E0|. `None` when E0 is 0.
fn drift_percent(total: &[f64]) -> Option<Vec<f64>> {
    let e0 = *total.first()?;
    if e0 == 0.0 || !e0.is_finite() {
        return None;
    }
    Some(total.iter().map(|e| (e - e0) / e0.abs() * 100.0).collect())
}

/// Summary of a plain trajectory (`TPTTRJ01`) from its header and file size.
fn trajectory_summary(path: &Path) -> Result<String, CliError> {
    let shown = path.display();
    let mut file =
        File::open(path).map_err(|e| input(format!("cannot open trajectory '{shown}': {e}")))?;
    let len = file.metadata().map_err(|e| input(format!("'{shown}': {e}")))?.len();
    let mut header = [0u8; 20];
    file.read_exact(&mut header)
        .map_err(|_| input(format!("trajectory '{shown}' is shorter than its header")))?;
    match &header[..8] {
        b"TPTTRJ01" => {
            let atoms = u64::from_le_bytes(header[12..20].try_into().unwrap_or([0; 8]));
            let frame_bytes = 8 + 8 + 24 + 12 * atoms;
            let frames = len.saturating_sub(20) / frame_bytes.max(1);
            Ok(format!("{frames} frames of {atoms} atoms (plain format)"))
        }
        b"TPTCTRJ1" => Ok("chunked format; frame count not summarised".to_string()),
        _ => Err(input(format!("'{shown}' is not a tptdyn trajectory"))),
    }
}

/// Escapes text for HTML and SVG content and attribute values.
fn esc(text: &str) -> String {
    let mut out = String::with_capacity(text.len());
    for c in text.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&#39;"),
            other => out.push(other),
        }
    }
    out
}

/// Axis tick label: fixed notation for ordinary magnitudes, scientific otherwise.
fn tick_label(v: f64) -> String {
    if v != 0.0 && (v.abs() >= 1.0e4 || v.abs() < 1.0e-3) {
        format!("{v:.2e}")
    } else {
        format!("{v:.3}")
    }
}

/// Inline SVG line plot of `ys` against `xs`. Non-finite points are skipped.
fn line_plot(title: &str, y_label: &str, xs: &[f64], ys: &[f64]) -> String {
    const W: f64 = 640.0;
    const H: f64 = 260.0;
    const LEFT: f64 = 78.0;
    const RIGHT: f64 = 16.0;
    const TOP: f64 = 16.0;
    const BOTTOM: f64 = 44.0;
    let points: Vec<(f64, f64)> = xs
        .iter()
        .zip(ys)
        .filter(|(x, y)| x.is_finite() && y.is_finite())
        .map(|(x, y)| (*x, *y))
        .collect();
    let mut svg = format!(
        "<svg viewBox=\"0 0 {W} {H}\" role=\"img\" aria-label=\"{}\" xmlns=\"http://www.w3.org/2000/svg\">\n",
        esc(title)
    );
    if points.is_empty() {
        svg.push_str("<text x=\"12\" y=\"24\">no finite data</text>\n</svg>");
        return svg;
    }
    let (mut x0, mut x1) = bounds(points.iter().map(|p| p.0));
    let (mut y0, mut y1) = bounds(points.iter().map(|p| p.1));
    if x1 <= x0 {
        x0 -= 0.5;
        x1 += 0.5;
    }
    if y1 <= y0 {
        y0 -= 0.5;
        y1 += 0.5;
    }
    let plot_w = W - LEFT - RIGHT;
    let plot_h = H - TOP - BOTTOM;
    let px = |x: f64| LEFT + (x - x0) / (x1 - x0) * plot_w;
    let py = |y: f64| TOP + plot_h - (y - y0) / (y1 - y0) * plot_h;

    svg.push_str("<g class=\"grid\">");
    for i in 0..=4 {
        let f = f64::from(i) / 4.0;
        let yv = y0 + f * (y1 - y0);
        let y = py(yv);
        let _ =
            writeln!(svg, "<line x1=\"{LEFT}\" x2=\"{}\" y1=\"{y:.2}\" y2=\"{y:.2}\"/>", W - RIGHT);
        let _ = writeln!(
            svg,
            "<text class=\"tick\" x=\"{:.2}\" y=\"{:.2}\" text-anchor=\"end\">{}</text>",
            LEFT - 6.0,
            y + 4.0,
            tick_label(yv)
        );
        let xv = x0 + f * (x1 - x0);
        let x = px(xv);
        let _ = writeln!(
            svg,
            "<line x1=\"{x:.2}\" x2=\"{x:.2}\" y1=\"{TOP}\" y2=\"{}\"/>",
            TOP + plot_h
        );
        let _ = writeln!(
            svg,
            "<text class=\"tick\" x=\"{x:.2}\" y=\"{:.2}\" text-anchor=\"middle\">{}</text>",
            TOP + plot_h + 16.0,
            tick_label(xv)
        );
    }
    svg.push_str("</g>\n");
    let _ = writeln!(
        svg,
        "<text class=\"label\" x=\"{:.2}\" y=\"{}\" text-anchor=\"middle\">time (ps)</text>",
        LEFT + plot_w / 2.0,
        H - 6.0
    );
    let _ = writeln!(svg,
        "<text class=\"label\" transform=\"translate(14 {:.2}) rotate(-90)\" text-anchor=\"middle\">{}</text>",
        TOP + plot_h / 2.0,
        esc(y_label)
    );
    let path: Vec<String> =
        points.iter().map(|(x, y)| format!("{:.2},{:.2}", px(*x), py(*y))).collect();
    let _ = writeln!(svg, "<polyline class=\"series\" points=\"{}\"/>", path.join(" "));
    svg.push_str("</svg>");
    svg
}

/// Minimum and maximum of a non-empty iterator (0, 0 when empty).
fn bounds(values: impl Iterator<Item = f64>) -> (f64, f64) {
    values.fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), v| (lo.min(v), hi.max(v)))
}

/// Builds the complete HTML document.
fn render_html(energy: &Energy, source: &Path, traj_note: Option<&str>) -> String {
    let t = &energy.time;
    let (t_mean, t_std) = mean_and_std(&energy.temperature);
    let (p_mean, p_std) = mean_and_std(&energy.pressure);
    let e0 = energy.total[0];
    let e_last = *energy.total.last().unwrap_or(&e0);
    let drift_final = if e0 == 0.0 {
        "n/a (initial energy is zero)".to_string()
    } else {
        format!("{:.4}", (e_last - e0) / e0.abs() * 100.0)
    };
    let drift_max = match drift_percent(&energy.total) {
        Some(d) => format!(
            "{:.4}",
            d.iter().fold(0.0_f64, |m, v| if v.is_finite() { m.max(v.abs()) } else { m })
        ),
        None => "n/a".to_string(),
    };

    let (drift_series, drift_label) = match drift_percent(&energy.total) {
        Some(d) => (d, "energy drift (% of initial)"),
        None => (energy.total.iter().map(|e| e - e0).collect(), "energy drift (kJ/mol)"),
    };

    let mut rows: Vec<(String, String)> = vec![
        ("Energy samples".into(), energy.time.len().to_string()),
        ("Time span (ps)".into(), format!("{:.4} to {:.4}", t[0], t[t.len() - 1])),
        ("Mean temperature (K)".into(), format!("{t_mean:.4}")),
        ("Std dev of temperature (K)".into(), format!("{t_std:.4}")),
        ("Mean pressure (kJ/mol/nm^3)".into(), format!("{p_mean:.4}")),
        ("Mean pressure (bar)".into(), format!("{:.4}", p_mean * BAR_PER_ENGINE_PRESSURE)),
        ("Std dev of pressure (kJ/mol/nm^3)".into(), format!("{p_std:.4}")),
        ("Initial total energy (kJ/mol)".into(), format!("{e0:.4}")),
        ("Final total energy (kJ/mol)".into(), format!("{e_last:.4}")),
        ("Energy drift, final (%)".into(), drift_final),
        ("Energy drift, largest excursion (%)".into(), drift_max),
    ];
    if let Some(note) = traj_note {
        rows.push(("Trajectory".into(), note.to_string()));
    }
    let mut table = String::new();
    for (k, v) in &rows {
        let _ = writeln!(table, "<tr><th scope=\"row\">{}</th><td>{}</td></tr>", esc(k), esc(v));
    }

    let plots = [
        line_plot("Total energy drift", drift_label, t, &drift_series),
        line_plot("Temperature", "temperature (K)", t, &energy.temperature),
        line_plot("Pressure", "pressure (kJ/mol/nm^3)", t, &energy.pressure),
    ]
    .join("\n");

    format!(
        "<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n\
<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n\
<title>tptdyn run report</title>\n<style>\n{STYLE}</style>\n</head>\n<body>\n\
<h1>tptdyn run report</h1>\n<p class=\"source\">Source: {}</p>\n\
<h2>Summary</h2>\n<table>\n{table}</table>\n\
<h2>Plots</h2>\n<div class=\"plots\">\n{plots}\n</div>\n\
<p class=\"foot\">Means and standard deviations are over the energy log samples. \
Drift is measured against the first sample.</p>\n</body>\n</html>\n",
        esc(&source.display().to_string()),
    )
}

/// Page styles: light and dark colour tokens, no external resources.
const STYLE: &str = "\
:root { color-scheme: light dark; --bg: #ffffff; --fg: #1d1f23; --muted: #5b616e; \
--line: #d7dbe2; --accent: #2b6cb0; }
@media (prefers-color-scheme: dark) {
  :root:not([data-theme=\"light\"]) { --bg: #16181d; --fg: #e6e8ec; --muted: #9aa3b2; \
--line: #343a46; --accent: #7fb2e8; }
}
:root[data-theme=\"dark\"] { --bg: #16181d; --fg: #e6e8ec; --muted: #9aa3b2; \
--line: #343a46; --accent: #7fb2e8; }
body { margin: 0; padding: 16px; background: var(--bg); color: var(--fg); \
font: 15px/1.5 system-ui, sans-serif; }
h1 { font-size: 1.5rem; margin: 0.4em 0; }
h2 { font-size: 1.15rem; margin: 1.2em 0 0.4em; }
.source, .foot { color: var(--muted); font-size: 0.9rem; }
table { border-collapse: collapse; max-width: 640px; width: 100%; }
th, td { text-align: left; padding: 4px 10px 4px 0; border-bottom: 1px solid var(--line); }
td { font-variant-numeric: tabular-nums; }
.plots { display: grid; gap: 12px; max-width: 720px; }
svg { width: 100%; height: auto; color: var(--fg); }
svg .grid line { stroke: var(--line); stroke-width: 1; }
svg .tick, svg .label { fill: var(--muted); font-size: 11px; }
svg .label { font-size: 12px; }
svg .series { fill: none; stroke: var(--accent); stroke-width: 1.5; }
";

#[cfg(test)]
mod tests {
    use super::*;

    fn sample() -> &'static str {
        "step,time,kinetic,potential,total,temperature,pressure\n\
0,0.0,10.0,-5.0,5.0,1.0,2.0\n\
10,0.02,10.0,-5.1,4.9,2.0,3.0\n\
20,0.04,10.0,-5.2,4.8,3.0,4.0\n"
    }

    #[test]
    fn parses_columns_by_name() {
        let e = parse_energy(sample(), Path::new("energy.csv")).unwrap();
        assert_eq!(e.time, vec![0.0, 0.02, 0.04]);
        assert_eq!(e.temperature, vec![1.0, 2.0, 3.0]);
        assert_eq!(e.total.last().copied(), Some(4.8));
    }

    #[test]
    fn missing_column_is_input_error() {
        let err = parse_energy("step,time,total\n0,0,1\n", Path::new("e.csv")).unwrap_err();
        assert!(matches!(err, CliError::Input(ref m) if m.contains("'temperature'")), "{err}");
    }

    #[test]
    fn mean_and_sample_std() {
        let (m, s) = mean_and_std(&[1.0, 2.0, 3.0]);
        assert!((m - 2.0).abs() < 1e-12);
        assert!((s - 1.0).abs() < 1e-12);
        assert_eq!(mean_and_std(&[4.0]), (4.0, 0.0));
    }

    #[test]
    fn drift_is_relative_to_first_sample() {
        let d = drift_percent(&[5.0, 4.9, 4.8]).unwrap();
        assert!((d[2] + 4.0).abs() < 1e-9);
        assert_eq!(drift_percent(&[0.0, 1.0]), None);
    }

    #[test]
    fn plot_is_svg_and_skips_non_finite_points() {
        let svg = line_plot("T", "K", &[0.0, 1.0, 2.0], &[1.0, f64::NAN, 3.0]);
        assert!(svg.starts_with("<svg"));
        assert!(svg.contains("points=\"78.00,"));
        assert!(line_plot("T", "K", &[0.0], &[f64::NAN]).contains("no finite data"));
    }

    #[test]
    fn html_escapes_and_includes_summary() {
        let e = parse_energy(sample(), Path::new("energy.csv")).unwrap();
        let html = render_html(&e, Path::new("a<b>/energy.csv"), None);
        assert!(html.contains("<svg"));
        assert!(html.contains("a&lt;b&gt;/energy.csv"));
        assert!(html.contains("<td>2.0000</td>"), "mean temperature row missing");
        assert!(!html.contains("<script"));
        assert!(!html.contains("<link"));
    }
}
