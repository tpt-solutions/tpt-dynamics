//! Live progress line for `tptdyn simulate`, drawn on stderr.
//!
//! [`ProgressObserver`] is an [`Observer`] that runs at every step (it only reads
//! the step counter and the clock) and redraws the line at most about four times a
//! second. The formatting functions are pure so they can be tested with fixed inputs.

use std::io::{self, IsTerminal, Write};
use std::time::{Duration, Instant};

use tpt_dynamics_sim::output::{Observation, Observer};
use tpt_dynamics_sim::SimError;

/// Minimum time between two redraws.
const REDRAW_EVERY: Duration = Duration::from_millis(250);

/// Width used to blank out the previous line when it gets shorter.
const LINE_WIDTH: usize = 79;

/// Draws `step done/total`, percent, ETA and ns/day to stderr while a run proceeds.
pub(crate) struct ProgressObserver {
    total: u64,
    started: Instant,
    last_draw: Option<Instant>,
    drawn: bool,
}

impl ProgressObserver {
    /// Creates the observer for a run of `total` steps, starting now.
    pub(crate) fn new(total: u64) -> Self {
        Self { total, started: Instant::now(), last_draw: None, drawn: false }
    }
}

/// True when progress lines should be drawn: stderr is a terminal and not `--quiet`.
pub(crate) fn should_show(quiet: bool) -> bool {
    !quiet && io::stderr().is_terminal()
}

impl Observer for ProgressObserver {
    fn interval(&self) -> u64 {
        1
    }

    fn observe(&mut self, obs: &Observation<'_>) -> Result<(), SimError> {
        let done = obs.energies.step;
        let now = Instant::now();
        let due = self.last_draw.is_none_or(|t| now.duration_since(t) >= REDRAW_EVERY);
        if !(due || done >= self.total) {
            return Ok(());
        }
        self.last_draw = Some(now);
        let elapsed = now.duration_since(self.started).as_secs_f64();
        let line = progress_line(done, self.total, elapsed, obs.energies.time);
        // Progress is best effort: a failed terminal write must not stop the simulation.
        let mut err = io::stderr();
        let _ = write!(err, "\r{line:<LINE_WIDTH$}");
        let _ = err.flush();
        self.drawn = true;
        Ok(())
    }

    fn finish(&mut self) -> Result<(), SimError> {
        if self.drawn {
            let _ = writeln!(io::stderr());
        }
        Ok(())
    }
}

/// Estimated seconds left, from the average rate so far. `None` before any step is done.
pub(crate) fn eta_seconds(done: u64, total: u64, elapsed_s: f64) -> Option<f64> {
    if done == 0 || elapsed_s <= 0.0 || done >= total {
        return None;
    }
    let rate = done as f64 / elapsed_s;
    Some((total - done) as f64 / rate)
}

/// Formats an ETA as `42s`, `3m 05s` or `2h 03m`; `--` when unknown.
pub(crate) fn format_eta(seconds: Option<f64>) -> String {
    let Some(s) = seconds.filter(|s| s.is_finite() && *s >= 0.0) else {
        return "--".into();
    };
    let total = s.round() as u64;
    if total < 60 {
        format!("{total}s")
    } else if total < 3600 {
        format!("{}m {:02}s", total / 60, total % 60)
    } else {
        format!("{}h {:02}m", total / 3600, (total % 3600) / 60)
    }
}

/// Simulated nanoseconds per wall-clock day, from simulated time in ps and wall time in s.
/// Returns `None` when no wall time has passed.
pub(crate) fn ns_per_day(sim_time_ps: f64, wall_s: f64) -> Option<f64> {
    if wall_s <= 0.0 || !wall_s.is_finite() || !sim_time_ps.is_finite() {
        return None;
    }
    Some(sim_time_ps / 1000.0 / wall_s * 86_400.0)
}

/// Formats a throughput as `12.35 ns/day`; `--` when unknown.
pub(crate) fn format_ns_per_day(value: Option<f64>) -> String {
    match value {
        Some(v) if v.is_finite() => format!("{v:.2} ns/day"),
        _ => "-- ns/day".into(),
    }
}

/// The full progress line, without a trailing newline.
pub(crate) fn progress_line(done: u64, total: u64, elapsed_s: f64, sim_time_ps: f64) -> String {
    let pct = if total == 0 { 100.0 } else { 100.0 * done as f64 / total as f64 };
    let eta = format_eta(eta_seconds(done, total, elapsed_s));
    let speed = format_ns_per_day(ns_per_day(sim_time_ps, elapsed_s));
    format!("step {done}/{total} ({pct:>5.1}%)  ETA {eta:<9} {speed}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn eta_formats_seconds_minutes_and_hours() {
        assert_eq!(format_eta(Some(0.4)), "0s");
        assert_eq!(format_eta(Some(42.2)), "42s");
        assert_eq!(format_eta(Some(185.0)), "3m 05s");
        assert_eq!(format_eta(Some(7380.0)), "2h 03m");
        assert_eq!(format_eta(None), "--");
        assert_eq!(format_eta(Some(f64::INFINITY)), "--");
        assert_eq!(format_eta(Some(-1.0)), "--");
    }

    #[test]
    fn eta_is_remaining_steps_at_average_rate() {
        // 250 of 1000 steps in 10 s: 25 steps/s, so 750 steps take 30 s.
        let eta = eta_seconds(250, 1000, 10.0).unwrap();
        assert!((eta - 30.0).abs() < 1e-12);
        assert_eq!(eta_seconds(0, 1000, 10.0), None);
        assert_eq!(eta_seconds(1000, 1000, 10.0), None);
    }

    #[test]
    fn ns_per_day_from_fixed_inputs() {
        // 2 ps (0.002 ns) simulated in 1 s of wall time: 0.002 * 86400 = 172.8 ns/day.
        let v = ns_per_day(2.0, 1.0).unwrap();
        assert!((v - 172.8).abs() < 1e-9);
        assert_eq!(format_ns_per_day(Some(v)), "172.80 ns/day");
        assert_eq!(ns_per_day(2.0, 0.0), None);
        assert_eq!(format_ns_per_day(None), "-- ns/day");
    }

    #[test]
    fn progress_line_combines_fields() {
        // 500 of 1000 steps after 10 s with dt 0.002 ps: 1 ps simulated, 8640 ns/day.
        let line = progress_line(500, 1000, 10.0, 1.0);
        assert_eq!(line, "step 500/1000 ( 50.0%)  ETA 10s       8.64 ns/day");
    }
}
