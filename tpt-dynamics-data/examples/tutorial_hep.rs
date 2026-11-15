//! Tutorial: dilepton mass analysis with the data crate.
//!
//! This is the code behind `docs/tutorials/hep-analysis.md`. It
//!
//! 1. generates a small synthetic event sample as CSV text,
//! 2. reads it back with the streaming CSV reader,
//! 3. builds a dilepton invariant mass for each event,
//! 4. fills a histogram, merges two partial histograms, and
//! 5. writes the result in the binary histogram format.
//!
//! Run with:
//!
//! ```text
//! cargo run --release -p tpt-dynamics-data --example tutorial_hep
//! ```

use std::fmt::Write as _;
use std::io::BufReader;

use tpt_dynamics_data::csv::{CsvOptions, CsvReader};
use tpt_dynamics_data::hist1d::Histogram1D;
use tpt_dynamics_data::physics::pair_mass;
use tpt_dynamics_data::rng::SplitMix64;

const NEVENTS: u64 = 20_000;

/// Writes a CSV with one row per particle: `event,pt,eta,phi,mass,charge`.
/// Two opposite-charge muons are produced in about half of the events.
fn synthetic_csv(seed: u64) -> String {
    let mut rng = SplitMix64::new(seed);
    let mut text = String::from("event,pt,eta,phi,mass,charge\n");
    for event in 0..NEVENTS {
        let n_bkg = rng.below(4);
        for _ in 0..n_bkg {
            let pt = rng.exponential(5.0) + 0.5;
            let eta = rng.uniform(-2.5, 2.5);
            let phi = rng.uniform(-std::f64::consts::PI, std::f64::consts::PI);
            let _ = writeln!(text, "{event},{pt:.5},{eta:.5},{phi:.5},0.10566,0");
        }
        if rng.below(2) == 0 {
            // Back-to-back muons from a resonance near 91 GeV. For equal pT and
            // opposite azimuth, m = 2 pT, so pT = m / 2 puts the peak at m.
            let m_res = 91.0 + (rng.uniform(-1.0, 1.0) + rng.uniform(-1.0, 1.0)) * 2.0;
            let pt = m_res * 0.5 * (1.0 + rng.uniform(-0.01, 0.01));
            let eta = rng.uniform(-2.0, 2.0);
            let phi = rng.uniform(-std::f64::consts::PI, std::f64::consts::PI);
            let _ = writeln!(text, "{event},{pt:.5},{eta:.5},{phi:.5},0.10566,1");
            let _ = writeln!(
                text,
                "{event},{pt:.5},{:.5},{:.5},0.10566,-1",
                eta,
                phi + std::f64::consts::PI
            );
        }
    }
    text
}

/// Reads the CSV and returns the dilepton mass of every event that has an
/// opposite-charge pair. Events are grouped by the `event` column; this assumes
/// rows of one event are contiguous, as the generator writes them.
fn dilepton_masses(csv: &str) -> Result<Vec<f64>, Box<dyn std::error::Error>> {
    let reader = CsvReader::new(BufReader::new(csv.as_bytes()), CsvOptions::default())?;
    let mut masses = Vec::new();
    let mut current: Option<u64> = None;
    // (charge, (pt, eta, phi, m)) for each lepton in the current event.
    let mut leptons: Vec<Lepton> = Vec::new();

    for record in reader {
        let rec = record?;
        let event = rec.u64(0)?;
        if current != Some(event) {
            push_pair_mass(&leptons, &mut masses);
            leptons.clear();
            current = Some(event);
        }
        let charge = rec.i64(5)?;
        if charge != 0 {
            leptons.push((charge, (rec.f64(1)?, rec.f64(2)?, rec.f64(3)?, rec.f64(4)?)));
        }
    }
    push_pair_mass(&leptons, &mut masses);
    Ok(masses)
}

/// Adds the mass of the first opposite-charge pair in `leptons`, if there is one.
/// A lepton as `(charge, (pt, eta, phi, m))`.
type Lepton = (i64, (f64, f64, f64, f64));

fn push_pair_mass(leptons: &[Lepton], masses: &mut Vec<f64>) {
    let positive = leptons.iter().find(|(q, _)| *q > 0);
    let negative = leptons.iter().find(|(q, _)| *q < 0);
    if let (Some((_, a)), Some((_, b))) = (positive, negative) {
        masses.push(pair_mass(*a, *b));
    }
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let csv = synthetic_csv(2026);
    let masses = dilepton_masses(&csv)?;
    println!("events with a dilepton pair: {}", masses.len());

    // Two partial histograms, as two workers would produce them.
    let mut first = Histogram1D::uniform(120, 60.0, 120.0)?;
    let mut second = Histogram1D::uniform(120, 60.0, 120.0)?;
    for (i, &m) in masses.iter().enumerate() {
        if i % 2 == 0 {
            first.fill(m);
        } else {
            second.fill(m);
        }
    }

    // Merge in either order; the result does not depend on the order.
    let mut total = Histogram1D::uniform(120, 60.0, 120.0)?;
    total.merge(&second)?;
    total.merge(&first)?;

    let peak_bin =
        total.contents().iter().enumerate().max_by(|a, b| a.1.total_cmp(b.1)).map_or(0, |(i, _)| i);
    println!(
        "entries: {}, peak bin centre: {:.1} GeV, mean: {:.2} GeV",
        total.entries(),
        total.bin_center(peak_bin),
        total.mean().unwrap_or(f64::NAN)
    );

    // Persist as the checksummed binary format and read it back.
    let bytes = total.to_bytes();
    let restored = Histogram1D::from_bytes(&bytes)?;
    assert_eq!(restored, total);
    println!("serialised histogram: {} bytes", bytes.len());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tutorial_runs_and_finds_a_peak_near_z() {
        let csv = synthetic_csv(2026);
        let masses = dilepton_masses(&csv).unwrap();
        assert!(masses.len() > 5_000);
        let mut h = Histogram1D::uniform(120, 60.0, 120.0).unwrap();
        for m in masses {
            h.fill(m);
        }
        let peak = h.bin_center(
            h.contents().iter().enumerate().max_by(|a, b| a.1.total_cmp(b.1)).unwrap().0,
        );
        assert!((peak - 91.0).abs() < 5.0, "peak at {peak}");
    }
}
