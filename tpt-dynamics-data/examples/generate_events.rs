//! Synthetic collision-event CSV generator.
//!
//! Writes one row per particle:
//!
//! ```text
//! event,particle,pdg,charge,pt,eta,phi,mass
//! ```
//!
//! About 30% of events contain a dilepton pair reconstructing a Z-like resonance
//! (mass near 91.2 `GeV`); the rest are soft background particles. The output is
//! deterministic for a given `--seed`.
//!
//! Usage:
//!
//! ```text
//! cargo run --release -p tpt-dynamics-data --example generate_events -- \
//!     --events 10000 --seed 42 --out events.csv
//! ```
//!
//! Use `--out -` to write to standard output.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::process::ExitCode;

use tpt_dynamics_data::rng::SplitMix64;

const MUON_MASS: f64 = 0.105_658;
const Z_MASS: f64 = 91.1876;
const Z_WIDTH: f64 = 2.4952;
const SIGNAL_FRACTION: f64 = 0.3;

/// Particle species used for background tracks: (pdg id, charge, mass in `GeV`).
const BACKGROUND: [(i32, i32, f64); 5] = [
    (211, 1, 0.139_570),
    (-211, -1, 0.139_570),
    (321, 1, 0.493_677),
    (2212, 1, 0.938_272),
    (22, 0, 0.0),
];

struct Options {
    events: u64,
    seed: u64,
    out: String,
}

fn parse_args() -> Result<Options, String> {
    let mut opts = Options { events: 10_000, seed: 1, out: "events.csv".to_owned() };
    let mut args = std::env::args().skip(1);
    while let Some(flag) = args.next() {
        let mut value = || args.next().ok_or_else(|| format!("{flag} needs a value"));
        match flag.as_str() {
            "--events" => {
                opts.events = value()?.parse().map_err(|e| format!("--events: {e}"))?;
            }
            "--seed" => {
                opts.seed = value()?.parse().map_err(|e| format!("--seed: {e}"))?;
            }
            "--out" => opts.out = value()?,
            "-h" | "--help" => {
                return Err("usage: generate_events [--events N] [--seed S] [--out PATH|-]".into())
            }
            other => return Err(format!("unknown argument {other}")),
        }
    }
    Ok(opts)
}

/// One generated particle in detector coordinates.
struct Particle {
    pdg: i32,
    charge: i32,
    pt: f64,
    eta: f64,
    phi: f64,
    mass: f64,
}

/// Fill `out` with the particles of one event. Reuses the buffer across calls.
fn generate_event(rng: &mut SplitMix64, out: &mut Vec<Particle>) {
    out.clear();
    let phi0 = rng.uniform(-std::f64::consts::PI, std::f64::consts::PI);

    if rng.next_f64() < SIGNAL_FRACTION {
        // Dilepton pair: back-to-back in phi, equal pt, with a Breit-Wigner-like
        // smear of the pair mass and a small rapidity separation.
        let mass = Z_MASS + 0.5 * Z_WIDTH * rng.gaussian();
        let deta = rng.gaussian() * 0.6;
        let pt = mass / (2.0 * (deta * 0.5).cosh());
        let eta_c = rng.gaussian();
        let sign = if rng.below(2) == 0 { 1 } else { -1 };
        out.push(Particle {
            pdg: 13 * sign,
            charge: sign,
            pt,
            eta: eta_c + deta * 0.5,
            phi: wrap(phi0),
            mass: MUON_MASS,
        });
        out.push(Particle {
            pdg: -13 * sign,
            charge: -sign,
            pt,
            eta: eta_c - deta * 0.5,
            phi: wrap(phi0 + std::f64::consts::PI),
            mass: MUON_MASS,
        });
    }

    // Background: 1 to 12 soft particles with an exponentially falling pt spectrum.
    let n_bg = 1 + rng.below(12);
    for _ in 0..n_bg {
        let (pdg, charge, mass) = BACKGROUND[rng.below(BACKGROUND.len() as u64) as usize];
        out.push(Particle {
            pdg,
            charge,
            pt: 0.3 + rng.exponential(0.8),
            eta: rng.uniform(-2.5, 2.5),
            phi: rng.uniform(-std::f64::consts::PI, std::f64::consts::PI),
            mass,
        });
    }
}

fn wrap(phi: f64) -> f64 {
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut p = phi % two_pi;
    if p <= -std::f64::consts::PI {
        p += two_pi;
    } else if p > std::f64::consts::PI {
        p -= two_pi;
    }
    p
}

fn run(opts: &Options) -> io::Result<()> {
    let sink: Box<dyn Write> = if opts.out == "-" {
        Box::new(io::stdout().lock())
    } else {
        Box::new(File::create(&opts.out)?)
    };
    let mut w = BufWriter::with_capacity(1 << 20, sink);
    writeln!(w, "event,particle,pdg,charge,pt,eta,phi,mass")?;

    let mut rng = SplitMix64::new(opts.seed);
    let mut particles = Vec::with_capacity(16);
    let mut rows = 0u64;
    for event in 0..opts.events {
        generate_event(&mut rng, &mut particles);
        for (i, p) in particles.iter().enumerate() {
            writeln!(
                w,
                "{event},{i},{},{},{:.6},{:.6},{:.6},{:.6}",
                p.pdg, p.charge, p.pt, p.eta, p.phi, p.mass
            )?;
            rows += 1;
        }
    }
    w.flush()?;
    eprintln!("wrote {} events, {rows} particle rows (seed {})", opts.events, opts.seed);
    Ok(())
}

fn main() -> ExitCode {
    match parse_args() {
        Ok(opts) => match run(&opts) {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("error writing {}: {e}", opts.out);
                ExitCode::FAILURE
            }
        },
        Err(msg) => {
            eprintln!("{msg}");
            ExitCode::from(2)
        }
    }
}
