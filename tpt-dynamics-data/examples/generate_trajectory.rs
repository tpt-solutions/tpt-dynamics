//! Writes a synthetic trajectory in the chunked binary format, one frame at a
//! time, so generating a multi-gigabyte file needs only one frame in memory.
//!
//! Atoms start on a random lattice in a cubic periodic box and take a random
//! walk of at most `--step` nm per frame. The output is deterministic for a
//! given `--seed`.
//!
//! ```text
//! cargo run --release -p tpt-dynamics-data --example generate_trajectory -- \
//!     --atoms 100000 --frames 2000 --precision f32 --out data/traj-2gb.trj
//! ```
//!
//! Options: `--atoms N` (default 100000), `--frames F` (default 100),
//! `--precision f32|f64` (default f32), `--chunk C` frames per chunk (default 16),
//! `--box L` edge in nm (default 10), `--step S` nm per frame (default 0.01),
//! `--seed K` (default 1), `--out PATH` (required).

use std::fs::File;
use std::io::BufWriter;
use std::process::ExitCode;

use tpt_dynamics_data::rng::SplitMix64;
use tpt_dynamics_data::trajectory::{Frame, FrameSink, Precision, TrajectoryWriter};

struct Options {
    atoms: usize,
    frames: u64,
    precision: Precision,
    chunk: u32,
    edge: f64,
    step: f64,
    seed: u64,
    out: String,
}

fn parse_args() -> Result<Options, String> {
    let mut opts = Options {
        atoms: 100_000,
        frames: 100,
        precision: Precision::F32,
        chunk: 16,
        edge: 10.0,
        step: 0.01,
        seed: 1,
        out: String::new(),
    };
    let mut args = std::env::args().skip(1);
    while let Some(flag) = args.next() {
        let mut value = || args.next().ok_or_else(|| format!("{flag} needs a value"));
        match flag.as_str() {
            "--atoms" => opts.atoms = value()?.parse().map_err(|e| format!("--atoms: {e}"))?,
            "--frames" => opts.frames = value()?.parse().map_err(|e| format!("--frames: {e}"))?,
            "--chunk" => opts.chunk = value()?.parse().map_err(|e| format!("--chunk: {e}"))?,
            "--box" => opts.edge = value()?.parse().map_err(|e| format!("--box: {e}"))?,
            "--step" => opts.step = value()?.parse().map_err(|e| format!("--step: {e}"))?,
            "--seed" => opts.seed = value()?.parse().map_err(|e| format!("--seed: {e}"))?,
            "--out" => opts.out = value()?,
            "--precision" => {
                opts.precision = match value()?.as_str() {
                    "f32" => Precision::F32,
                    "f64" => Precision::F64,
                    other => return Err(format!("--precision must be f32 or f64, got '{other}'")),
                };
            }
            other => return Err(format!("unknown option '{other}'")),
        }
    }
    if opts.out.is_empty() {
        return Err("--out PATH is required".into());
    }
    if opts.atoms == 0 || opts.frames == 0 || opts.edge.is_nan() || opts.edge <= 0.0 {
        return Err("--atoms, --frames and --box must be positive".into());
    }
    Ok(opts)
}

fn main() -> ExitCode {
    let opts = match parse_args() {
        Ok(o) => o,
        Err(e) => {
            eprintln!("error: {e}");
            return ExitCode::FAILURE;
        }
    };
    match run(&opts) {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}

fn run(opts: &Options) -> Result<(), String> {
    let file = File::create(&opts.out).map_err(|e| format!("cannot create '{}': {e}", opts.out))?;
    let mut writer =
        TrajectoryWriter::new(BufWriter::new(file), opts.atoms, opts.precision, opts.chunk)
            .map_err(|e| e.to_string())?;

    let mut rng = SplitMix64::new(opts.seed);
    let mut positions: Vec<[f64; 3]> = (0..opts.atoms)
        .map(|i| {
            let side = (opts.atoms as f64).cbrt().ceil() as usize;
            let (x, y, z) = (i % side, (i / side) % side, i / (side * side));
            let cell = opts.edge / side as f64;
            [(x as f64 + 0.5) * cell, (y as f64 + 0.5) * cell, (z as f64 + 0.5) * cell]
        })
        .collect();

    for step in 0..opts.frames {
        if step > 0 {
            for p in &mut positions {
                for c in p.iter_mut() {
                    let d = rng.uniform(-opts.step, opts.step);
                    *c = (*c + d).rem_euclid(opts.edge);
                }
            }
        }
        let frame = Frame {
            step,
            time: step as f64 * 0.002,
            box_lengths: [opts.edge; 3],
            positions: positions.clone(),
        };
        writer.push_frame(&frame).map_err(|e| e.to_string())?;
    }
    let inner = writer.finish().map_err(|e| e.to_string())?;
    let file = inner.into_inner().map_err(|e| format!("flush failed: {}", e.error()))?;
    file.sync_all().map_err(|e| format!("sync failed: {e}"))?;
    let bytes = std::fs::metadata(&opts.out).map_or(0, |m| m.len());
    println!(
        "wrote {} frame(s) of {} atom(s), {bytes} bytes, to {}",
        opts.frames, opts.atoms, opts.out
    );
    Ok(())
}
