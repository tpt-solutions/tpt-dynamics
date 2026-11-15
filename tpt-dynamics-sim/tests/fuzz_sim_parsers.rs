//! Fuzz-style robustness tests for the PDB and GRO readers.
//!
//! A deterministic xoshiro PRNG (seeded per test) generates random, mutated and
//! truncated inputs. The readers must return `Ok` or `Err` for each input and
//! never panic. Each test runs several thousand iterations, so a failure can be
//! reproduced by re-running the same test.

use tpt_dynamics_sim::gro::{read_gro, read_gro_with_limits, write_gro};
use tpt_dynamics_sim::pbc::SimBox;
use tpt_dynamics_sim::pdb::{read_pdb, read_pdb_with_limits, write_pdb};
use tpt_dynamics_sim::rng::Rng;
use tpt_dynamics_sim::solvate::place_water_grid;
use tpt_dynamics_sim::system::System;
use tpt_dynamics_sim::textio::ReadLimits;

/// Iterations per test. Override with `FUZZ_ITERS`.
fn iterations() -> usize {
    std::env::var("FUZZ_ITERS").ok().and_then(|v| v.parse().ok()).unwrap_or(4000)
}
/// Characters used to overwrite one byte of a valid GRO file.
const ALPHABET: &[u8] = b"0123456789-+. eEnaNi\n";

/// A valid 24-atom water system, rendered as PDB and as GRO.
fn valid_inputs() -> (Vec<u8>, Vec<u8>, usize) {
    let sys: System = place_water_grid(SimBox::new(2.0, 2.0, 2.0).unwrap(), 1.0).unwrap();
    let mut pdb_buf = Vec::new();
    write_pdb(&mut pdb_buf, &sys).unwrap();
    let mut gro_buf = Vec::new();
    write_gro(&mut gro_buf, &sys, "fuzz corpus").unwrap();
    (pdb_buf, gro_buf, sys.len())
}

fn random_bytes(rng: &mut Rng, max_len: u64) -> Vec<u8> {
    let len = (rng.next_u64() % max_len) as usize;
    (0..len).map(|_| (rng.next_u64() & 0xFF) as u8).collect()
}

/// Applies 1 to 8 random byte overwrites, insertions or deletions.
fn mutate(rng: &mut Rng, base: &[u8]) -> Vec<u8> {
    let mut out = base.to_vec();
    let edits = 1 + rng.next_u64() % 8;
    for _ in 0..edits {
        if out.is_empty() {
            out.push((rng.next_u64() & 0xFF) as u8);
            continue;
        }
        let at = (rng.next_u64() as usize) % out.len();
        match rng.next_u64() % 3 {
            0 => out[at] = (rng.next_u64() & 0xFF) as u8,
            1 => out.insert(at, (rng.next_u64() & 0xFF) as u8),
            _ => {
                out.remove(at);
            }
        }
    }
    out
}

#[test]
fn random_bytes_return_errors_and_never_panic() {
    let mut rng = Rng::seed_from_u64(0xF022_0001);
    for i in 0..iterations() {
        let bytes = random_bytes(&mut rng, 600);
        // Random bytes cannot form a valid structure with a box and atoms.
        assert!(read_pdb(&bytes[..]).is_err(), "PDB iteration {i} unexpectedly parsed");
        assert!(read_gro(&bytes[..]).is_err(), "GRO iteration {i} unexpectedly parsed");
    }
}

#[test]
fn mutated_valid_inputs_never_panic() {
    let (pdb_text, gro_text, n_atoms) = valid_inputs();
    let mut rng = Rng::seed_from_u64(0xF022_0002);
    let tight = ReadLimits { max_line_bytes: 60, max_atoms: 10 };
    for _ in 0..iterations() {
        let p = mutate(&mut rng, &pdb_text);
        if let Ok(sys) = read_pdb(&p[..]) {
            assert!(sys.len() <= n_atoms.max(1) + p.len());
        }
        let _ = read_pdb_with_limits(&p[..], tight);
        let g = mutate(&mut rng, &gro_text);
        let _ = read_gro(&g[..]);
        let _ = read_gro_with_limits(&g[..], tight);
    }
}

#[test]
fn truncated_valid_inputs_never_panic_and_gro_truncation_is_an_error() {
    let (pdb_text, gro_text, n_atoms) = valid_inputs();
    // The GRO box line is the last line; every cut before it leaves atoms missing.
    let box_line_start =
        gro_text[..gro_text.len() - 1].iter().rposition(|&b| b == b'\n').map_or(0, |p| p + 1);
    for n in 0..=gro_text.len() {
        let res = read_gro(&gro_text[..n]);
        if n < box_line_start {
            assert!(res.is_err(), "GRO cut at {n} parsed");
        }
        if let Ok(sys) = res {
            assert_eq!(sys.len(), n_atoms);
        }
    }
    for n in 0..=pdb_text.len() {
        if let Ok(sys) = read_pdb(&pdb_text[..n]) {
            assert!(sys.len() <= n_atoms);
        }
    }
}

#[test]
fn byte_flips_in_numeric_fields_are_rejected_or_parsed_without_panic() {
    let (_, gro_text, _) = valid_inputs();
    let mut rng = Rng::seed_from_u64(0xF022_0003);
    for _ in 0..iterations() {
        let mut g = gro_text.clone();
        let at = (rng.next_u64() as usize) % g.len();
        // Replace one byte with a digit, sign, dot, letter or space.
        g[at] = ALPHABET[(rng.next_u64() as usize) % ALPHABET.len()];
        let _ = read_gro(&g[..]);
    }
}
