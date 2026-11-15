//! Robustness tests for the TOML inputs: the force-field parameter file
//! (`ffparams`) and the simulation config. Each input must parse or return an
//! error, never panic, and a valid input must survive truncation and mutation
//! without a panic.
//!
//! The iteration count is `FUZZ_ITERS` (default 4000). For a longer run:
//! `FUZZ_ITERS=200000 cargo test -p tpt-dynamics-sim --release --test fuzz_config`.

use tpt_dynamics_sim::config::SimConfig;
use tpt_dynamics_sim::ffparams::ParamSet;
use tpt_dynamics_sim::rng::Rng;

fn iterations() -> usize {
    std::env::var("FUZZ_ITERS").ok().and_then(|v| v.parse().ok()).unwrap_or(4000)
}

const FFPARAMS: &str = include_str!("../data/synthetic-test.toml");

const CONFIG: &str = r#"
units = "reduced"
seed = 1
dt = 0.004
nsteps = 100
temperature = 0.72
integrator = "velocity-verlet"

[neighbour]
cutoff = 2.5
skin = 0.3

[lennard_jones]
epsilon = 1.0
sigma = 1.0
modifier = "shift"

[output]
energy_interval = 10
energy_file = "e.csv"
"#;

/// Characters that appear in TOML documents, so random text reaches deeper
/// parser states than random bytes do.
const ALPHABET: &[u8] = b"[]{}=\"'\n .,#-_+0123456789eE:tfnu()abcdefghijklmnopqrstuvwxyz";

fn random_text(rng: &mut Rng, max_len: u64) -> String {
    let len = (rng.next_u64() % max_len) as usize;
    (0..len).map(|_| ALPHABET[(rng.next_u64() % ALPHABET.len() as u64) as usize] as char).collect()
}

/// Applies a few random edits to `base`: byte overwrites, insertions of TOML
/// punctuation, deletions, and large numbers written over digit runs.
fn mutate(rng: &mut Rng, base: &str) -> String {
    let mut bytes = base.as_bytes().to_vec();
    let edits = 1 + rng.next_u64() % 6;
    for _ in 0..edits {
        if bytes.is_empty() {
            break;
        }
        let at = (rng.next_u64() % bytes.len() as u64) as usize;
        match rng.next_u64() % 4 {
            0 => bytes[at] = ALPHABET[(rng.next_u64() % ALPHABET.len() as u64) as usize],
            1 => bytes.insert(at, ALPHABET[(rng.next_u64() % ALPHABET.len() as u64) as usize]),
            2 => {
                bytes.remove(at);
            }
            _ => {
                let big = b"99999999999999999999999";
                let end = (at + big.len()).min(bytes.len());
                bytes.splice(at..end, big.iter().copied());
            }
        }
    }
    // Mutations may break UTF-8 only inside multi-byte characters; keep valid text.
    String::from_utf8_lossy(&bytes).into_owned()
}

fn exercise_paramset(text: &str) {
    let _ = ParamSet::from_toml_str(text);
}

fn exercise_config(text: &str) {
    let _ = SimConfig::from_toml_str(text);
}

#[test]
fn valid_inputs_parse() {
    ParamSet::from_toml_str(FFPARAMS).expect("shipped parameter file parses");
    SimConfig::from_toml_str(CONFIG).expect("sample config parses");
}

#[test]
fn paramset_truncations_never_panic() {
    for cut in 0..FFPARAMS.len() {
        if FFPARAMS.is_char_boundary(cut) {
            exercise_paramset(&FFPARAMS[..cut]);
        }
    }
}

#[test]
fn config_truncations_never_panic() {
    for cut in 0..CONFIG.len() {
        if CONFIG.is_char_boundary(cut) {
            exercise_config(&CONFIG[..cut]);
        }
    }
}

#[test]
fn random_text_never_panics() {
    let mut rng = Rng::seed_from_u64(0xF00D);
    for _ in 0..iterations() {
        let text = random_text(&mut rng, 300);
        exercise_paramset(&text);
        exercise_config(&text);
    }
}

#[test]
fn mutated_valid_inputs_never_panic() {
    let mut rng = Rng::seed_from_u64(0xBEEF);
    for _ in 0..iterations() {
        exercise_paramset(&mutate(&mut rng, FFPARAMS));
        exercise_config(&mutate(&mut rng, CONFIG));
    }
}
