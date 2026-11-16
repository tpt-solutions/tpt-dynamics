//! Cross-check against OpenMM (todo Phase 5, validation).
//!
//! `validation/openmm_chain_reference.py` builds the same residue chains, with
//! its own geometry and topology code, and evaluates them with OpenMM. Its
//! outputs are the fixtures in `validation/fixtures/`: positions and the total
//! potential energy in kJ/mol. This test rebuilds each chain with
//! `ParamSet::assemble_chains` and checks that positions agree exactly and that
//! the potential agrees to a tight relative tolerance.
//!
//! The parameters are the synthetic test set (ours), not a published force field.
//! The check therefore tests the implementation (topology, bonded terms, exclusions,
//! 1-4 scaling, Lennard-Jones mixing), not the accuracy of any parameter set.
//!
//! Regenerate the fixtures with `python tpt-dynamics-sim/validation/openmm_chain_reference.py`.

use tpt_dynamics_core::Vec3;
use tpt_dynamics_sim::ffparams::{FfOptions, ParamSet};
use tpt_dynamics_sim::forcefield::{ForceField, Modifier};
use tpt_dynamics_sim::neighbour::VerletList;
use tpt_dynamics_sim::nonbonded::Electrostatics;
use tpt_dynamics_sim::pbc::SimBox;
use tpt_dynamics_sim::system::UnitSystem;

const SYNTHETIC: &str = include_str!("../data/synthetic-test.toml");
const TEMPLATES: &str = include_str!("../validation/chain-residues.toml");

/// Matches the constants in the Python script.
const BOND_LENGTH: f64 = 0.15;
const ORIGIN: Vec3 = Vec3::new(0.5, 0.5, 0.5);
const BOX: f64 = 12.0;
const CUTOFF: f64 = 4.0;

/// The fixture as written by the Python script.
struct Reference {
    potential: f64,
    positions: Vec<Vec3>,
    forces: Vec<Vec3>,
}

fn parse_fixture(text: &str) -> Reference {
    let mut lines = text.lines();
    let potential: f64 = lines
        .next()
        .and_then(|l| l.strip_prefix("potential "))
        .and_then(|v| v.trim().parse().ok())
        .expect("potential line");
    let _components = lines.next().expect("components line");
    let n: usize = lines
        .next()
        .and_then(|l| l.strip_prefix("atoms "))
        .and_then(|v| v.trim().parse().ok())
        .expect("atoms line");
    let parse_vec = |l: &str| -> Vec3 {
        let v: Vec<f64> = l.split_whitespace().map(|x| x.parse().expect("number")).collect();
        Vec3::new(v[0], v[1], v[2])
    };
    let positions: Vec<Vec3> = lines.by_ref().take(n).map(parse_vec).collect();
    assert_eq!(positions.len(), n, "fixture has fewer positions than atoms");
    assert_eq!(lines.next(), Some("forces"), "missing forces header");
    let forces: Vec<Vec3> = lines.by_ref().take(n).map(parse_vec).collect();
    assert_eq!(forces.len(), n, "fixture has fewer forces than atoms");
    Reference { potential, positions, forces }
}

fn params() -> ParamSet {
    ParamSet::from_toml_str(&format!("{SYNTHETIC}\n{TEMPLATES}")).expect("parameters parse")
}

fn check_case(name: &str, residues: &[&str], fixture: &str) {
    let reference = parse_fixture(fixture);
    let p = params();
    let asm = p
        .assemble_chains(
            &[(residues, ORIGIN)],
            BOND_LENGTH,
            SimBox::new(BOX, BOX, BOX).expect("box"),
        )
        .expect("assemble");
    let system = asm.system;

    // Positions: the two implementations place residues identically.
    assert_eq!(system.len(), reference.positions.len(), "{name}: atom count");
    let worst = system
        .positions
        .iter()
        .zip(&reference.positions)
        .map(|(a, b)| (*a - *b).norm())
        .fold(0.0_f64, f64::max);
    assert!(worst < 1e-12, "{name}: positions differ by {worst}");

    // Potential energy in kJ/mol, Lennard-Jones only (no charges in these templates).
    let opts = FfOptions {
        units: UnitSystem::Md,
        cutoff: CUTOFF,
        modifier: Modifier::None,
        electrostatics: Electrostatics::None,
    };
    let ff = p.force_field(&system, &opts).expect("force field");
    let mut neighbours = VerletList::new(CUTOFF, 0.2).expect("neighbour list");
    neighbours.rebuild(&system.positions, &system.sim_box).expect("rebuild");
    let mut forces = vec![Vec3::ZERO; system.len()];
    let ours = ff
        .compute(
            &system.positions,
            &system.atom_types,
            &system.charges,
            &system.sim_box,
            &neighbours,
            &mut forces,
        )
        .potential;

    // Forces, compared with the largest force as the scale.
    let scale = reference.forces.iter().map(|f| f.norm()).fold(1.0_f64, f64::max);
    let worst_force = forces
        .iter()
        .zip(&reference.forces)
        .map(|(a, b)| (*a - *b).norm())
        .fold(0.0_f64, f64::max);
    assert!(
        worst_force / scale < 1e-9,
        "{name}: forces differ by {worst_force} (largest reference force {scale})"
    );

    let relative = (ours - reference.potential).abs() / reference.potential.abs().max(1.0);
    assert!(
        relative < 1e-9,
        "{name}: potential {ours} vs OpenMM {} (relative {relative:e})",
        reference.potential
    );
}

#[test]
fn five_residue_rsa_chain_matches_openmm() {
    check_case("rsa5", &["RSA"; 5], include_str!("../validation/fixtures/rsa5.txt"));
}

#[test]
fn three_lipid_like_chain_matches_openmm() {
    check_case("lip3", &["LIP"; 3], include_str!("../validation/fixtures/lip3.txt"));
}

#[test]
fn mixed_chain_matches_openmm() {
    check_case(
        "mixed",
        &["RSA", "RSA", "LIP"],
        include_str!("../validation/fixtures/mixed.txt"),
    );
}
