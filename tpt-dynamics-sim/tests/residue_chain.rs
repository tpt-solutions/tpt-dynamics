//! Residue-chain assembly (todo Phase 5): linear chains of residue templates,
//! linked head to tail, with the topology generated across residue boundaries.
//!
//! The parameters are the synthetic test set shipped with the crate, extended
//! here with residue templates. They are not a published force field, so these
//! tests check the topology builder and the integrator. They do not check that
//! any real protein, lipid or nucleic acid is described correctly.

use tpt_dynamics_core::Vec3;
use tpt_dynamics_sim::ffparams::{FfOptions, ParamSet};
use tpt_dynamics_sim::forcefield::Modifier;
use tpt_dynamics_sim::integrator::{ForceEvaluator, VelocityVerlet};
use tpt_dynamics_sim::nonbonded::Electrostatics;
use tpt_dynamics_sim::pbc::SimBox;
use tpt_dynamics_sim::simulation::{maxwell_boltzmann_velocities, Simulation};
use tpt_dynamics_sim::system::UnitSystem;

const SYNTHETIC: &str = include_str!("../data/synthetic-test.toml");

/// Residue templates for the tests. `RSA` is a four-atom backbone-like unit,
/// `LIP` a six-atom unit standing in for a lipid tail, and `NUC` a three-atom
/// unit with a hydrogen-type end.
const RESIDUES: &str = r#"
[[molecule]]
name = "RSA"
atoms = [
  { name = "N",  type = "CT", charge = 0.0, position = [0.0,  0.0, 0.0] },
  { name = "CA", type = "CT", charge = 0.0, position = [0.12, 0.10, 0.0] },
  { name = "C",  type = "CT", charge = 0.0, position = [0.25, 0.1, 0.0] },
  { name = "O",  type = "CT", charge = 0.0, position = [0.38, 0.0, 0.0] },
]
bonds = [[0, 1], [1, 2], [2, 3]]
head = 0
tail = 3

[[molecule]]
name = "LIP"
atoms = [
  { name = "C1", type = "CT", charge = 0.0, position = [0.00, 0.00, 0.0] },
  { name = "C2", type = "CT", charge = 0.0, position = [0.13, 0.05, 0.0] },
  { name = "C3", type = "CT", charge = 0.0, position = [0.26, 0.00, 0.0] },
  { name = "C4", type = "CT", charge = 0.0, position = [0.39, 0.05, 0.0] },
  { name = "C5", type = "CT", charge = 0.0, position = [0.52, 0.00, 0.0] },
  { name = "C6", type = "CT", charge = 0.0, position = [0.65, 0.05, 0.0] },
]
bonds = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5]]
head = 0
tail = 5

[[molecule]]
name = "NUC"
atoms = [
  { name = "P",  type = "CT", charge = 0.0, position = [0.0,  0.0, 0.0] },
  { name = "S",  type = "CT", charge = 0.0, position = [0.14, 0.05, 0.0] },
  { name = "H1", type = "HC", charge = 0.0, position = [0.20, 0.12, 0.0] },
]
bonds = [[0, 1], [1, 2]]
head = 0
tail = 1
"#;

fn params() -> ParamSet {
    ParamSet::from_toml_str(&format!("{SYNTHETIC}\n{RESIDUES}")).expect("residue templates parse")
}

fn opts() -> FfOptions {
    FfOptions {
        units: UnitSystem::Reduced,
        cutoff: 1.0,
        modifier: Modifier::Shift,
        electrostatics: Electrostatics::None,
    }
}

#[test]
fn linear_chain_has_the_expected_topology_counts() {
    let p = params();
    let residues = ["RSA"; 5];
    let asm = p
        .assemble_chains(
            &[(&residues, Vec3::new(0.5, 0.5, 0.5))],
            0.15,
            SimBox::new(6.0, 6.0, 6.0).unwrap(),
        )
        .unwrap();
    let s = &asm.system;
    // 5 residues of 4 atoms, 4 links.
    assert_eq!(s.len(), 20);
    // 3 bonds per residue plus 4 links.
    assert_eq!(s.topology.bonds.len(), 19);
    // A path of n atoms has n-2 angles, n-3 proper dihedrals and n-3 1-4 pairs.
    assert_eq!(s.topology.angles.len(), 18);
    assert_eq!(s.topology.dihedrals.len(), 17);
    assert_eq!(s.topology.pairs14.len(), 17);
}

#[test]
fn links_are_placed_at_the_requested_bond_length_and_residues_are_numbered() {
    let p = params();
    let chain = ["RSA", "RSA", "LIP", "RSA"];
    let asm = p
        .assemble_chains(&[(&chain, Vec3::ZERO)], 0.15, SimBox::new(8.0, 8.0, 8.0).unwrap())
        .unwrap();
    let s = &asm.system;
    // Each link joins the tail of residue k to the head of residue k+1.
    let mut offset = 0;
    let mut link_lengths = Vec::new();
    let sizes = [4, 4, 6, 4];
    for (k, &size) in sizes.iter().enumerate().take(chain.len() - 1) {
        let tail = offset + size - 1;
        let head = offset + size;
        let d = s.positions[head] - s.positions[tail];
        link_lengths.push(d.norm());
        offset += sizes[k];
    }
    for d in &link_lengths {
        assert!((d - 0.15).abs() < 1e-12, "link length {d}");
    }
    // Residue ids count residues, not atoms.
    assert_eq!(s.residue_ids[0], 1);
    assert_eq!(s.residue_ids[4], 2);
    assert_eq!(s.residue_ids[8], 3);
    assert_eq!(s.residue_ids[14], 4);
}

#[test]
fn mixed_residues_in_several_chains_are_independent_molecules() {
    let p = params();
    let first = ["RSA", "RSA", "LIP"];
    let second = ["NUC", "NUC"];
    let asm = p
        .assemble_chains(
            &[(&first, Vec3::new(0.5, 0.5, 0.5)), (&second, Vec3::new(0.5, 3.0, 0.5))],
            0.15,
            SimBox::new(6.0, 6.0, 6.0).unwrap(),
        )
        .unwrap();
    let s = &asm.system;
    assert_eq!(s.len(), 4 + 4 + 6 + 3 + 3);
    // Template bonds plus one link per residue boundary.
    let template_bonds = 3 + 3 + 5 + 2 + 2;
    let links = 2 + 1;
    assert_eq!(s.topology.bonds.len(), template_bonds + links);
    // No bonded term may join the two chains.
    let chain_of = |a: usize| usize::from(a >= 14);
    for &(i, j) in &s.topology.bonds {
        assert_eq!(chain_of(i), chain_of(j), "bond {i}-{j} joins two chains");
    }
}

#[test]
fn unlinkable_residues_are_rejected() {
    let p = params();
    let box_ = SimBox::new(6.0, 6.0, 6.0).unwrap();
    // A residue template without a tail cannot be followed by another residue.
    let no_tail = format!(
        "{SYNTHETIC}\n{RESIDUES}\n[[molecule]]\nname = \"CAP\"\natoms = [{{ name = \"X\", type = \"CT\", charge = 0.0, position = [0.0, 0.0, 0.0] }}]\nhead = 0\n"
    );
    let p_no_tail = ParamSet::from_toml_str(&no_tail).unwrap();
    let err = p_no_tail.assemble_chains(&[(&["CAP", "RSA"], Vec3::ZERO)], 0.15, box_).unwrap_err();
    assert!(err.to_string().contains("no tail"), "got: {err}");

    // A residue after another must have a head.
    let no_head = format!(
        "{SYNTHETIC}\n{RESIDUES}\n[[molecule]]\nname = \"TAILONLY\"\natoms = [{{ name = \"X\", type = \"CT\", charge = 0.0, position = [0.0, 0.0, 0.0] }}]\ntail = 0\n"
    );
    let p_no_head = ParamSet::from_toml_str(&no_head).unwrap();
    let err =
        p_no_head.assemble_chains(&[(&["RSA", "TAILONLY"], Vec3::ZERO)], 0.15, box_).unwrap_err();
    assert!(err.to_string().contains("no head"), "got: {err}");

    assert!(p.assemble_chains(&[(&["NOPE"], Vec3::ZERO)], 0.15, box_).is_err());
    assert!(p.assemble_chains(&[(&["RSA"], Vec3::ZERO)], 0.0, box_).is_err());
    assert!(p.assemble_chains(&[(&[], Vec3::ZERO)], 0.15, box_).is_err());
}

#[test]
fn out_of_range_head_or_tail_is_rejected_at_load() {
    let bad = format!(
        "{SYNTHETIC}\n[[molecule]]\nname = \"BAD\"\natoms = [{{ name = \"X\", type = \"CT\", charge = 0.0, position = [0.0, 0.0, 0.0] }}]\ntail = 4\n"
    );
    assert!(ParamSet::from_toml_str(&bad).is_err());
}

#[test]
fn energy_is_conserved_for_a_residue_chain() {
    let p = params();
    let residues = ["RSA"; 5];
    let asm = p
        .assemble_chains(
            &[(&residues, Vec3::new(0.5, 0.5, 0.5))],
            0.15,
            SimBox::new(6.0, 6.0, 6.0).unwrap(),
        )
        .unwrap();
    let mut system = asm.system;
    maxwell_boltzmann_velocities(&mut system, 0.3, 4, UnitSystem::Reduced);
    let ff = p.force_field(&system, &opts()).unwrap();
    let eval = ForceEvaluator::new(Box::new(ff), 0.2).unwrap();
    let mut sim =
        Simulation::new(system, eval, Box::new(VelocityVerlet), None, 0.0005, UnitSystem::Reduced)
            .unwrap();
    let summary = sim.run(4000, &mut []).unwrap();
    println!("relative energy drift: {:.3e}", summary.relative_energy_drift());
    assert!(
        summary.relative_energy_drift().abs() < 2e-3,
        "drift {}",
        summary.relative_energy_drift()
    );
}

#[test]
fn straight_links_are_rejected() {
    // A residue whose bonds run along the link direction makes a 180 degree angle at
    // the link. The angle force is singular there, so the builder must refuse it.
    let straight = format!(
        "{SYNTHETIC}\n[[molecule]]\nname = \"LINE\"\natoms = [\n  {{ name = \"A\", type = \"CT\", charge = 0.0, position = [0.0, 0.0, 0.0] }},\n  {{ name = \"B\", type = \"CT\", charge = 0.0, position = [0.15, 0.0, 0.0] }},\n  {{ name = \"C\", type = \"CT\", charge = 0.0, position = [0.30, 0.0, 0.0] }},\n]\nbonds = [[0, 1], [1, 2]]\nhead = 0\ntail = 2\n"
    );
    let p = ParamSet::from_toml_str(&straight).unwrap();
    let err = p
        .assemble_chains(
            &[(&["LINE", "LINE"], Vec3::ZERO)],
            0.15,
            SimBox::new(6.0, 6.0, 6.0).unwrap(),
        )
        .unwrap_err();
    assert!(err.to_string().contains("nearly straight"), "got: {err}");
}
