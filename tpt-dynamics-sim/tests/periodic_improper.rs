//! Periodic improper torsions (`periodicimpropertype`), the form Amber uses.
//!
//! An improper matched by a periodic type uses its terms and is not counted as a
//! harmonic improper. An improper matched by no periodic type falls back to the
//! harmonic `impropertype`, as before.

use tpt_dynamics_core::Vec3;
use tpt_dynamics_sim::ffparams::ParamSet;
use tpt_dynamics_sim::pbc::SimBox;

const BASE: &str = r#"
name = "periodic-improper-test"
format_version = 1

[[atomtype]]
name = "C"
mass = 12.0
sigma = 0.34
epsilon = 0.36

[[atomtype]]
name = "H"
mass = 1.0
sigma = 0.25
epsilon = 0.1

[[bondtype]]
types = ["C", "C"]
k = 1000.0
r0 = 0.15

[[bondtype]]
types = ["C", "H"]
k = 1000.0
r0 = 0.11

[[angletype]]
types = ["H", "C", "H"]
k = 100.0
theta0 = 109.5

[[molecule]]
name = "CENTRE"
atoms = [
  { name = "C0", type = "C", charge = 0.0, position = [0.0, 0.0, 0.0] },
  { name = "A", type = "H", charge = 0.0, position = [0.11, 0.0, 0.0] },
  { name = "B", type = "H", charge = 0.0, position = [-0.05, 0.10, 0.0] },
  { name = "D", type = "H", charge = 0.0, position = [-0.05, -0.05, 0.09] },
]
bonds = [[0, 1], [0, 2], [0, 3]]
impropers = [[0, 1, 2, 3]]
"#;

fn with_periodic(periodic: &str) -> ParamSet {
    ParamSet::from_toml_str(&format!("{BASE}\n{periodic}")).expect("parameters parse")
}

fn placed(p: &ParamSet) -> tpt_dynamics_sim::system::System {
    p.assemble(&[("CENTRE", Vec3::new(1.0, 1.0, 1.0))], SimBox::new(4.0, 4.0, 4.0).unwrap())
        .expect("assemble")
        .system
}

#[test]
fn matched_improper_uses_periodic_terms() {
    let p = with_periodic(
        r#"
[[periodicimpropertype]]
types = ["C", "H", "H", "H"]
terms = [{ k = 4.0, n = 2, phase = 180.0 }]
"#,
    );
    let system = placed(&p);
    assert_eq!(system.topology.impropers.len(), 1, "one improper in the template");
    let terms = p.bonded_terms(&system).expect("bonded terms");
    assert!(terms.impropers.is_empty(), "periodic match must not produce a harmonic improper");
    let periodic: Vec<_> =
        terms.dihedrals.iter().filter(|d| d.atoms == system.topology.impropers[0]).collect();
    assert_eq!(periodic.len(), 1, "the improper must appear once, with periodic terms");
    assert_eq!(periodic[0].terms.len(), 1);
    assert_eq!(periodic[0].terms[0].n, 2);
    assert!((periodic[0].terms[0].phase - std::f64::consts::PI).abs() < 1e-12);
    assert!((periodic[0].terms[0].k - 4.0).abs() < 1e-12);
}

#[test]
fn unmatched_improper_falls_back_to_harmonic() {
    let p = with_periodic(
        r#"
[[impropertype]]
types = ["*", "*", "*", "*"]
k = 50.0
xi0 = 0.0
"#,
    );
    let system = placed(&p);
    let terms = p.bonded_terms(&system).expect("bonded terms");
    assert_eq!(terms.impropers.len(), 1);
    assert!(terms.dihedrals.iter().all(|d| d.atoms != system.topology.impropers[0]));
}

#[test]
fn periodic_improper_types_are_validated() {
    let bad = format!(
        "{BASE}\n[[periodicimpropertype]]\ntypes = [\"C\", \"H\", \"H\", \"XX\"]\nterms = [{{ k = 4.0, n = 2, phase = 180.0 }}]\n"
    );
    assert!(ParamSet::from_toml_str(&bad).is_err(), "unknown atom type must be rejected");
    let zero_n = format!(
        "{BASE}\n[[periodicimpropertype]]\ntypes = [\"C\", \"H\", \"H\", \"H\"]\nterms = [{{ k = 4.0, n = 0, phase = 180.0 }}]\n"
    );
    assert!(ParamSet::from_toml_str(&zero_n).is_err(), "multiplicity must be at least 1");
}
