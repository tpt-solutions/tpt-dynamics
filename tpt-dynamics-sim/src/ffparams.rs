//! Force-field parameter files (TOML), molecule templates and the topology
//! builder.
//!
//! The format describes atom types, bond, angle, dihedral and improper types,
//! and molecule templates (atoms with types, charges and reference coordinates,
//! bonds and explicit impropers). Angles, proper dihedrals and 1-4 pairs are
//! generated from the bond graph. Parameters are looked up by type name; `"*"` is
//! a wildcard allowed in dihedral and improper patterns.
//!
//! Units: lengths in nm, energies in kJ/mol, masses in amu, charges in e, and
//! angles in **degrees** in the file (converted to radians on load). Dihedral
//! multiplicities are integers. Lennard-Jones parameters are `sigma` (nm) and
//! `epsilon` (kJ/mol) per atom type, combined with Lorentz-Berthelot rules.
//!
//! Example (abbreviated):
//!
//! ```toml
//! name = "my-field"
//! format_version = 1
//! fudge_lj = 0.5          # 1-4 Lennard-Jones scale
//! fudge_qq = 0.8333       # 1-4 Coulomb scale
//!
//! [[atomtype]]
//! name = "CT"
//! mass = 12.011
//! sigma = 0.34
//! epsilon = 0.36
//!
//! [[bondtype]]
//! types = ["CT", "CT"]
//! k = 250000.0            # kJ mol⁻¹ nm⁻²  (V = ½ k (r − r0)²)
//! r0 = 0.153
//!
//! [[angletype]]
//! types = ["CT", "CT", "CT"]
//! k = 500.0               # kJ mol⁻¹ rad⁻²
//! theta0 = 112.0         # degrees
//!
//! [[dihedraltype]]
//! types = ["*", "CT", "CT", "*"]
//! terms = [ { k = 1.2, n = 3, phase = 0.0 } ]   # phase in degrees
//!
//! [[molecule]]
//! name = "ethane"
//! atoms = [
//!   { name = "C1", type = "CT", charge = 0.0, position = [0.0, 0.0, 0.0] },
//!   { name = "C2", type = "CT", charge = 0.0, position = [0.153, 0.0, 0.0] },
//! ]
//! bonds = [[0, 1]]
//! ```
//!
//! The `[[angletype]]` table may also carry `urey_bradley_k` (kJ mol⁻¹ nm⁻²)
//! and `urey_bradley_r0` (nm) for a Urey-Bradley 1-3 term.
//!
//! The crate ships only a synthetic test file (`data/synthetic-test.toml`),
//! written for these tests and not derived from any published parameter set.
//! Published AMBER, CHARMM or GROMOS parameter data is not bundled.

use serde::Deserialize;
use tpt_dynamics_core::Vec3;

use crate::bonded::CombinedForceField;
use crate::bonded::{
    BondedForce, BondedTerms, DihedralTerm, HarmonicAngle, HarmonicBond, HarmonicImproper,
    PeriodicDihedral, UreyBradley,
};
use crate::error::SimError;
use crate::forcefield::{ForceField, LjType, Modifier};
use crate::nonbonded::{Electrostatics, Nonbonded, NonbondedSettings, COULOMB_FACTOR_MD};
use crate::pbc::SimBox;
use crate::system::{NewAtom, System, UnitSystem};
use crate::topology::Topology;

/// Supported file format version.
pub const FORMAT_VERSION: u32 = 1;

fn default_fudge_lj() -> f64 {
    0.5
}

fn default_fudge_qq() -> f64 {
    0.8333
}

/// One atom type.
#[derive(Clone, Debug, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AtomTypeDef {
    /// Unique type name.
    pub name: String,
    /// Atomic mass (amu).
    pub mass: f64,
    /// Lennard-Jones σ (nm).
    pub sigma: f64,
    /// Lennard-Jones ε (kJ/mol).
    pub epsilon: f64,
}

/// Harmonic bond type.
#[derive(Clone, Debug, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct BondTypeDef {
    /// The two type names.
    pub types: [String; 2],
    /// Force constant (kJ mol⁻¹ nm⁻²).
    pub k: f64,
    /// Equilibrium length (nm).
    pub r0: f64,
}

/// Harmonic angle type, with an optional Urey-Bradley term.
#[derive(Clone, Debug, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AngleTypeDef {
    /// The three type names (central type second).
    pub types: [String; 3],
    /// Force constant (kJ mol⁻¹ rad⁻²).
    pub k: f64,
    /// Equilibrium angle (degrees).
    pub theta0: f64,
    /// Urey-Bradley force constant (kJ mol⁻¹ nm⁻²).
    pub urey_bradley_k: Option<f64>,
    /// Urey-Bradley equilibrium 1-3 distance (nm).
    pub urey_bradley_r0: Option<f64>,
}

/// One periodic term of a dihedral type.
#[derive(Clone, Copy, Debug, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DihedralTermDef {
    /// Force constant (kJ/mol).
    pub k: f64,
    /// Multiplicity.
    pub n: u32,
    /// Phase (degrees).
    pub phase: f64,
}

/// Proper dihedral type (pattern of four type names, `"*"` allowed).
#[derive(Clone, Debug, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DihedralTypeDef {
    /// Type pattern.
    pub types: [String; 4],
    /// Periodic terms.
    pub terms: Vec<DihedralTermDef>,
}

/// Periodic improper torsion type (pattern of four type names, `"*"` allowed).
///
/// Each term is `k (1 + cos(n φ − phase))`, the form Amber uses for impropers. An
/// improper matched by one of these types uses its terms; an improper that matches
/// none falls back to [`ImproperTypeDef`].
#[derive(Clone, Debug, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PeriodicImproperTypeDef {
    /// Type pattern.
    pub types: [String; 4],
    /// Periodic terms.
    pub terms: Vec<DihedralTermDef>,
}

/// Improper dihedral type (pattern of four type names, `"*"` allowed).
#[derive(Clone, Debug, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ImproperTypeDef {
    /// Type pattern.
    pub types: [String; 4],
    /// Force constant (kJ mol⁻¹ rad⁻²).
    pub k: f64,
    /// Reference angle (degrees).
    pub xi0: f64,
}

/// One atom of a molecule template.
#[derive(Clone, Debug, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TemplateAtom {
    /// Atom name (for output).
    pub name: String,
    /// Atom type name.
    #[serde(rename = "type")]
    pub atom_type: String,
    /// Partial charge (e).
    pub charge: f64,
    /// Reference position within the molecule (nm).
    pub position: [f64; 3],
}

/// A molecule template.
#[derive(Clone, Debug, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct MoleculeDef {
    /// Template name, used when placing molecules.
    pub name: String,
    /// Atoms in template order.
    pub atoms: Vec<TemplateAtom>,
    /// Bonds as pairs of template atom indices.
    #[serde(default)]
    pub bonds: Vec<[usize; 2]>,
    /// Explicit impropers as quadruples of template atom indices.
    #[serde(default)]
    pub impropers: Vec<[usize; 4]>,
    /// Atom that bonds to the previous residue of a chain (residue templates only).
    #[serde(default)]
    pub head: Option<usize>,
    /// Atom that bonds to the next residue of a chain (residue templates only).
    #[serde(default)]
    pub tail: Option<usize>,
}

/// A parsed, validated force-field parameter set.
#[derive(Clone, Debug, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ParamSet {
    /// Name of the parameter set.
    pub name: String,
    /// Format version (must equal [`FORMAT_VERSION`]).
    pub format_version: u32,
    /// 1-4 Lennard-Jones scale.
    #[serde(default = "default_fudge_lj")]
    pub fudge_lj: f64,
    /// 1-4 Coulomb scale.
    #[serde(default = "default_fudge_qq")]
    pub fudge_qq: f64,
    /// Atom types; a system's `atom_types` index into this list.
    #[serde(default)]
    pub atomtype: Vec<AtomTypeDef>,
    /// Bond types.
    #[serde(default)]
    pub bondtype: Vec<BondTypeDef>,
    /// Angle types.
    #[serde(default)]
    pub angletype: Vec<AngleTypeDef>,
    /// Proper dihedral types.
    #[serde(default)]
    pub dihedraltype: Vec<DihedralTypeDef>,
    /// Improper dihedral types.
    #[serde(default)]
    pub impropertype: Vec<ImproperTypeDef>,
    /// Periodic improper torsion types, tried before `impropertype`.
    #[serde(default)]
    pub periodicimpropertype: Vec<PeriodicImproperTypeDef>,
    /// Molecule templates.
    #[serde(default)]
    pub molecule: Vec<MoleculeDef>,
}

/// Options for building a force field from a parameter set and a system.
#[derive(Clone, Debug, PartialEq)]
pub struct FfOptions {
    /// Unit system (sets the Coulomb factor).
    pub units: UnitSystem,
    /// Non-bonded cutoff (nm).
    pub cutoff: f64,
    /// Lennard-Jones cutoff treatment.
    pub modifier: Modifier,
    /// Electrostatic treatment.
    pub electrostatics: Electrostatics,
}

/// Rejects chain links whose neighbouring angles are nearly straight.
///
/// The harmonic angle force carries a `1/sin θ` factor, so a straight angle makes
/// the forces singular and the energy non-smooth there. A residue template whose
/// bonds at a link line up with the link direction hits this. Angles with
/// `|cos θ| > 0.996` (within about 5° of straight) are refused.
fn check_link_angles(
    positions: &[Vec3],
    bonds: &[(usize, usize)],
    links: &[(usize, usize)],
) -> Result<(), SimError> {
    let mut adjacency: Vec<Vec<usize>> = vec![Vec::new(); positions.len()];
    for &(i, j) in bonds {
        adjacency[i].push(j);
        adjacency[j].push(i);
    }
    for &(a, b) in links {
        for (centre, other) in [(a, b), (b, a)] {
            for &n in &adjacency[centre] {
                if n == other {
                    continue;
                }
                let u = positions[n] - positions[centre];
                let v = positions[other] - positions[centre];
                let cos = u.dot(v) / (u.norm() * v.norm());
                if cos.abs() > 0.996 {
                    return Err(SimError::Config(format!(
                        "link {a}-{b} makes a nearly straight angle at atom {centre}; \
                         give the residue templates a bent geometry at their head and tail"
                    )));
                }
            }
        }
    }
    Ok(())
}

/// Reference position of a template atom, as a vector.
fn template_pos_of(p: &[f64; 3]) -> Vec3 {
    Vec3::new(p[0], p[1], p[2])
}

/// Reference position of template atom `i` (caller checks the index).
fn template_pos(m: &MoleculeDef, i: usize) -> Vec3 {
    template_pos_of(&m.atoms[i].position)
}

/// A system and its bonded topology, produced by [`ParamSet::assemble`].
#[derive(Clone, Debug)]
pub struct Assembled {
    /// The system, with `topology` filled in.
    pub system: System,
}

impl ParamSet {
    /// Parses and validates a parameter file.
    ///
    /// # Errors
    /// Returns [`SimError::Config`] for syntax errors, unknown keys, an
    /// unsupported format version, duplicate or undefined type names, or
    /// non-finite parameters.
    pub fn from_toml_str(text: &str) -> Result<Self, SimError> {
        let p: Self = toml::from_str(text).map_err(|e| SimError::Config(e.to_string()))?;
        p.validate()?;
        Ok(p)
    }

    /// Checks internal consistency.
    ///
    /// # Errors
    /// Returns [`SimError::Config`] describing the first problem found.
    pub fn validate(&self) -> Result<(), SimError> {
        let bad = |m: String| Err(SimError::Config(m));
        if self.format_version != FORMAT_VERSION {
            return bad(format!("unsupported format_version {}", self.format_version));
        }
        if !(self.fudge_lj.is_finite() && self.fudge_qq.is_finite()) {
            return bad("fudge factors must be finite".into());
        }
        for (i, a) in self.atomtype.iter().enumerate() {
            if self.atomtype[..i].iter().any(|b| b.name == a.name) {
                return bad(format!("duplicate atom type {}", a.name));
            }
            if !(a.mass.is_finite() && a.mass > 0.0) {
                return bad(format!("atom type {}: mass must be > 0", a.name));
            }
            if !(a.sigma.is_finite() && a.sigma > 0.0 && a.epsilon.is_finite() && a.epsilon >= 0.0)
            {
                return bad(format!("atom type {}: invalid Lennard-Jones parameters", a.name));
            }
        }
        let known = |t: &str| t == "*" || self.atomtype.iter().any(|a| a.name == t);
        for b in &self.bondtype {
            if !(b.types.iter().all(|t| known(t)) && b.types.iter().all(|t| t != "*")) {
                return bad(format!("bond type {:?} uses an unknown atom type", b.types));
            }
            if !(b.k.is_finite() && b.k >= 0.0 && b.r0.is_finite() && b.r0 > 0.0) {
                return bad(format!("bond type {:?}: invalid parameters", b.types));
            }
        }
        for a in &self.angletype {
            if !a.types.iter().all(|t| known(t) && t != "*") {
                return bad(format!("angle type {:?} uses an unknown atom type", a.types));
            }
            if !(a.k.is_finite()
                && a.k >= 0.0
                && a.theta0.is_finite()
                && a.theta0 > 0.0
                && a.theta0 <= 180.0)
            {
                return bad(format!("angle type {:?}: invalid parameters", a.types));
            }
            match (a.urey_bradley_k, a.urey_bradley_r0) {
                (None, None) => {}
                (Some(k), Some(r0)) if k.is_finite() && k >= 0.0 && r0.is_finite() && r0 > 0.0 => {}
                _ => {
                    return bad(format!(
                        "angle type {:?}: invalid Urey-Bradley parameters",
                        a.types
                    ))
                }
            }
        }
        for d in &self.dihedraltype {
            if !d.types.iter().all(|t| known(t)) {
                return bad(format!("dihedral type {:?} uses an unknown atom type", d.types));
            }
            for t in &d.terms {
                if !(t.k.is_finite() && t.phase.is_finite() && t.n >= 1) {
                    return bad(format!("dihedral type {:?}: invalid term", d.types));
                }
            }
        }
        for d in &self.periodicimpropertype {
            if !d.types.iter().all(|t| known(t)) {
                return bad(format!(
                    "periodic improper type {:?} uses an unknown atom type",
                    d.types
                ));
            }
            for t in &d.terms {
                if !(t.k.is_finite() && t.phase.is_finite() && t.n >= 1) {
                    return bad(format!("periodic improper type {:?}: invalid term", d.types));
                }
            }
        }
        for d in &self.impropertype {
            if !d.types.iter().all(|t| known(t)) {
                return bad(format!("improper type {:?} uses an unknown atom type", d.types));
            }
            if !(d.k.is_finite() && d.k >= 0.0 && d.xi0.is_finite()) {
                return bad(format!("improper type {:?}: invalid parameters", d.types));
            }
        }
        self.validate_molecules()?;
        Ok(())
    }

    fn atom_type_index(&self, name: &str) -> Option<usize> {
        self.atomtype.iter().position(|a| a.name == name)
    }

    fn molecule(&self, name: &str) -> Result<&MoleculeDef, SimError> {
        self.molecule
            .iter()
            .find(|m| m.name == name)
            .ok_or_else(|| SimError::Config(format!("no molecule template named {name}")))
    }

    /// Checks every molecule template: atom types, finite values and index ranges.
    fn validate_molecules(&self) -> Result<(), SimError> {
        let bad = |m: String| Err(SimError::Config(m));
        for m in &self.molecule {
            let n = m.atoms.len();
            if n == 0 {
                return bad(format!("molecule {} has no atoms", m.name));
            }
            for a in &m.atoms {
                if !self.atomtype.iter().any(|t| t.name == a.atom_type) {
                    return bad(format!(
                        "molecule {}: atom {} has unknown type {}",
                        m.name, a.name, a.atom_type
                    ));
                }
                if !(a.charge.is_finite() && a.position.iter().all(|x| x.is_finite())) {
                    return bad(format!(
                        "molecule {}: atom {} has non-finite values",
                        m.name, a.name
                    ));
                }
            }
            for &[i, j] in &m.bonds {
                if i >= n || j >= n || i == j {
                    return bad(format!("molecule {}: bond ({i}, {j}) is invalid", m.name));
                }
            }
            for q in &m.impropers {
                if q.iter().any(|&x| x >= n) {
                    return bad(format!("molecule {}: improper {q:?} is invalid", m.name));
                }
            }
            for (label, idx) in [("head", m.head), ("tail", m.tail)] {
                if idx.is_some_and(|i| i >= n) {
                    return bad(format!("molecule {}: {label} atom index is out of range", m.name));
                }
            }
        }
        Ok(())
    }

    /// Builds a system containing one copy of each placed molecule. Each entry
    /// of `placements` is a template name and a translation (nm). Velocities are
    /// zero and the topology lists are filled in.
    ///
    /// # Errors
    /// Returns [`SimError::Config`] for unknown templates, and
    /// [`SimError::Invalid`] if the generated topology is inconsistent.
    pub fn assemble(
        &self,
        placements: &[(&str, Vec3)],
        sim_box: SimBox,
    ) -> Result<Assembled, SimError> {
        let mut system = System::new(sim_box);
        let mut topo = Topology::default();
        for (resid, (name, offset)) in placements.iter().enumerate() {
            let m = self.molecule(name)?;
            let base = system.len();
            for a in &m.atoms {
                let mass = self.atomtype[self.atom_type_index(&a.atom_type).unwrap_or(0)].mass;
                system.push_atom(NewAtom {
                    name: a.name.clone(),
                    residue_name: m.name.chars().take(3).collect(),
                    residue_id: i32::try_from(resid + 1).unwrap_or(i32::MAX),
                    element: String::new(),
                    position: Vec3::new(a.position[0], a.position[1], a.position[2]) + *offset,
                    velocity: Vec3::ZERO,
                    mass,
                    charge: a.charge,
                    atom_type: self.atom_type_index(&a.atom_type).unwrap_or(0) as u16,
                });
            }
            let bonds: Vec<(usize, usize)> = m.bonds.iter().map(|&[i, j]| (i, j)).collect();
            // Angles, dihedrals and pairs come from the bonds; impropers are explicit.
            let mut mol_topo = Topology::from_bonds(m.atoms.len(), &bonds)?;
            mol_topo.impropers.clone_from(&m.impropers);
            topo.append_shifted(&mol_topo, base);
        }
        system.topology = topo;
        system.validate()?;
        Ok(Assembled { system })
    }

    /// Builds linear chains of residues: proteins, nucleic acids and lipid tails.
    ///
    /// Each chain is a sequence of template names and an origin. Residue `k+1`
    /// bonds from its `head` atom to the `tail` atom of residue `k`. Each residue
    /// after the first is translated so that its head lies `bond_length` nm along
    /// +x from the previous tail, so every template keeps its reference geometry.
    /// The first residue of a chain is translated to the chain's origin. Angles,
    /// dihedrals and 1-4 pairs come from the full bond graph, so they span residue
    /// boundaries. Several chains may be combined in one system (for example a
    /// membrane of lipid chains).
    ///
    /// # Errors
    /// Returns [`SimError::Config`] for unknown templates, a non-positive
    /// `bond_length`, an empty chain, or a residue that cannot be linked (no
    /// `tail` before another residue, or no `head` after one). Returns
    /// [`SimError::Invalid`] if the generated topology is inconsistent.
    pub fn assemble_chains(
        &self,
        chains: &[(&[&str], Vec3)],
        bond_length: f64,
        sim_box: SimBox,
    ) -> Result<Assembled, SimError> {
        if !(bond_length.is_finite() && bond_length > 0.0) {
            return Err(SimError::Config("bond_length must be positive and finite".into()));
        }
        let mut system = System::new(sim_box);
        let mut bonds: Vec<(usize, usize)> = Vec::new();
        let mut links: Vec<(usize, usize)> = Vec::new();
        let mut impropers: Vec<[usize; 4]> = Vec::new();
        let mut resid = 0usize;
        for (residues, origin) in chains {
            if residues.is_empty() {
                return Err(SimError::Config("chain has no residues".into()));
            }
            // Global index and world position of the tail atom of the previous residue.
            let mut link: Option<(usize, Vec3)> = None;
            for (k, name) in residues.iter().enumerate() {
                let m = self.molecule(name)?;
                let base = system.len();
                let shift = match link {
                    None => *origin,
                    Some((_, tail_world)) => {
                        let head = m.head.ok_or_else(|| {
                            SimError::Config(format!(
                                "residue {name} follows another residue but has no head atom"
                            ))
                        })?;
                        let h = template_pos(m, head);
                        tail_world + Vec3::new(bond_length, 0.0, 0.0) - h
                    }
                };
                if let (Some((tail_global, _)), Some(head)) = (link, m.head) {
                    bonds.push((tail_global, base + head));
                    links.push((tail_global, base + head));
                }
                for a in &m.atoms {
                    let mass = self.atomtype[self.atom_type_index(&a.atom_type).unwrap_or(0)].mass;
                    system.push_atom(NewAtom {
                        name: a.name.clone(),
                        residue_name: m.name.chars().take(3).collect(),
                        residue_id: i32::try_from(resid + 1).unwrap_or(i32::MAX),
                        element: String::new(),
                        position: template_pos_of(&a.position) + shift,
                        velocity: Vec3::ZERO,
                        mass,
                        charge: a.charge,
                        atom_type: self.atom_type_index(&a.atom_type).unwrap_or(0) as u16,
                    });
                }
                bonds.extend(m.bonds.iter().map(|&[i, j]| (base + i, base + j)));
                impropers.extend(m.impropers.iter().map(|q| q.map(|x| base + x)));
                resid += 1;
                link = if k + 1 < residues.len() {
                    let t = m.tail.ok_or_else(|| {
                        SimError::Config(format!(
                            "residue {name} is followed by another residue but has no tail atom"
                        ))
                    })?;
                    Some((base + t, template_pos(m, t) + shift))
                } else {
                    None
                };
            }
        }
        check_link_angles(&system.positions, &bonds, &links)?;
        let mut topo = Topology::from_bonds(system.len(), &bonds)?;
        topo.impropers = impropers;
        system.topology = topo;
        system.validate()?;
        Ok(Assembled { system })
    }

    /// Binds parameters to the bonded terms of `system`.
    ///
    /// # Errors
    /// Returns [`SimError::Config`] if a bonded term has no matching parameters.
    pub fn bonded_terms(&self, system: &System) -> Result<BondedTerms, SimError> {
        let type_name = |a: usize| -> &str {
            self.atomtype.get(system.atom_types[a] as usize).map_or("", |t| t.name.as_str())
        };
        let mut terms = BondedTerms::default();
        for &(i, j) in &system.topology.bonds {
            let def = self
                .bondtype
                .iter()
                .find(|b| same_pair(&b.types, type_name(i), type_name(j)))
                .ok_or_else(|| missing("bond", &[type_name(i), type_name(j)]))?;
            terms.bonds.push(HarmonicBond { i, j, k: def.k, r0: def.r0 });
        }
        for &[i, j, k] in &system.topology.angles {
            let names = [type_name(i), type_name(j), type_name(k)];
            let def = self
                .angletype
                .iter()
                .find(|a| same_triple(&a.types, names))
                .ok_or_else(|| missing("angle", &names))?;
            let urey_bradley = match (def.urey_bradley_k, def.urey_bradley_r0) {
                (Some(k), Some(r0)) => Some(UreyBradley { k, r0 }),
                _ => None,
            };
            terms.angles.push(HarmonicAngle {
                i,
                j,
                k,
                k_theta: def.k,
                theta0: def.theta0.to_radians(),
                urey_bradley,
            });
        }
        for &atoms in &system.topology.dihedrals {
            let names = atoms.map(type_name);
            let def = self
                .dihedraltype
                .iter()
                .find(|d| pattern_matches(&d.types, names))
                .ok_or_else(|| missing("dihedral", &names))?;
            let terms_out = def
                .terms
                .iter()
                .map(|t| DihedralTerm { k: t.k, n: t.n, phase: t.phase.to_radians() })
                .collect();
            terms.dihedrals.push(PeriodicDihedral { atoms, terms: terms_out });
        }
        for &atoms in &system.topology.impropers {
            let names = atoms.map(type_name);
            if let Some(def) =
                self.periodicimpropertype.iter().find(|d| pattern_matches(&d.types, names))
            {
                let terms_out = def
                    .terms
                    .iter()
                    .map(|t| DihedralTerm { k: t.k, n: t.n, phase: t.phase.to_radians() })
                    .collect();
                terms.dihedrals.push(PeriodicDihedral { atoms, terms: terms_out });
                continue;
            }
            let def = self
                .impropertype
                .iter()
                .find(|d| pattern_matches(&d.types, names))
                .ok_or_else(|| missing("improper", &names))?;
            terms.impropers.push(HarmonicImproper { atoms, k: def.k, xi0: def.xi0.to_radians() });
        }
        Ok(terms)
    }

    /// Non-bonded settings for this parameter set.
    pub fn nonbonded_settings(&self, opts: &FfOptions) -> NonbondedSettings {
        NonbondedSettings {
            lj_types: self
                .atomtype
                .iter()
                .map(|a| LjType { epsilon: a.epsilon, sigma: a.sigma })
                .collect(),
            cutoff: opts.cutoff,
            modifier: opts.modifier,
            electrostatics: opts.electrostatics,
            coulomb_factor: match opts.units {
                UnitSystem::Md => COULOMB_FACTOR_MD,
                UnitSystem::Reduced => 1.0,
            },
            fudge_lj: self.fudge_lj,
            fudge_qq: self.fudge_qq,
        }
    }

    /// The complete force field (bonded plus non-bonded) for `system`, whose
    /// topology must have been produced by [`ParamSet::assemble`].
    ///
    /// # Errors
    /// Returns [`SimError::Config`] for missing parameters or invalid options.
    pub fn force_field(
        &self,
        system: &System,
        opts: &FfOptions,
    ) -> Result<CombinedForceField, SimError> {
        let n = system.len();
        let terms = self.bonded_terms(system)?;
        let bonded = BondedForce::new(terms, n)?;
        let excl = system.topology.exclusion_partners(n);
        let nb =
            Nonbonded::new(self.nonbonded_settings(opts), excl, system.topology.pairs14.clone())?;
        Ok(CombinedForceField::new(vec![Box::new(nb) as Box<dyn ForceField>, Box::new(bonded)]))
    }
}

fn missing(kind: &str, names: &[&str]) -> SimError {
    SimError::Config(format!("no {kind} parameters for atom types {names:?}"))
}

fn same_pair(types: &[String; 2], a: &str, b: &str) -> bool {
    (types[0] == a && types[1] == b) || (types[0] == b && types[1] == a)
}

fn same_triple(types: &[String; 3], n: [&str; 3]) -> bool {
    (types[0] == n[0] && types[1] == n[1] && types[2] == n[2])
        || (types[0] == n[2] && types[1] == n[1] && types[2] == n[0])
}

/// Matches a type pattern (with `*` wildcards) against a quadruple, forward or reversed.
fn pattern_matches(pat: &[String; 4], n: [&str; 4]) -> bool {
    let fwd = pat.iter().zip(n).all(|(p, t)| p == "*" || p == t);
    let rev = pat.iter().zip(n.iter().rev()).all(|(p, t)| p == "*" || p == t);
    fwd || rev
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::integrator::{ForceEvaluator, VelocityVerlet};
    use crate::neighbour::VerletList;
    use crate::simulation::Simulation;

    /// The synthetic test parameter set shipped with the crate.
    pub(crate) const SYNTHETIC: &str = include_str!("../data/synthetic-test.toml");

    fn params() -> ParamSet {
        ParamSet::from_toml_str(SYNTHETIC).unwrap()
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
    fn synthetic_file_parses_and_is_labelled() {
        let p = params();
        assert!(p.name.to_lowercase().contains("synthetic"));
        assert!(SYNTHETIC.contains("SYNTHETIC"));
        assert!(p.molecule.iter().any(|m| m.name == "butane"));
    }

    #[test]
    fn butane_topology_counts() {
        let p = params();
        let asm = p
            .assemble(&[("butane", Vec3::new(0.5, 0.5, 0.5))], SimBox::new(4.0, 4.0, 4.0).unwrap())
            .unwrap();
        let t = &asm.system.topology;
        // Linear C4 chain: 3 bonds, 2 angles, 1 dihedral, 1 pair at distance 3.
        assert_eq!(t.bonds.len(), 3);
        assert_eq!(t.angles.len(), 2);
        assert_eq!(t.dihedrals.len(), 1);
        assert_eq!(t.pairs14.len(), 1);
        assert_eq!(t.exclusions.len(), 3 + 2 + 1);
        assert_eq!(t.impropers, Vec::<[usize; 4]>::new());
    }

    #[test]
    fn two_molecules_have_independent_indices() {
        let p = params();
        let asm = p
            .assemble(
                &[("butane", Vec3::new(0.5, 0.5, 0.5)), ("butane", Vec3::new(2.5, 0.5, 0.5))],
                SimBox::new(4.0, 4.0, 4.0).unwrap(),
            )
            .unwrap();
        assert_eq!(asm.system.len(), 8);
        assert_eq!(asm.system.topology.bonds.len(), 6);
        assert!(asm.system.topology.bonds.contains(&(4, 5)));
        // No cross-molecule exclusions.
        assert!(!asm.system.topology.exclusions.iter().any(|&(i, j)| i < 4 && j >= 4));
    }

    #[test]
    fn bonded_parameters_bind_by_type() {
        let p = params();
        let asm =
            p.assemble(&[("butane", Vec3::ZERO)], SimBox::new(4.0, 4.0, 4.0).unwrap()).unwrap();
        let terms = p.bonded_terms(&asm.system).unwrap();
        assert_eq!(terms.bonds.len(), 3);
        // 112.7 degrees in the synthetic file.
        assert!(terms.angles.iter().all(|a| (a.theta0 - 112.7f64.to_radians()).abs() < 1e-12));
        assert_eq!(terms.dihedrals.len(), 1);
        assert_eq!(terms.dihedrals[0].terms.len(), 1);
    }

    #[test]
    fn missing_parameters_are_reported() {
        let p = params();
        let mut asm =
            p.assemble(&[("butane", Vec3::ZERO)], SimBox::new(4.0, 4.0, 4.0).unwrap()).unwrap();
        // Give atom 1 a type that has no bonded parameters.
        let mut q = p.clone();
        q.atomtype.push(AtomTypeDef { name: "XX".into(), mass: 1.0, sigma: 0.3, epsilon: 0.1 });
        asm.system.atom_types[1] = (q.atomtype.len() - 1) as u16;
        assert!(q.bonded_terms(&asm.system).is_err());
    }

    #[test]
    fn rejects_bad_files() {
        assert!(ParamSet::from_toml_str("name = \"x\"\nformat_version = 2\n").is_err());
        assert!(ParamSet::from_toml_str("name = \"x\"\nformat_version = 1\nunknown = 1\n").is_err());
        let dup = format!(
            "{SYNTHETIC}\n[[atomtype]]\nname = \"CT\"\nmass = 1.0\nsigma = 0.3\nepsilon = 0.1\n"
        );
        assert!(ParamSet::from_toml_str(&dup).is_err());
        let bad_angle = format!(
            "{SYNTHETIC}\n[[angletype]]\ntypes = [\"CT\",\"CT\",\"CT\"]\nk = 1.0\ntheta0 = 200.0\n"
        );
        assert!(ParamSet::from_toml_str(&bad_angle).is_err());
    }

    /// Analytic energy of a butane whose bonds, angles and dihedral are at
    /// their equilibrium values is only the non-bonded part; the bonded part
    /// must vanish when the geometry matches the reference template.
    #[test]
    fn reference_geometry_has_small_bonded_energy() {
        let p = params();
        let asm = p
            .assemble(&[("butane", Vec3::new(0.3, 0.3, 0.3))], SimBox::new(4.0, 4.0, 4.0).unwrap())
            .unwrap();
        let bonded =
            BondedForce::new(p.bonded_terms(&asm.system).unwrap(), asm.system.len()).unwrap();
        let list = VerletList::new(1.0, 0.1).unwrap();
        let mut f = vec![Vec3::ZERO; asm.system.len()];
        let e = bonded
            .compute(
                &asm.system.positions,
                &asm.system.atom_types,
                &asm.system.charges,
                &asm.system.sim_box,
                &list,
                &mut f,
            )
            .potential;
        // The template is only approximately at equilibrium, so require a small value.
        assert!(e < 5.0, "bonded energy {e}");
    }

    /// Energy conservation under NVE for a small bonded system, with the full
    /// force field (bonded, exclusions, 1-4 pairs and LJ).
    #[test]
    fn nve_conserves_energy_for_a_bonded_chain() {
        let p = params();
        let mut asm = p
            .assemble(&[("butane", Vec3::new(0.5, 0.5, 0.5))], SimBox::new(4.0, 4.0, 4.0).unwrap())
            .unwrap();
        crate::simulation::maxwell_boltzmann_velocities(
            &mut asm.system,
            0.3,
            4,
            UnitSystem::Reduced,
        );
        let ff = p.force_field(&asm.system, &opts()).unwrap();
        let eval = ForceEvaluator::new(Box::new(ff), 0.2).unwrap();
        let mut sim = Simulation::new(
            asm.system,
            eval,
            Box::new(VelocityVerlet),
            None,
            0.0005,
            UnitSystem::Reduced,
        )
        .unwrap();
        let summary = sim.run(4000, &mut []).unwrap();
        assert!(
            summary.relative_energy_drift().abs() < 2e-3,
            "drift {}",
            summary.relative_energy_drift()
        );
    }

    /// Random text must produce errors, never panics.
    #[test]
    fn parser_never_panics_on_random_input() {
        use crate::rng::Rng;
        let mut rng = Rng::seed_from_u64(31337);
        let alphabet = b"[]{}=\"'\n ,.#-+*0123456789abcxyztypemoleculeatomtype";
        for _ in 0..3000 {
            let len = (rng.next_u64() % 300) as usize;
            let s: String = (0..len)
                .map(|_| alphabet[(rng.next_u64() % alphabet.len() as u64) as usize] as char)
                .collect();
            let _ = ParamSet::from_toml_str(&s);
        }
        // Every prefix of the valid file.
        for cut in (0..SYNTHETIC.len()).step_by(7) {
            let _ = ParamSet::from_toml_str(&SYNTHETIC[..cut]);
        }
    }
}
