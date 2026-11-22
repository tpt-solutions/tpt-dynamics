//! The molecular system: per-atom state stored in structure-of-arrays (`SoA`)
//! layout, the periodic box and the bonded topology.
//!
//! Units follow the GROMACS convention: nm, ps, amu, elementary charge, kJ/mol
//! and Kelvin. With these units `½ m v²` is directly in kJ/mol. Reduced
//! Lennard-Jones units (σ = ε = m = 1) are equally valid; the engine does not
//! assume a particular unit system beyond this bookkeeping.

use tpt_dynamics_core::{constants::GAS_CONSTANT, Vec3};

use crate::error::SimError;
use crate::pbc::SimBox;
pub use crate::topology::Topology;

/// Removes the entries of `v` whose flag in `keep` is false.
fn compact<T>(v: &mut Vec<T>, keep: &[bool]) {
    let mut flags = keep.iter();
    v.retain(|_| *flags.next().unwrap_or(&false));
}

/// Boltzmann constant in kJ mol⁻¹ K⁻¹.
pub const BOLTZMANN_KJ: f64 = GAS_CONSTANT / 1000.0;

/// Unit system used for temperature conversion.
///
/// `Md` uses kJ mol⁻¹ K⁻¹ (GROMACS-style units); `Reduced` sets `k_B` = 1 for
/// Lennard-Jones reduced units (σ = ε = m = 1).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "kebab-case")]
pub enum UnitSystem {
    /// nm, ps, amu, kJ/mol, K.
    #[default]
    Md,
    /// Reduced units with `k_B` = 1.
    Reduced,
}

impl UnitSystem {
    /// Boltzmann constant in this unit system.
    pub fn boltzmann(self) -> f64 {
        match self {
            Self::Md => BOLTZMANN_KJ,
            Self::Reduced => 1.0,
        }
    }
}

/// Descriptive data for one atom, supplied when adding it to a [`System`].
#[derive(Clone, Debug, PartialEq)]
pub struct NewAtom {
    /// Atom name, e.g. `"CA"` or `"OW"` (at most 4 characters for PDB output).
    pub name: String,
    /// Residue name, e.g. `"ALA"` or `"SOL"` (at most 3 characters for PDB output).
    pub residue_name: String,
    /// Residue number.
    pub residue_id: i32,
    /// Element symbol, e.g. `"C"` or `"Cl"`; may be empty for coarse-grained beads.
    pub element: String,
    /// Position in nm.
    pub position: Vec3,
    /// Velocity in nm/ps.
    pub velocity: Vec3,
    /// Mass in amu (must be finite and > 0 before dynamics).
    pub mass: f64,
    /// Partial charge in elementary charge units.
    pub charge: f64,
    /// Index into the force field's per-type parameter table.
    pub atom_type: u16,
}

/// A molecular system in structure-of-arrays layout.
///
/// Every per-atom vector has length [`System::len`]. Positions, velocities and
/// forces are kept as separate `Vec<Vec3>` so each array can be streamed
/// independently by the integrators and force kernels.
#[derive(Clone, Debug, PartialEq)]
pub struct System {
    /// Positions in nm.
    pub positions: Vec<Vec3>,
    /// Velocities in nm/ps.
    pub velocities: Vec<Vec3>,
    /// Forces in kJ mol⁻¹ nm⁻¹ (valid after a force evaluation).
    pub forces: Vec<Vec3>,
    /// Masses in amu.
    pub masses: Vec<f64>,
    /// Partial charges in elementary charge units.
    pub charges: Vec<f64>,
    /// Atom type indices.
    pub atom_types: Vec<u16>,
    /// Atom names.
    pub names: Vec<String>,
    /// Residue names.
    pub residue_names: Vec<String>,
    /// Residue numbers.
    pub residue_ids: Vec<i32>,
    /// Element symbols.
    pub elements: Vec<String>,
    /// Periodic box.
    pub sim_box: SimBox,
    /// Bonded and pair topology (empty for PDB/GRO-only systems).
    pub topology: Topology,
    /// Number of holonomic constraints (bond or SETTLE distances) removed from
    /// the degrees of freedom. Set by the constraint algorithms.
    pub n_constraints: usize,
}

impl System {
    /// Creates an empty system in the given periodic box.
    pub fn new(sim_box: SimBox) -> Self {
        Self {
            positions: Vec::new(),
            velocities: Vec::new(),
            forces: Vec::new(),
            masses: Vec::new(),
            charges: Vec::new(),
            atom_types: Vec::new(),
            names: Vec::new(),
            residue_names: Vec::new(),
            residue_ids: Vec::new(),
            elements: Vec::new(),
            sim_box,
            topology: Topology::default(),
            n_constraints: 0,
        }
    }

    /// Number of atoms.
    pub fn len(&self) -> usize {
        self.positions.len()
    }

    /// True if the system has no atoms.
    pub fn is_empty(&self) -> bool {
        self.positions.is_empty()
    }

    /// Appends an atom. The position is wrapped into the primary box.
    pub fn push_atom(&mut self, atom: NewAtom) {
        self.positions.push(self.sim_box.wrap(atom.position));
        self.velocities.push(atom.velocity);
        self.forces.push(Vec3::ZERO);
        self.masses.push(atom.mass);
        self.charges.push(atom.charge);
        self.atom_types.push(atom.atom_type);
        self.names.push(atom.name);
        self.residue_names.push(atom.residue_name);
        self.residue_ids.push(atom.residue_id);
        self.elements.push(atom.element);
    }

    /// Checks that all per-atom arrays have equal length and that masses are
    /// positive and finite, as required before running dynamics.
    pub fn validate(&self) -> Result<(), SimError> {
        let n = self.len();
        let lengths = [
            self.velocities.len(),
            self.forces.len(),
            self.masses.len(),
            self.charges.len(),
            self.atom_types.len(),
            self.names.len(),
            self.residue_names.len(),
            self.residue_ids.len(),
            self.elements.len(),
        ];
        if lengths.iter().any(|&l| l != n) {
            return Err(SimError::Invalid("per-atom arrays have inconsistent lengths".into()));
        }
        if let Some(i) = self.masses.iter().position(|&m| !(m.is_finite() && m > 0.0)) {
            return Err(SimError::Invalid(format!("atom {i} has a non-positive mass")));
        }
        if self.positions.iter().any(|p| !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite())) {
            return Err(SimError::Invalid("non-finite coordinate".into()));
        }
        Ok(())
    }

    /// Degrees of freedom after removing the three centre-of-mass translations.
    ///
    /// Holonomic constraints recorded in [`System::n_constraints`] are also
    /// subtracted.
    pub fn degrees_of_freedom(&self) -> f64 {
        let n = self.len() as f64;
        let base = if n > 1.0 { 3.0 * n - 3.0 } else { 3.0 * n };
        (base - self.n_constraints as f64).max(0.0)
    }

    /// Keeps only the atoms whose flag in `keep` is true, compacting every
    /// per-atom array. The topology and constraint count are cleared because
    /// their indices would no longer be valid.
    ///
    /// # Panics
    /// Panics if `keep.len() != self.len()`.
    pub fn retain(&mut self, keep: &[bool]) {
        assert_eq!(keep.len(), self.len(), "keep mask length must equal the atom count");
        compact(&mut self.positions, keep);
        compact(&mut self.velocities, keep);
        compact(&mut self.forces, keep);
        compact(&mut self.masses, keep);
        compact(&mut self.charges, keep);
        compact(&mut self.atom_types, keep);
        compact(&mut self.names, keep);
        compact(&mut self.residue_names, keep);
        compact(&mut self.residue_ids, keep);
        compact(&mut self.elements, keep);
        self.topology = Topology::default();
        self.n_constraints = 0;
    }

    /// Total kinetic energy `Σ ½ m v²` in kJ/mol.
    pub fn kinetic_energy(&self) -> f64 {
        self.masses.iter().zip(&self.velocities).map(|(&m, v)| 0.5 * m * v.norm_sq()).sum()
    }

    /// Instantaneous temperature in K from the kinetic energy and the
    /// centre-of-mass-corrected degree-of-freedom count.
    pub fn temperature(&self) -> f64 {
        self.temperature_in(UnitSystem::Md)
    }

    /// Instantaneous temperature in the given unit system.
    pub fn temperature_in(&self, units: UnitSystem) -> f64 {
        let dof = self.degrees_of_freedom();
        if dof <= 0.0 {
            return 0.0;
        }
        2.0 * self.kinetic_energy() / (dof * units.boltzmann())
    }

    /// Total linear momentum `Σ m v` in amu nm/ps.
    pub fn total_momentum(&self) -> Vec3 {
        self.masses.iter().zip(&self.velocities).fold(Vec3::ZERO, |acc, (&m, &v)| acc + v * m)
    }

    /// Subtracts the centre-of-mass velocity so the total momentum is zero.
    pub fn remove_com_velocity(&mut self) {
        let total_mass: f64 = self.masses.iter().sum();
        if total_mass <= 0.0 {
            return;
        }
        let vcom = self.total_momentum() * (1.0 / total_mass);
        for v in &mut self.velocities {
            *v -= vcom;
        }
    }
}

/// Standard atomic masses (amu) for common elements. Unknown symbols return `None`.
pub fn element_mass(symbol: &str) -> Option<f64> {
    Some(match symbol {
        "H" => 1.008,
        "He" => 4.0026,
        "B" => 10.81,
        "C" => 12.011,
        "N" => 14.007,
        "O" => 15.999,
        "F" => 18.998,
        "Ne" => 20.180,
        "Na" => 22.990,
        "Mg" => 24.305,
        "Si" => 28.086,
        "P" => 30.974,
        "S" => 32.06,
        "Cl" => 35.45,
        "Ar" => 39.948,
        "K" => 39.098,
        "Ca" => 40.078,
        "Fe" => 55.845,
        "Cu" => 63.546,
        "Zn" => 65.38,
        "Br" => 79.904,
        "I" => 126.904,
        _ => return None,
    })
}

/// Normalises an element symbol to the `Xx` capitalisation used in tables.
pub fn normalise_element(symbol: &str) -> String {
    let mut chars = symbol.chars();
    match chars.next() {
        Some(first) => {
            let mut s: String = first.to_uppercase().collect();
            s.extend(chars.map(|c| c.to_ascii_lowercase()));
            s
        }
        None => String::new(),
    }
}

/// Guesses an element symbol from an atom name and residue name.
///
/// Heuristic: if the residue name equals the atom name and the name is a known
/// two-letter element (e.g. ion residues `NA`/`CL`), that element is used;
/// otherwise the first alphabetic character of the atom name is used.
pub fn guess_element(atom_name: &str, residue_name: &str) -> Option<String> {
    let name = atom_name.trim();
    let res = residue_name.trim();
    if name.eq_ignore_ascii_case(res) && name.len() == 2 {
        let cand = normalise_element(name);
        if element_mass(&cand).is_some() {
            return Some(cand);
        }
    }
    name.chars().find(char::is_ascii_alphabetic).map(|c| c.to_ascii_uppercase().to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atom(name: &str, m: f64, v: Vec3) -> NewAtom {
        NewAtom {
            name: name.into(),
            residue_name: "RES".into(),
            residue_id: 1,
            element: "C".into(),
            position: Vec3::ZERO,
            velocity: v,
            mass: m,
            charge: 0.0,
            atom_type: 0,
        }
    }

    #[test]
    fn momentum_removed_and_temperature_positive() {
        let mut s = System::new(SimBox::new(5.0, 5.0, 5.0).unwrap());
        s.push_atom(atom("A", 1.0, Vec3::new(1.0, 2.0, 3.0)));
        s.push_atom(atom("B", 3.0, Vec3::new(-0.5, 0.0, 1.0)));
        s.push_atom(atom("C", 2.0, Vec3::new(0.0, -1.0, 0.0)));
        s.remove_com_velocity();
        let p = s.total_momentum();
        assert!(p.norm() < 1e-12);
        assert!(s.temperature() > 0.0);
        assert!(s.validate().is_ok());
    }

    #[test]
    fn validate_catches_bad_mass_and_lengths() {
        let mut s = System::new(SimBox::new(5.0, 5.0, 5.0).unwrap());
        s.push_atom(atom("A", 0.0, Vec3::ZERO));
        assert!(s.validate().is_err());
        let mut t = System::new(SimBox::new(5.0, 5.0, 5.0).unwrap());
        t.push_atom(atom("A", 1.0, Vec3::ZERO));
        t.charges.clear();
        assert!(t.validate().is_err());
    }

    #[test]
    fn element_guessing() {
        assert_eq!(guess_element("CA", "ALA").as_deref(), Some("C"));
        assert_eq!(guess_element("NA", "NA").as_deref(), Some("Na"));
        assert_eq!(guess_element("CL", "CL").as_deref(), Some("Cl"));
        assert_eq!(guess_element("HW1", "SOL").as_deref(), Some("H"));
        assert_eq!(guess_element("1HB", "ALA").as_deref(), Some("H"));
        assert_eq!(guess_element("123", "X"), None);
        assert_eq!(element_mass("Cl"), Some(35.45));
    }
}
