//! Rigid water model presets: TIP3P, SPC/E and TIP4P geometry, charges and
//! Lennard-Jones parameters, with SETTLE constraints for the rigid geometry.
//!
//! Units follow the crate convention: nm, kJ/mol and elementary charge. Source
//! values in kcal/mol or K are converted with 4.184 kJ/mol per kcal/mol and
//! `k_B` = 0.0083145 kJ mol⁻¹ K⁻¹.
//!
//! Verification status of each constant (checked against web sources on the
//! date of writing, see [`WaterModel::params`] for the per-field notes):
//!
//! * TIP3P: Jorgensen, Chandrasekhar, Madura, Impey & Klein, J. Chem. Phys.
//!   79, 926 (1983). `r_OH` 0.9572 Å, θ 104.52°, `q_O` −0.834, `q_H` 0.417,
//!   `σ_O` 3.15061 Å, `ε_O` 0.6364 kJ/mol (0.1521 kcal/mol). Sources: the LAMMPS
//!   "Howto tip3p" page (`docs.lammps.org/Howto_tip3p.html`) and the
//!   pysisyphus TIP3P module.
//! * SPC/E: Berendsen, Grigera & Straatsma, J. Phys. Chem. 91, 6269 (1987).
//!   `r_OH` 1.0 Å, θ 109.47°, `q_O` −0.8476, `q_H` 0.4238, `σ_O` 3.166 Å,
//!   `ε_O` 0.65 kJ/mol. Sources: the water-model table at
//!   `pmc.ncbi.nlm.nih.gov/articles/PMC3929574/table/tab1` and the
//!   pysisyphus/LAMMPS summaries; the ε value is rounded in those sources.
//! * TIP4P: Jorgensen et al. 1983 (as above). `r_OH` 0.9572 Å, θ 104.52°,
//!   `q_H` 0.520, M-site 0.15 Å from O on the bisector, ε/`k_B` = 78.0 K
//!   (0.6485 kJ/mol), `σ_O` 3.154 Å. The oxygen charge −1.04 e is **unverified**
//!   (inferred from charge neutrality and recalled from the original paper; the
//!   sources checked did not list it).
//!
//! The M-site of TIP4P is a virtual particle with no mass. Its position is
//! `M = O + a (H1 - O) + a (H2 - O)` with `a` from [`WaterModel::m_site_weight`],
//! which puts it on the H-O-H bisector at the model's M-site distance. Its force
//! is spread onto O (`1 - 2a`) and each H (`a`) by [`spread_m_force`]; this
//! preserves both the total force and the torque. The helper [`Tip4pSites`]
//! applies these operations to a whole system.
//!
//! The simulation loop does not yet place M-sites, evaluate their charge or
//! spread their force, so [`WaterModel::Tip4p`] still reports
//! `is_supported() == false` and the builders in [`crate::config`] refuse it.

use tpt_dynamics_core::Vec3;

use crate::constraints::Settle;
use crate::error::SimError;

/// Supported rigid water models.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WaterModel {
    /// TIP3P (Jorgensen et al. 1983).
    Tip3p,
    /// SPC/E (Berendsen et al. 1987).
    SpcE,
    /// TIP4P (Jorgensen et al. 1983), with a virtual M-site.
    Tip4p,
}

/// Geometry, charges and Lennard-Jones parameters of one water model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterParams {
    /// O–H bond length in nm.
    pub r_oh: f64,
    /// H–O–H angle in degrees.
    pub theta_deg: f64,
    /// Distance of the virtual M-site from O in nm (TIP4P only).
    pub m_site: Option<f64>,
    /// Oxygen partial charge in e.
    pub q_o: f64,
    /// Hydrogen partial charge in e.
    pub q_h: f64,
    /// M-site charge in e (TIP4P only).
    pub q_m: Option<f64>,
    /// Oxygen Lennard-Jones σ in nm. Hydrogens carry no Lennard-Jones term.
    pub sigma_o: f64,
    /// Oxygen Lennard-Jones ε in kJ/mol.
    pub epsilon_o: f64,
}

impl WaterModel {
    /// Parameters of this model. Each field is documented in the module docs;
    /// the TIP4P oxygen charge is the one unverified constant.
    pub fn params(self) -> WaterParams {
        match self {
            Self::Tip3p => WaterParams {
                r_oh: 0.09572,
                theta_deg: 104.52,
                m_site: None,
                q_o: -0.834,
                q_h: 0.417,
                q_m: None,
                sigma_o: 0.315_061,
                epsilon_o: 0.6364,
            },
            Self::SpcE => WaterParams {
                r_oh: 0.1,
                theta_deg: 109.47,
                m_site: None,
                q_o: -0.8476,
                q_h: 0.4238,
                q_m: None,
                sigma_o: 0.3166,
                epsilon_o: 0.65,
            },
            Self::Tip4p => WaterParams {
                r_oh: 0.09572,
                theta_deg: 104.52,
                m_site: Some(0.015),
                // UNVERIFIED: see module docs.
                q_o: 0.0,
                q_h: 0.520,
                q_m: Some(-1.04),
                sigma_o: 0.3154,
                epsilon_o: 0.6485,
            },
        }
    }

    /// True if the model can be simulated by this crate.
    pub fn is_supported(self) -> bool {
        !matches!(self, Self::Tip4p)
    }

    /// H–H distance in nm, `2 r_OH sin(θ/2)`.
    pub fn d_hh(self) -> f64 {
        let p = self.params();
        2.0 * p.r_oh * (0.5 * p.theta_deg.to_radians()).sin()
    }

    /// Weight `a` of each hydrogen in the M-site position
    /// `M = O + a (H1 - O) + a (H2 - O)`, chosen so that `|M - O|` equals the
    /// model's M-site distance. `None` for models without an M-site.
    pub fn m_site_weight(self) -> Option<f64> {
        let p = self.params();
        p.m_site.map(|d| d / (2.0 * p.r_oh * (0.5 * p.theta_deg.to_radians()).cos()))
    }

    /// A SETTLE solver for the given molecules `[O, H1, H2]` in a system of
    /// `n_atoms` atoms.
    ///
    /// # Errors
    /// Returns [`SimError::Config`] for the unsupported TIP4P model, and
    /// propagates [`Settle::new`] errors for invalid indices or tolerance.
    pub fn settle(
        self,
        molecules: Vec<[usize; 3]>,
        n_atoms: usize,
        tolerance: f64,
    ) -> Result<Settle, SimError> {
        if !self.is_supported() {
            return Err(SimError::Config(
                "TIP4P needs virtual-site force redistribution, which is not implemented".into(),
            ));
        }
        let p = self.params();
        Settle::new(molecules, n_atoms, p.r_oh, self.d_hh(), tolerance)
    }
}

/// Position of a TIP4P M-site from its water oxygen and hydrogens, with
/// hydrogen weight `a` from [`WaterModel::m_site_weight`].
pub fn m_site_position(o: Vec3, h1: Vec3, h2: Vec3, a: f64) -> Vec3 {
    o + a * (h1 - o) + a * (h2 - o)
}

/// Spreads a force `f_m` acting on the M-site onto `[O, H1, H2]`: O receives
/// `1 - 2a` of it and each hydrogen `a`. The three forces sum to `f_m`, and
/// their torque about any origin equals the torque of `f_m` at the M-site,
/// because the M-site is an affine combination of the three atoms.
pub fn spread_m_force(f_m: Vec3, a: f64) -> [Vec3; 3] {
    [(1.0 - 2.0 * a) * f_m, a * f_m, a * f_m]
}

/// The TIP4P M-sites of a set of water molecules, each listed as
/// `[oxygen, hydrogen1, hydrogen2]`.
#[derive(Clone, Debug)]
pub struct Tip4pSites {
    molecules: Vec<[usize; 3]>,
    a: f64,
}

impl Tip4pSites {
    /// Creates the M-site helper for `model`, which must have an M-site.
    ///
    /// # Errors
    /// Returns [`SimError::Config`] for a model without an M-site and
    /// [`SimError::Invalid`] for an atom index outside `0..n_atoms`.
    pub fn new(
        model: WaterModel,
        molecules: Vec<[usize; 3]>,
        n_atoms: usize,
    ) -> Result<Self, SimError> {
        let a = model
            .m_site_weight()
            .ok_or_else(|| SimError::Config(format!("{model:?} has no M-site")))?;
        if molecules.iter().flatten().any(|&i| i >= n_atoms) {
            return Err(SimError::Invalid("TIP4P molecule refers to a missing atom".into()));
        }
        Ok(Self { molecules, a })
    }

    /// Hydrogen weight `a` of the M-site construction.
    pub fn weight(&self) -> f64 {
        self.a
    }

    /// M-site position of each molecule, in the order of the molecule list.
    pub fn positions(&self, positions: &[Vec3]) -> Vec<Vec3> {
        self.molecules
            .iter()
            .map(|m| m_site_position(positions[m[0]], positions[m[1]], positions[m[2]], self.a))
            .collect()
    }

    /// Adds the spread of each M-site force onto the oxygen and hydrogens in
    /// `forces`. `m_forces` holds one force per molecule. The M-site forces
    /// are not stored in `forces`, so the caller must handle them separately.
    pub fn spread_forces(&self, m_forces: &[Vec3], forces: &mut [Vec3]) {
        for (m, f_m) in self.molecules.iter().zip(m_forces) {
            for (&atom, f) in m.iter().zip(spread_m_force(*f_m, self.a)) {
                forces[atom] += f;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hh_distance_matches_geometry() {
        // TIP3P: 2 · 0.9572 · sin(52.26°) ≈ 1.5139 Å.
        assert!((WaterModel::Tip3p.d_hh() - 0.15139).abs() < 1e-4);
        // SPC/E: 2 · 1.0 · sin(54.735°) ≈ 1.6330 Å.
        assert!((WaterModel::SpcE.d_hh() - 0.16330).abs() < 1e-4);
    }

    #[test]
    fn tip3p_and_spce_charges_are_neutral() {
        for m in [WaterModel::Tip3p, WaterModel::SpcE] {
            let p = m.params();
            assert!((p.q_o + 2.0 * p.q_h).abs() < 1e-12, "{m:?}");
        }
    }

    #[test]
    fn tip4p_is_refused_by_the_settle_builder() {
        assert!(!WaterModel::Tip4p.is_supported());
        assert!(WaterModel::Tip4p.settle(vec![[0, 1, 2]], 3, 1e-10).is_err());
        assert!(WaterModel::SpcE.settle(vec![[0, 1, 2]], 3, 1e-10).is_ok());
    }

    /// Deterministic pseudo-random numbers in [-1, 1) for the tests.
    fn lcg(state: &mut u64) -> f64 {
        *state =
            state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        ((*state >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
    }

    fn random_vec(state: &mut u64) -> Vec3 {
        Vec3::new(lcg(state), lcg(state), lcg(state))
    }

    /// A TIP4P molecule with O at `origin`, rotated about z and then about x so
    /// that the bisector is not along a coordinate axis.
    fn tip4p_molecule(origin: Vec3) -> [Vec3; 3] {
        let p = WaterModel::Tip4p.params();
        let half = 0.5 * p.theta_deg.to_radians();
        let raw = [
            Vec3::ZERO,
            Vec3::new(p.r_oh * half.sin(), p.r_oh * half.cos(), 0.0),
            Vec3::new(-p.r_oh * half.sin(), p.r_oh * half.cos(), 0.0),
        ];
        let (c, s) = (0.6f64.cos(), 0.6f64.sin());
        let (c2, s2) = (0.3f64.cos(), 0.3f64.sin());
        raw.map(|v| {
            let x = c * v.x - s * v.y;
            let y = s * v.x + c * v.y;
            let z = v.z;
            Vec3::new(x, c2 * y - s2 * z, s2 * y + c2 * z) + origin
        })
    }

    #[test]
    fn tip4p_weight_gives_m_site_distance() {
        let a = WaterModel::Tip4p.m_site_weight().unwrap();
        // 0.015 / (2 · 0.09572 · cos 52.26°) ≈ 0.1278.
        assert!((a - 0.1278).abs() < 1e-3, "a = {a}");
        assert!(WaterModel::Tip3p.m_site_weight().is_none());
    }

    #[test]
    fn m_site_lies_on_bisector_at_tip4p_distance() {
        let [o, h1, h2] = tip4p_molecule(Vec3::new(1.0, 2.0, 3.0));
        let a = WaterModel::Tip4p.m_site_weight().unwrap();
        let m = m_site_position(o, h1, h2, a);
        let bisector = (h1 - o) + (h2 - o);
        let d = WaterModel::Tip4p.params().m_site.unwrap();
        assert!(((m - o).norm() - d).abs() < 1e-12, "|OM| = {}", (m - o).norm());
        assert!((m - o).cross(bisector).norm() < 1e-12, "M-site is off the bisector");
        // The M-site points from O towards the hydrogens.
        assert!((m - o).dot(bisector) > 0.0);
    }

    #[test]
    fn spread_conserves_force_and_torque_for_random_geometry() {
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        let a = WaterModel::Tip4p.m_site_weight().unwrap();
        for _ in 0..200 {
            let atoms = [random_vec(&mut state), random_vec(&mut state), random_vec(&mut state)];
            let f_m = random_vec(&mut state);
            let m = m_site_position(atoms[0], atoms[1], atoms[2], a);
            let spread = spread_m_force(f_m, a);

            let total = spread.iter().fold(Vec3::ZERO, |acc, f| acc + *f);
            assert!((total - f_m).norm() < 1e-12, "force not conserved");

            let torque =
                atoms.iter().zip(&spread).fold(Vec3::ZERO, |acc, (x, f)| acc + x.cross(*f));
            assert!((torque - m.cross(f_m)).norm() < 1e-12, "torque not conserved");
        }
    }

    #[test]
    fn tip4p_sites_spread_conserves_totals_in_a_small_box() {
        // Two TIP4P waters in a 3 nm box. The M-sites repel each other with a
        // 1/r^2 force, standing in for the M-site Coulomb term.
        let mut positions = Vec::new();
        let mut molecules = Vec::new();
        for (k, origin) in
            [Vec3::new(0.5, 0.5, 0.5), Vec3::new(1.5, 1.2, 1.1)].into_iter().enumerate()
        {
            positions.extend(tip4p_molecule(origin));
            molecules.push([3 * k, 3 * k + 1, 3 * k + 2]);
        }
        let sites = Tip4pSites::new(WaterModel::Tip4p, molecules, positions.len()).unwrap();
        assert!((sites.weight() - WaterModel::Tip4p.m_site_weight().unwrap()).abs() < 1e-15);
        let m_pos = sites.positions(&positions);
        assert_eq!(m_pos.len(), 2);

        let r = m_pos[0] - m_pos[1];
        let f_on_0 = r * (1.0 / r.norm().powi(3));
        let m_forces = [f_on_0, -f_on_0];
        let mut forces = vec![Vec3::ZERO; positions.len()];
        sites.spread_forces(&m_forces, &mut forces);

        let net = forces.iter().fold(Vec3::ZERO, |acc, f| acc + *f);
        assert!(net.norm() < 1e-12, "net force {net:?}");
        let torque =
            positions.iter().zip(&forces).fold(Vec3::ZERO, |acc, (x, f)| acc + x.cross(*f));
        let torque_m =
            m_pos.iter().zip(&m_forces).fold(Vec3::ZERO, |acc, (x, f)| acc + x.cross(*f));
        assert!((torque - torque_m).norm() < 1e-12);
    }

    #[test]
    fn tip4p_sites_reject_models_without_m_site_and_bad_indices() {
        assert!(Tip4pSites::new(WaterModel::SpcE, vec![[0, 1, 2]], 3).is_err());
        assert!(Tip4pSites::new(WaterModel::Tip4p, vec![[0, 1, 3]], 3).is_err());
        assert!(Tip4pSites::new(WaterModel::Tip4p, vec![[0, 1, 2]], 3).is_ok());
    }
}
