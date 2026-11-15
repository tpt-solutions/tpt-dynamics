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
//! The M-site of TIP4P is a virtual particle whose force must be distributed
//! onto the oxygen and hydrogens. That redistribution is not implemented, so
//! [`WaterModel::Tip4p`] reports `is_supported() == false` and the builders in
//! [`crate::config`] refuse it.

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
}
