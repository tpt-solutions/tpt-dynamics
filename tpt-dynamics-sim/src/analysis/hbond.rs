//! Geometric hydrogen-bond counts.
//!
//! A donor is a pair of atom indices `(D, H)`: the heavy donor atom and the
//! hydrogen bonded to it. Acceptors are single atom indices. A donor–acceptor
//! pair `(D–H, A)` counts as a hydrogen bond when both criteria hold:
//!
//! 1. the donor–acceptor distance `|A − D|` (minimum image) is at most
//!    [`HBondCriteria::max_distance`], and
//! 2. the angle `D–H···A` at the hydrogen (minimum image) is at least
//!    [`HBondCriteria::min_angle_deg`].
//!
//! The defaults (0.35 nm, 150°) follow the usual GROMACS `gmx hbond` geometry.
//! A pair with `A` equal to `D` or `H` is never counted.

use tpt_dynamics_core::Vec3;

use super::{AnalysisError, Snapshot};

/// A hydrogen-bond donor: the heavy atom and the hydrogen attached to it.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Donor {
    /// Index of the donor heavy atom (D).
    pub donor: usize,
    /// Index of the hydrogen bonded to the donor (H).
    pub hydrogen: usize,
}

/// Geometric thresholds for a hydrogen bond.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HBondCriteria {
    /// Maximum donor–acceptor distance in nm.
    pub max_distance: f64,
    /// Minimum donor–hydrogen–acceptor angle in degrees.
    pub min_angle_deg: f64,
}

impl Default for HBondCriteria {
    fn default() -> Self {
        Self { max_distance: 0.35, min_angle_deg: 150.0 }
    }
}

impl HBondCriteria {
    fn validate(&self) -> Result<(), AnalysisError> {
        if !(self.max_distance.is_finite() && self.max_distance > 0.0) {
            return Err(AnalysisError::InvalidParameter("max_distance must be finite and > 0"));
        }
        if !(self.min_angle_deg.is_finite() && (0.0..=180.0).contains(&self.min_angle_deg)) {
            return Err(AnalysisError::InvalidParameter("min_angle_deg must be in [0, 180]"));
        }
        Ok(())
    }
}

/// Counts hydrogen bonds in one snapshot.
///
/// The result is the number of (donor, acceptor) pairs that satisfy both
/// criteria. A hydrogen bonded to two acceptors counts twice.
///
/// # Errors
/// Returns an error for invalid criteria or an index outside the snapshot.
pub fn count_hbonds(
    snapshot: &Snapshot,
    donors: &[Donor],
    acceptors: &[usize],
    criteria: &HBondCriteria,
) -> Result<usize, AnalysisError> {
    criteria.validate()?;
    let n = snapshot.positions.len();
    let check = |i: usize| {
        if i < n {
            Ok(())
        } else {
            Err(AnalysisError::IndexOutOfRange { index: i, n_atoms: n })
        }
    };
    for d in donors {
        check(d.donor)?;
        check(d.hydrogen)?;
    }
    for &a in acceptors {
        check(a)?;
    }

    let sim_box = snapshot.sim_box;
    let pos = &snapshot.positions;
    let cos_limit = criteria.min_angle_deg.to_radians().cos();
    let max_sq = criteria.max_distance * criteria.max_distance;
    let mut count = 0usize;
    for d in donors {
        let dp = pos[d.donor];
        let hp = pos[d.hydrogen];
        // Vector from H to D, used as the first arm of the angle at H.
        let hd = sim_box.min_image(dp - hp);
        if hd.norm_sq() == 0.0 {
            continue;
        }
        for &a in acceptors {
            if a == d.donor || a == d.hydrogen {
                continue;
            }
            let da = sim_box.min_image(pos[a] - dp);
            if da.norm_sq() > max_sq {
                continue;
            }
            let ha = sim_box.min_image(pos[a] - hp);
            if is_aligned(hd, ha, cos_limit) {
                count += 1;
            }
        }
    }
    Ok(count)
}

/// True when the angle between `hd` and `ha` is at least the limit, i.e.
/// `cos(angle) <= cos_limit`. Degenerate zero-length arms never align.
fn is_aligned(hd: Vec3, ha: Vec3, cos_limit: f64) -> bool {
    let denom = hd.norm() * ha.norm();
    if denom == 0.0 {
        return false;
    }
    let cos = (hd.dot(ha) / denom).clamp(-1.0, 1.0);
    cos <= cos_limit
}

/// Hydrogen-bond count for each frame, in frame order.
///
/// # Errors
/// Returns the first error from [`count_hbonds`].
pub fn hbond_series<I>(
    frames: I,
    donors: &[Donor],
    acceptors: &[usize],
    criteria: &HBondCriteria,
) -> Result<Vec<usize>, AnalysisError>
where
    I: IntoIterator<Item = Snapshot>,
{
    frames.into_iter().map(|f| count_hbonds(&f, donors, acceptors, criteria)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pbc::SimBox;

    fn snap(points: &[Vec3], edge: f64) -> Snapshot {
        Snapshot { positions: points.to_vec(), sim_box: SimBox::new(edge, edge, edge).unwrap() }
    }

    #[test]
    fn linear_bond_counts_and_bent_or_long_do_not() {
        let d = Donor { donor: 0, hydrogen: 1 };
        let crit = HBondCriteria::default();
        // D at origin, H at 0.1 nm along +x. Acceptor straight beyond H (angle 180°).
        let linear = snap(&[Vec3::ZERO, Vec3::new(0.1, 0.0, 0.0), Vec3::new(0.3, 0.0, 0.0)], 4.0);
        assert_eq!(count_hbonds(&linear, &[d], &[2], &crit).unwrap(), 1);

        // Acceptor 90° off at H: D–A distance 0.316 nm passes but the angle fails.
        let bent = snap(&[Vec3::ZERO, Vec3::new(0.1, 0.0, 0.0), Vec3::new(0.1, 0.3, 0.0)], 4.0);
        assert_eq!(count_hbonds(&bent, &[d], &[2], &crit).unwrap(), 0);

        // Linear but too far: 0.5 nm exceeds the 0.35 nm cutoff.
        let far = snap(&[Vec3::ZERO, Vec3::new(0.1, 0.0, 0.0), Vec3::new(0.5, 0.0, 0.0)], 4.0);
        assert_eq!(count_hbonds(&far, &[d], &[2], &crit).unwrap(), 0);
    }

    #[test]
    fn angle_threshold_is_respected() {
        let d = Donor { donor: 0, hydrogen: 1 };
        // Angle at H of 160°: H at origin-ish; D at (-0.1, 0, 0) relative to H,
        // A at 0.2 nm in a direction 20° away from +x.
        let h = Vec3::new(0.1, 0.0, 0.0);
        let th = 20f64.to_radians();
        let a = h + Vec3::new(th.cos(), th.sin(), 0.0) * 0.2;
        let s = snap(&[Vec3::ZERO, h, a], 4.0);
        let lax = HBondCriteria { max_distance: 0.35, min_angle_deg: 150.0 };
        assert_eq!(count_hbonds(&s, &[d], &[2], &lax).unwrap(), 1, "160° passes 150°");
        let strict = HBondCriteria { max_distance: 0.35, min_angle_deg: 170.0 };
        assert_eq!(count_hbonds(&s, &[d], &[2], &strict).unwrap(), 0, "160° fails 170°");
    }

    #[test]
    fn uses_minimum_image() {
        let d = Donor { donor: 0, hydrogen: 1 };
        // D at x=0.05, H at x=0.15, acceptor at x=3.95 in a 4 nm box. The acceptor
        // image sits at x=-0.05, behind D on the far side from H, so the angle at H
        // is 0 and no bond is counted.
        let s = snap(
            &[Vec3::new(0.05, 2.0, 2.0), Vec3::new(0.15, 2.0, 2.0), Vec3::new(3.95, 2.0, 2.0)],
            4.0,
        );
        assert_eq!(
            count_hbonds(&s, &[d], &[2], &HBondCriteria::default()).unwrap(),
            0,
            "acceptor image is on the same side of H as D, so the angle at H is 0"
        );
        let s2 = snap(
            &[Vec3::new(0.15, 2.0, 2.0), Vec3::new(0.05, 2.0, 2.0), Vec3::new(3.95, 2.0, 2.0)],
            4.0,
        );
        let d2 = Donor { donor: 0, hydrogen: 1 };
        assert_eq!(
            count_hbonds(&s2, &[d2], &[2], &HBondCriteria::default()).unwrap(),
            1,
            "H at x=0.05, D at 0.15: acceptor image at x=-0.05 lies beyond H, linear"
        );
    }

    #[test]
    fn excludes_self_and_validates_input() {
        let d = Donor { donor: 0, hydrogen: 1 };
        let s = snap(&[Vec3::ZERO, Vec3::new(0.1, 0.0, 0.0)], 4.0);
        assert_eq!(count_hbonds(&s, &[d], &[0, 1], &HBondCriteria::default()).unwrap(), 0);
        assert!(count_hbonds(&s, &[d], &[5], &HBondCriteria::default()).is_err());
        let bad = HBondCriteria { max_distance: -1.0, min_angle_deg: 150.0 };
        assert!(count_hbonds(&s, &[d], &[], &bad).is_err());
    }

    #[test]
    fn series_counts_each_frame() {
        let d = [Donor { donor: 0, hydrogen: 1 }];
        let linear = snap(&[Vec3::ZERO, Vec3::new(0.1, 0.0, 0.0), Vec3::new(0.3, 0.0, 0.0)], 4.0);
        let none = snap(&[Vec3::ZERO, Vec3::new(0.1, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)], 4.0);
        let series = hbond_series(vec![linear, none], &d, &[2], &HBondCriteria::default()).unwrap();
        assert_eq!(series, vec![1, 0]);
    }
}
