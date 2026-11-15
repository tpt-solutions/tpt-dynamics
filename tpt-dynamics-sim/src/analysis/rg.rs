//! Radius of gyration.

use tpt_dynamics_core::Vec3;

use super::{check_weights, AnalysisError, Snapshot};

/// Radius of gyration `Rg = sqrt(Σ mᵢ |rᵢ − r_c|² / Σ mᵢ)` in nm, where `r_c`
/// is the same weighted centre.
///
/// `masses` gives the weights; without them every atom has weight 1 (geometric
/// radius of gyration). Positions are used as given, so molecules must be
/// unwrapped beforehand.
///
/// # Errors
/// Returns an error for an empty selection, a mass slice of the wrong length,
/// or masses that are negative, non-finite or sum to zero.
pub fn radius_of_gyration(
    positions: &[Vec3],
    masses: Option<&[f64]>,
) -> Result<f64, AnalysisError> {
    let n = positions.len();
    if n == 0 {
        return Err(AnalysisError::EmptyInput("radius of gyration needs at least one atom"));
    }
    check_weights(masses, n)?;
    let w_of = |i: usize| masses.map_or(1.0, |m| m[i]);
    let total: f64 = (0..n).map(w_of).sum();
    if total <= 0.0 {
        return Err(AnalysisError::EmptyInput("total mass is zero"));
    }
    let mut centre = Vec3::ZERO;
    for (i, &p) in positions.iter().enumerate() {
        centre += p * w_of(i);
    }
    centre = centre / total;
    let sum: f64 =
        positions.iter().enumerate().map(|(i, &p)| w_of(i) * (p - centre).norm_sq()).sum();
    Ok((sum / total).sqrt())
}

/// Radius of gyration of each frame, in frame order.
///
/// # Errors
/// Returns the first error from [`radius_of_gyration`].
pub fn rg_series<I>(frames: I, masses: Option<&[f64]>) -> Result<Vec<f64>, AnalysisError>
where
    I: IntoIterator<Item = Snapshot>,
{
    frames.into_iter().map(|f| radius_of_gyration(&f.positions, masses)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regular_tetrahedron_matches_analytic_value() {
        // Vertices of a cube's alternate corners: each lies at distance sqrt(3)
        // from the centroid at the origin, so Rg = sqrt(3).
        let pts = [
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(1.0, -1.0, -1.0),
            Vec3::new(-1.0, 1.0, -1.0),
            Vec3::new(-1.0, -1.0, 1.0),
        ];
        let rg = radius_of_gyration(&pts, None).unwrap();
        assert!((rg - 3f64.sqrt()).abs() < 1e-12, "rg = {rg}");
    }

    #[test]
    fn uniform_ring_matches_radius() {
        let n = 360;
        let r = 0.7;
        let pts: Vec<Vec3> = (0..n)
            .map(|k| {
                let th = std::f64::consts::TAU * f64::from(k) / f64::from(n);
                Vec3::new(r * th.cos(), r * th.sin(), 0.0)
            })
            .collect();
        let rg = radius_of_gyration(&pts, None).unwrap();
        assert!((rg - r).abs() < 1e-12, "rg = {rg}");
    }

    #[test]
    fn mass_weighting_moves_centre_of_mass() {
        // Two atoms 2 nm apart with masses 1 and 3: centre at 1.5, so
        // Rg² = (1·1.5² + 3·0.5²)/4 = (2.25 + 0.75)/4 = 0.75.
        let pts = [Vec3::new(0.0, 0.0, 0.0), Vec3::new(2.0, 0.0, 0.0)];
        let rg = radius_of_gyration(&pts, Some(&[1.0, 3.0])).unwrap();
        assert!((rg * rg - 0.75).abs() < 1e-12);
        let geometric = radius_of_gyration(&pts, None).unwrap();
        assert!((geometric - 1.0).abs() < 1e-12);
    }

    #[test]
    fn series_and_errors() {
        let sim_box = crate::pbc::SimBox::new(5.0, 5.0, 5.0).unwrap();
        let frames = vec![
            Snapshot { positions: vec![Vec3::ZERO, Vec3::new(2.0, 0.0, 0.0)], sim_box },
            Snapshot { positions: vec![Vec3::ZERO, Vec3::new(4.0, 0.0, 0.0)], sim_box },
        ];
        let series = rg_series(frames, None).unwrap();
        assert!((series[0] - 1.0).abs() < 1e-12);
        assert!((series[1] - 2.0).abs() < 1e-12);
        assert!(radius_of_gyration(&[], None).is_err());
        assert!(radius_of_gyration(&[Vec3::ZERO], Some(&[0.0])).is_err());
        assert!(radius_of_gyration(&[Vec3::ZERO], Some(&[1.0, 2.0])).is_err());
    }
}
