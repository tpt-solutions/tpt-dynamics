//! Root-mean-square deviation after optimal rigid-body superposition.
//!
//! The optimal rotation is found with Horn's (1987) quaternion method. After
//! both point sets are centred on their (weighted) centroids, the largest
//! eigenvalue `λ` of a 4×4 symmetric matrix built from the cross-covariance
//! gives the best achievable overlap, and the minimal squared deviation is
//! `(Σw|x|² + Σw|y|² − 2λ) / Σw`. Eigenvalues are computed with cyclic Jacobi
//! rotations, which are simple and accurate for a matrix this small.

use tpt_dynamics_core::Vec3;

use super::{check_weights, AnalysisError, Snapshot};

/// RMSD between `mobile` and `reference` after the best rigid superposition.
///
/// `weights` (one per atom, non-negative) gives the weighting, for example
/// masses. Without weights every atom counts equally. Translation and rotation
/// are removed; reflections are not.
///
/// # Errors
/// Returns an error if the inputs are empty, have different lengths, or the
/// weights are invalid or sum to zero.
pub fn rmsd(
    reference: &[Vec3],
    mobile: &[Vec3],
    weights: Option<&[f64]>,
) -> Result<f64, AnalysisError> {
    let n = reference.len();
    if n == 0 {
        return Err(AnalysisError::EmptyInput("RMSD needs at least one atom"));
    }
    if mobile.len() != n {
        return Err(AnalysisError::AtomCountMismatch { expected: n, found: mobile.len() });
    }
    check_weights(weights, n)?;
    let w_of = |i: usize| weights.map_or(1.0, |w| w[i]);
    let total: f64 = (0..n).map(w_of).sum();
    if total <= 0.0 {
        return Err(AnalysisError::EmptyInput("weights sum to zero"));
    }

    let centroid = |pts: &[Vec3]| {
        let mut c = Vec3::ZERO;
        for (i, &p) in pts.iter().enumerate() {
            c += p * w_of(i);
        }
        c / total
    };
    let cr = centroid(reference);
    let cm = centroid(mobile);

    // Cross-covariance S[a][b] = Σ w x_a y_b with x from reference, y from mobile.
    let mut s = [[0.0f64; 3]; 3];
    let mut sum_sq = 0.0;
    for i in 0..n {
        let w = w_of(i);
        let x = reference[i] - cr;
        let y = mobile[i] - cm;
        sum_sq += w * (x.norm_sq() + y.norm_sq());
        let xa = [x.x, x.y, x.z];
        let ya = [y.x, y.y, y.z];
        for (a, row) in s.iter_mut().enumerate() {
            for (b, cell) in row.iter_mut().enumerate() {
                *cell += w * xa[a] * ya[b];
            }
        }
    }

    let lambda_max = max_eigenvalue(&horn_matrix(&s));
    let msd = ((sum_sq - 2.0 * lambda_max) / total).max(0.0);
    Ok(msd.sqrt())
}

/// RMSD of each frame against a fixed reference (same atom count and selection).
///
/// Frames are consumed one at a time. The result holds one value per frame.
///
/// # Errors
/// Returns the first error from [`rmsd`] or a snapshot whose atom count
/// differs from the reference.
pub fn rmsd_series<I>(
    reference: &[Vec3],
    frames: I,
    weights: Option<&[f64]>,
) -> Result<Vec<f64>, AnalysisError>
where
    I: IntoIterator<Item = Snapshot>,
{
    frames.into_iter().map(|f| rmsd(reference, &f.positions, weights)).collect()
}

/// Horn's 4×4 symmetric matrix whose top eigenvalue equals the maximum of
/// `tr(R S)` over rotations `R`.
fn horn_matrix(s: &[[f64; 3]; 3]) -> [[f64; 4]; 4] {
    let (sxx, sxy, sxz) = (s[0][0], s[0][1], s[0][2]);
    let (syx, syy, syz) = (s[1][0], s[1][1], s[1][2]);
    let (szx, szy, szz) = (s[2][0], s[2][1], s[2][2]);
    [
        [sxx + syy + szz, syz - szy, szx - sxz, sxy - syx],
        [syz - szy, sxx - syy - szz, sxy + syx, szx + sxz],
        [szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy],
        [sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz],
    ]
}

/// Largest eigenvalue of a symmetric matrix by cyclic Jacobi rotations.
fn max_eigenvalue(input: &[[f64; 4]; 4]) -> f64 {
    const N: usize = 4;
    let mut a = *input;
    for _sweep in 0..100 {
        let off: f64 = (0..N)
            .flat_map(|p| (0..N).filter(move |&q| q != p).map(move |q| (p, q)))
            .map(|(p, q)| a[p][q] * a[p][q])
            .sum();
        let diag: f64 = (0..N).map(|i| a[i][i] * a[i][i]).sum();
        if off <= f64::EPSILON * f64::EPSILON * diag.max(f64::MIN_POSITIVE) {
            break;
        }
        for p in 0..N {
            for q in (p + 1)..N {
                if a[p][q] == 0.0 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                // A <- J^T A J with J = [[c, s], [-s, c]] acting on (p, q).
                for row in &mut a {
                    let (akp, akq) = (row[p], row[q]);
                    row[p] = c * akp - s * akq;
                    row[q] = s * akp + c * akq;
                }
                let (rp, rq) = (a[p], a[q]);
                a[p] = std::array::from_fn(|k| c * rp[k] - s * rq[k]);
                a[q] = std::array::from_fn(|k| s * rp[k] + c * rq[k]);
            }
        }
    }
    (0..N).map(|i| a[i][i]).fold(f64::NEG_INFINITY, f64::max)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    /// Rotation about an axis (unit vector `u`) by `angle` radians (Rodrigues).
    fn rotate(p: Vec3, u: Vec3, angle: f64) -> Vec3 {
        let (s, c) = angle.sin_cos();
        p * c + u.cross(p) * s + u * (u.dot(p) * (1.0 - c))
    }

    fn asymmetric_cloud(rng: &mut Rng, n: usize) -> Vec<Vec3> {
        (0..n)
            .map(|_| Vec3::new(2.0 * rng.uniform() - 1.0, rng.uniform() - 0.5, 2.0 * rng.uniform()))
            .collect()
    }

    #[test]
    fn rotated_and_translated_copy_has_zero_rmsd() {
        let mut rng = Rng::seed_from_u64(11);
        let reference = asymmetric_cloud(&mut rng, 25);
        let axis = Vec3::new(1.0, 2.0, -0.5).normalized().unwrap();
        let t = Vec3::new(3.0, -4.0, 0.25);
        let mobile: Vec<Vec3> = reference.iter().map(|&p| rotate(p, axis, 1.234) + t).collect();
        let d = rmsd(&reference, &mobile, None).unwrap();
        assert!(d < 1e-7, "rmsd of an exact rigid copy was {d}");
    }

    #[test]
    fn rigid_copy_with_weights_has_zero_rmsd() {
        let mut rng = Rng::seed_from_u64(5);
        let reference = asymmetric_cloud(&mut rng, 12);
        let weights: Vec<f64> = (0..12u32).map(|i| 1.0 + f64::from(i)).collect();
        let mobile: Vec<Vec3> = reference
            .iter()
            .map(|&p| rotate(p, Vec3::new(0.0, 0.0, 1.0), 2.0) + Vec3::new(-1.0, 0.5, 0.0))
            .collect();
        let d = rmsd(&reference, &mobile, Some(&weights)).unwrap();
        assert!(d < 1e-7, "weighted rmsd was {d}");
    }

    #[test]
    fn rmsd_matches_known_displacement() {
        // Moving one point of a square by 2 units: after superposition the
        // minimum is not simply 2/sqrt(4), so compare to a brute-force scan
        // over the rotation angle in the plane.
        let reference = vec![
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(0.0, -1.0, 0.0),
        ];
        let mut mobile = reference.clone();
        mobile[0] = Vec3::new(3.0, 0.0, 0.0);
        let d = rmsd(&reference, &mobile, None).unwrap();
        let mut best = f64::INFINITY;
        for k in 0..20_000 {
            let th = std::f64::consts::TAU * f64::from(k) / 20_000.0;
            // Centroid of the mobile set is (0.5, 0, 0) shifted; superposition
            // includes translation, so compare centred sets.
            let cm = Vec3::new(0.5, 0.0, 0.0);
            let cr = Vec3::ZERO;
            let mut sum = 0.0;
            for (r, m) in reference.iter().zip(&mobile) {
                let x = *r - cr;
                let y = rotate(*m - cm, Vec3::new(0.0, 0.0, 1.0), th);
                sum += (x - y).norm_sq();
            }
            best = best.min((sum / 4.0).sqrt());
        }
        assert!((d - best).abs() < 1e-3, "rmsd {d} vs brute force {best}");
        assert!(d > 0.0);
    }

    #[test]
    fn series_matches_single_calls() {
        let mut rng = Rng::seed_from_u64(3);
        let reference = asymmetric_cloud(&mut rng, 8);
        let sim_box = crate::pbc::SimBox::new(10.0, 10.0, 10.0).unwrap();
        let frames: Vec<Snapshot> = (0..4)
            .map(|k| Snapshot {
                positions: reference
                    .iter()
                    .map(|&p| p + Vec3::new(0.01 * f64::from(k), 0.0, 0.0))
                    .collect(),
                sim_box,
            })
            .collect();
        let series = rmsd_series(&reference, frames.clone(), None).unwrap();
        assert_eq!(series.len(), 4);
        for (f, s) in frames.iter().zip(series) {
            assert_eq!(s, rmsd(&reference, &f.positions, None).unwrap());
        }
    }

    #[test]
    fn rejects_bad_input() {
        let a = vec![Vec3::ZERO; 3];
        assert!(rmsd(&[], &[], None).is_err());
        assert!(matches!(
            rmsd(&a, &a[..2], None),
            Err(AnalysisError::AtomCountMismatch { expected: 3, found: 2 })
        ));
        assert!(rmsd(&a, &a, Some(&[0.0, 0.0, 0.0])).is_err());
        assert!(rmsd(&a, &a, Some(&[1.0, -1.0, 1.0])).is_err());
    }
}
