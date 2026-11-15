//! Radial distribution function g(r).
//!
//! For a selection of `N` atoms in a frame with volume `V`, the ideal-gas
//! expectation for the number of distinct pairs whose separation falls in a
//! shell `[r_k, r_{k+1})` is
//!
//! ```text
//! E_k = (N(N−1)/2) · (4π/3)(r_{k+1}³ − r_k³) / V
//! ```
//!
//! and `g(r_k) = Σ_frames count_k / Σ_frames E_k`. The separations use the
//! minimum-image convention, so `r_max` may be at most half the smallest box
//! edge in every frame. Each distinct unordered pair is counted once.

use tpt_dynamics_core::Vec3;

use super::{AnalysisError, Snapshot};

/// A binned radial distribution function.
#[derive(Clone, Debug, PartialEq)]
pub struct Rdf {
    /// Bin centres in nm.
    pub r: Vec<f64>,
    /// g(r) per bin. Bins with no ideal-gas expectation are 0.
    pub g: Vec<f64>,
    /// Raw pair counts per bin, summed over frames.
    pub counts: Vec<u64>,
    /// Frames accumulated.
    pub n_frames: usize,
}

/// Accumulates g(r) over the frames of a trajectory.
///
/// `selection` restricts the atoms to these indices (`None` means all atoms).
/// Only pairs inside the selection are counted.
///
/// # Errors
/// Returns an error if there are no frames, `n_bins` is zero, `r_max` is not
/// positive and finite, `r_max` exceeds half the smallest box edge of a frame,
/// the selection is empty, or an index is out of range.
pub fn rdf<I>(
    frames: I,
    selection: Option<&[usize]>,
    r_max: f64,
    n_bins: usize,
) -> Result<Rdf, AnalysisError>
where
    I: IntoIterator<Item = Snapshot>,
{
    if n_bins == 0 {
        return Err(AnalysisError::InvalidParameter("n_bins must be at least 1"));
    }
    if !(r_max.is_finite() && r_max > 0.0) {
        return Err(AnalysisError::InvalidParameter("r_max must be finite and > 0"));
    }
    let dr = r_max / n_bins as f64;
    let shell_volume: Vec<f64> = (0..n_bins)
        .map(|k| {
            let lo = k as f64 * dr;
            let hi = lo + dr;
            4.0 * std::f64::consts::PI / 3.0 * (hi * hi * hi - lo * lo * lo)
        })
        .collect();

    let mut counts = vec![0u64; n_bins];
    let mut expected = vec![0.0f64; n_bins];
    let mut n_frames = 0usize;
    for frame in frames {
        let n_atoms = frame.positions.len();
        let sim_box = frame.sim_box;
        if r_max > 0.5 * sim_box.min_edge() {
            return Err(AnalysisError::InvalidParameter(
                "r_max exceeds half the smallest box edge",
            ));
        }
        let members: Vec<usize> = match selection {
            Some(sel) => {
                if let Some(&bad) = sel.iter().find(|&&i| i >= n_atoms) {
                    return Err(AnalysisError::IndexOutOfRange { index: bad, n_atoms });
                }
                sel.to_vec()
            }
            None => (0..n_atoms).collect(),
        };
        let n = members.len();
        if n < 2 {
            return Err(AnalysisError::EmptyInput("RDF needs at least two atoms in the selection"));
        }
        let pts: Vec<Vec3> = members.iter().map(|&i| frame.positions[i]).collect();
        for (i, &pi) in pts.iter().enumerate() {
            for &pj in &pts[i + 1..] {
                let r = sim_box.min_image(pj - pi).norm();
                if r < r_max {
                    let k = ((r / dr) as usize).min(n_bins - 1);
                    counts[k] += 1;
                }
            }
        }
        let pairs = n as f64 * (n as f64 - 1.0) / 2.0;
        let volume = sim_box.volume();
        for (e, &sv) in expected.iter_mut().zip(&shell_volume) {
            *e += pairs * sv / volume;
        }
        n_frames += 1;
    }
    if n_frames == 0 {
        return Err(AnalysisError::EmptyInput("RDF needs at least one frame"));
    }

    let g = counts
        .iter()
        .zip(&expected)
        .map(|(&c, &e)| if e > 0.0 { c as f64 / e } else { 0.0 })
        .collect();
    let r = (0..n_bins).map(|k| (k as f64 + 0.5) * dr).collect();
    Ok(Rdf { r, g, counts, n_frames })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pbc::SimBox;
    use crate::rng::Rng;

    #[test]
    fn ideal_gas_tends_to_one_at_large_r() {
        let edge = 4.0;
        let sim_box = SimBox::new(edge, edge, edge).unwrap();
        let mut rng = Rng::seed_from_u64(2024);
        let frames: Vec<Snapshot> = (0..40)
            .map(|_| {
                let positions = (0..300)
                    .map(|_| {
                        Vec3::new(edge * rng.uniform(), edge * rng.uniform(), edge * rng.uniform())
                    })
                    .collect();
                Snapshot { positions, sim_box }
            })
            .collect();
        let out = rdf(frames, None, 2.0, 20).unwrap();
        assert_eq!(out.n_frames, 40);
        let tail: Vec<f64> =
            out.g.iter().zip(&out.r).filter(|(_, &r)| r > 1.2).map(|(&g, _)| g).collect();
        let mean = tail.iter().sum::<f64>() / tail.len() as f64;
        assert!((mean - 1.0).abs() < 0.03, "mean g(r>1.2) = {mean}");
    }

    #[test]
    fn single_pair_lands_in_expected_bin() {
        let sim_box = SimBox::new(4.0, 4.0, 4.0).unwrap();
        // Pair separated by 0.5 nm; bins of width 0.1 over [0, 1.0).
        let s = Snapshot {
            positions: vec![Vec3::new(1.0, 1.0, 1.0), Vec3::new(1.5, 1.0, 1.0)],
            sim_box,
        };
        let out = rdf(vec![s], None, 1.0, 10).unwrap();
        assert_eq!(out.counts.iter().sum::<u64>(), 1);
        assert_eq!(out.counts[5], 1);
        assert!((out.r[5] - 0.55).abs() < 1e-12);
    }

    #[test]
    fn minimum_image_counts_across_boundary() {
        let sim_box = SimBox::new(4.0, 4.0, 4.0).unwrap();
        let s = Snapshot {
            positions: vec![Vec3::new(0.05, 2.0, 2.0), Vec3::new(3.95, 2.0, 2.0)],
            sim_box,
        };
        let out = rdf(vec![s], None, 0.5, 5).unwrap();
        assert_eq!(out.counts[0], 1, "0.1 nm apart through the boundary");
    }

    #[test]
    fn selection_restricts_pairs() {
        let sim_box = SimBox::new(4.0, 4.0, 4.0).unwrap();
        let s = Snapshot {
            positions: vec![
                Vec3::new(1.0, 1.0, 1.0),
                Vec3::new(1.2, 1.0, 1.0),
                Vec3::new(3.0, 3.0, 3.0),
            ],
            sim_box,
        };
        let out = rdf(vec![s.clone()], Some(&[0, 2]), 1.0, 10).unwrap();
        assert_eq!(out.counts.iter().sum::<u64>(), 0);
        let all = rdf(vec![s], None, 1.0, 10).unwrap();
        assert_eq!(all.counts.iter().sum::<u64>(), 1);
    }

    #[test]
    fn rejects_bad_parameters() {
        let sim_box = SimBox::new(4.0, 4.0, 4.0).unwrap();
        let s = Snapshot { positions: vec![Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0)], sim_box };
        assert!(rdf(vec![s.clone()], None, 2.5, 10).is_err(), "r_max beyond L/2");
        assert!(rdf(vec![s.clone()], None, 1.0, 0).is_err());
        assert!(rdf(vec![s.clone()], Some(&[7]), 1.0, 4).is_err());
        assert!(rdf(Vec::<Snapshot>::new(), None, 1.0, 4).is_err());
    }
}
