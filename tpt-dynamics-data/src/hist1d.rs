//! One-dimensional histogram with weights and merge support.
//!
//! Bins are half-open intervals `[edge[i], edge[i+1])`. Values below the first
//! edge go to the underflow counter; values at or above the last edge go to the
//! overflow counter. `NaN` fills are rejected and counted separately.

use core::fmt;

/// Errors from histogram construction or merging.
#[derive(Clone, Debug, PartialEq)]
pub enum HistogramError {
    /// Zero bins were requested.
    NoBins,
    /// Bin edges are not finite, or are not strictly increasing.
    InvalidEdges,
    /// Two histograms with different binning were merged.
    BinningMismatch,
    /// Bin contents or variances do not match the edges, or are not finite.
    InvalidContents,
}

impl fmt::Display for HistogramError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoBins => write!(f, "histogram must have at least one bin"),
            Self::InvalidEdges => {
                write!(f, "bin edges must be finite and strictly increasing")
            }
            Self::BinningMismatch => {
                write!(f, "cannot merge histograms with different binning")
            }
            Self::InvalidContents => {
                write!(f, "bin contents and variances must match the edges and be finite")
            }
        }
    }
}

impl std::error::Error for HistogramError {}

/// A 1D histogram with weighted fills.
#[derive(Clone, Debug, PartialEq)]
pub struct Histogram1D {
    pub(crate) edges: Vec<f64>,
    pub(crate) sumw: Vec<f64>,
    pub(crate) sumw2: Vec<f64>,
    pub(crate) underflow: f64,
    pub(crate) overflow: f64,
    pub(crate) underflow_w2: f64,
    pub(crate) overflow_w2: f64,
    pub(crate) entries: u64,
    pub(crate) nan_fills: u64,
    // Moments over in-range fills, used for mean / std dev.
    pub(crate) m_w: f64,
    pub(crate) m_wx: f64,
    pub(crate) m_wx2: f64,
}

/// Position of a value relative to a binned axis (used by `Histogram2D::content`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Slot {
    /// Below the first edge.
    Under,
    /// Inside bin `i`.
    In(usize),
    /// At or above the last edge.
    Over,
}

/// Locate `x` on the axis defined by `edges` (ascending, length >= 2).
/// `x` must not be NaN.
#[inline]
pub(crate) fn locate(edges: &[f64], x: f64) -> Slot {
    if x < edges[0] {
        Slot::Under
    } else if x >= edges[edges.len() - 1] {
        Slot::Over
    } else {
        // Index of the first edge strictly greater than x, minus one.
        Slot::In(edges.partition_point(|&e| e <= x) - 1)
    }
}

/// Check that an axis definition has at least one bin and strictly increasing,
/// finite edges.
pub(crate) fn validate_edges(edges: &[f64]) -> Result<(), HistogramError> {
    if edges.len() < 2 {
        return Err(HistogramError::NoBins);
    }
    if edges.iter().any(|e| !e.is_finite()) || edges.windows(2).any(|w| w[0] >= w[1]) {
        return Err(HistogramError::InvalidEdges);
    }
    Ok(())
}

/// Build `nbins` uniform edges over `[lo, hi]` with the last edge pinned to `hi`.
pub(crate) fn uniform_edges(nbins: usize, lo: f64, hi: f64) -> Result<Vec<f64>, HistogramError> {
    if nbins == 0 {
        return Err(HistogramError::NoBins);
    }
    if !(lo.is_finite() && hi.is_finite() && lo < hi) {
        return Err(HistogramError::InvalidEdges);
    }
    let width = (hi - lo) / nbins as f64;
    let mut edges: Vec<f64> = (0..=nbins).map(|i| lo + width * i as f64).collect();
    // Pin the last edge exactly so that values equal to `hi` overflow cleanly.
    edges[nbins] = hi;
    Ok(edges)
}

impl Histogram1D {
    /// Create a histogram with `nbins` uniform bins over `[lo, hi)`.
    pub fn uniform(nbins: usize, lo: f64, hi: f64) -> Result<Self, HistogramError> {
        Ok(Self::from_edges_unchecked(uniform_edges(nbins, lo, hi)?))
    }

    /// Create a histogram with variable bin edges.
    ///
    /// `edges` must contain at least two finite, strictly increasing values.
    pub fn variable(edges: &[f64]) -> Result<Self, HistogramError> {
        validate_edges(edges)?;
        Ok(Self::from_edges_unchecked(edges.to_vec()))
    }

    /// Build a histogram from bin contents and per-bin variances (`sum of w^2`).
    ///
    /// Used to import histograms from other libraries. Flow, fill counts and
    /// moments start at zero, so [`Self::mean`] and [`Self::std_dev`] return `None`
    /// until new fills arrive. `entries` is zero for imported data.
    ///
    /// # Errors
    /// Returns [`HistogramError::InvalidEdges`] for bad edges, and
    /// [`HistogramError::InvalidContents`] if the slices do not have `edges.len() - 1`
    /// entries, or if a content is not finite or a variance is negative or not finite.
    pub fn from_bins(
        edges: &[f64],
        contents: &[f64],
        variances: &[f64],
    ) -> Result<Self, HistogramError> {
        validate_edges(edges)?;
        let n = edges.len() - 1;
        if contents.len() != n
            || variances.len() != n
            || contents.iter().any(|c| !c.is_finite())
            || variances.iter().any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(HistogramError::InvalidContents);
        }
        Self::from_parts(
            edges.to_vec(),
            contents.to_vec(),
            variances.to_vec(),
            [0.0; 4],
            0,
            0,
            [0.0; 3],
        )
    }

    /// Build a histogram from raw state. Used by serialisation and projections.
    ///
    /// Lengths of `sumw` and `sumw2` must equal `edges.len() - 1`; the caller is
    /// responsible for the consistency of the remaining values.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn from_parts(
        edges: Vec<f64>,
        sumw: Vec<f64>,
        sumw2: Vec<f64>,
        flow: [f64; 4],
        entries: u64,
        nan_fills: u64,
        moments: [f64; 3],
    ) -> Result<Self, HistogramError> {
        validate_edges(&edges)?;
        if sumw.len() != edges.len() - 1 || sumw2.len() != sumw.len() {
            return Err(HistogramError::InvalidEdges);
        }
        Ok(Self {
            edges,
            sumw,
            sumw2,
            underflow: flow[0],
            overflow: flow[1],
            underflow_w2: flow[2],
            overflow_w2: flow[3],
            entries,
            nan_fills,
            m_w: moments[0],
            m_wx: moments[1],
            m_wx2: moments[2],
        })
    }

    /// Raw moments `(Σw, Σwx, Σwx²)` over in-range fills.
    pub(crate) fn moments(&self) -> [f64; 3] {
        [self.m_w, self.m_wx, self.m_wx2]
    }

    /// Weighted under/overflow and their squared-weight sums, in the order
    /// `[underflow, overflow, underflow_w2, overflow_w2]`.
    pub(crate) fn flow(&self) -> [f64; 4] {
        [self.underflow, self.overflow, self.underflow_w2, self.overflow_w2]
    }

    fn from_edges_unchecked(edges: Vec<f64>) -> Self {
        let n = edges.len() - 1;
        Self {
            edges,
            sumw: vec![0.0; n],
            sumw2: vec![0.0; n],
            underflow: 0.0,
            overflow: 0.0,
            underflow_w2: 0.0,
            overflow_w2: 0.0,
            entries: 0,
            nan_fills: 0,
            m_w: 0.0,
            m_wx: 0.0,
            m_wx2: 0.0,
        }
    }

    /// Number of bins (excluding under/overflow).
    #[inline]
    pub fn nbins(&self) -> usize {
        self.sumw.len()
    }

    /// Bin edges, length `nbins() + 1`.
    #[inline]
    pub fn edges(&self) -> &[f64] {
        &self.edges
    }

    /// Fill with unit weight.
    #[inline]
    pub fn fill(&mut self, x: f64) {
        self.fill_weighted(x, 1.0);
    }

    /// Fill with weight `w`. Non-finite weights are treated like `NaN` values and
    /// counted in [`Self::nan_fills`].
    pub fn fill_weighted(&mut self, x: f64, w: f64) {
        if x.is_nan() || !w.is_finite() {
            self.nan_fills += 1;
            return;
        }
        self.entries += 1;
        let w2 = w * w;
        match locate(&self.edges, x) {
            Slot::Under => {
                self.underflow += w;
                self.underflow_w2 += w2;
            }
            Slot::Over => {
                self.overflow += w;
                self.overflow_w2 += w2;
            }
            Slot::In(idx) => {
                self.sumw[idx] += w;
                self.sumw2[idx] += w2;
                self.m_w += w;
                self.m_wx += w * x;
                self.m_wx2 += w * x * x;
            }
        }
    }

    /// Weighted content of bin `i`.
    #[inline]
    pub fn bin(&self, i: usize) -> f64 {
        self.sumw[i]
    }

    /// Statistical error of bin `i`: `sqrt(Σw²)`.
    #[inline]
    pub fn bin_error(&self, i: usize) -> f64 {
        self.sumw2[i].sqrt()
    }

    /// Centre of bin `i`.
    #[inline]
    pub fn bin_center(&self, i: usize) -> f64 {
        f64::midpoint(self.edges[i], self.edges[i + 1])
    }

    /// Weighted content of all bins, in order.
    #[inline]
    pub fn contents(&self) -> &[f64] {
        &self.sumw
    }

    /// Sum of squared weights per bin, in order.
    #[inline]
    pub fn sum_weights_squared(&self) -> &[f64] {
        &self.sumw2
    }

    /// Weighted underflow.
    #[inline]
    pub fn underflow(&self) -> f64 {
        self.underflow
    }

    /// Weighted overflow.
    #[inline]
    pub fn overflow(&self) -> f64 {
        self.overflow
    }

    /// Sum of squared weights in underflow.
    #[inline]
    pub fn underflow_sum_w2(&self) -> f64 {
        self.underflow_w2
    }

    /// Sum of squared weights in overflow.
    #[inline]
    pub fn overflow_sum_w2(&self) -> f64 {
        self.overflow_w2
    }

    /// Number of fill calls accepted (excludes NaN / non-finite weights).
    #[inline]
    pub fn entries(&self) -> u64 {
        self.entries
    }

    /// Number of fills rejected because the value or weight was not finite.
    #[inline]
    pub fn nan_fills(&self) -> u64 {
        self.nan_fills
    }

    /// Sum of weights inside the histogram range (excludes under/overflow).
    #[inline]
    pub fn integral(&self) -> f64 {
        self.sumw.iter().sum()
    }

    /// Weighted mean of in-range fills, or `None` if there are none.
    pub fn mean(&self) -> Option<f64> {
        (self.m_w != 0.0).then(|| self.m_wx / self.m_w)
    }

    /// Weighted standard deviation of in-range fills, or `None` if fewer than
    /// two effective fills.
    pub fn std_dev(&self) -> Option<f64> {
        if self.m_w == 0.0 {
            return None;
        }
        let mean = self.m_wx / self.m_w;
        let var = self.m_wx2 / self.m_w - mean * mean;
        // Guard against tiny negative values from floating-point cancellation.
        Some(var.max(0.0).sqrt())
    }

    /// Approximate quantile `p ∈ [0, 1]` of the in-range distribution, using
    /// linear interpolation within bins. Returns `None` if the histogram is empty
    /// or `p` is outside `[0, 1]`.
    pub fn quantile(&self, p: f64) -> Option<f64> {
        if !(0.0..=1.0).contains(&p) {
            return None;
        }
        let total = self.integral();
        if total <= 0.0 {
            return None;
        }
        let target = p * total;
        let mut cum = 0.0;
        for (i, &c) in self.sumw.iter().enumerate() {
            if c > 0.0 && cum + c >= target {
                let frac = (target - cum) / c;
                return Some(self.edges[i] + frac * (self.edges[i + 1] - self.edges[i]));
            }
            cum += c;
        }
        Some(self.edges[self.edges.len() - 1])
    }

    /// Add the contents of `other` into `self`. Requires identical binning.
    ///
    /// Merging is associative and commutative (up to floating-point rounding), so
    /// partial histograms from parallel or distributed workers can be combined in
    /// any order.
    pub fn merge(&mut self, other: &Self) -> Result<(), HistogramError> {
        if self.edges != other.edges {
            return Err(HistogramError::BinningMismatch);
        }
        for i in 0..self.sumw.len() {
            self.sumw[i] += other.sumw[i];
            self.sumw2[i] += other.sumw2[i];
        }
        self.underflow += other.underflow;
        self.overflow += other.overflow;
        self.underflow_w2 += other.underflow_w2;
        self.overflow_w2 += other.overflow_w2;
        self.entries += other.entries;
        self.nan_fills += other.nan_fills;
        self.m_w += other.m_w;
        self.m_wx += other.m_wx;
        self.m_wx2 += other.m_wx2;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-12 * a.abs().max(b.abs()).max(1.0)
    }

    #[test]
    fn rejects_bad_binning() {
        assert_eq!(Histogram1D::uniform(0, 0.0, 1.0), Err(HistogramError::NoBins));
        assert_eq!(Histogram1D::uniform(4, 1.0, 1.0), Err(HistogramError::InvalidEdges));
        assert_eq!(Histogram1D::uniform(4, 0.0, f64::NAN), Err(HistogramError::InvalidEdges));
        assert_eq!(Histogram1D::variable(&[0.0]), Err(HistogramError::NoBins));
        assert_eq!(Histogram1D::variable(&[0.0, 2.0, 1.0]), Err(HistogramError::InvalidEdges));
    }

    #[test]
    fn fixed_binning_and_bounds() {
        let mut h = Histogram1D::uniform(10, 0.0, 10.0).unwrap();
        h.fill(0.0); // first bin, inclusive lower edge
        h.fill(9.999); // last bin
        h.fill(10.0); // overflow (upper edge exclusive)
        h.fill(-0.001); // underflow
        assert_eq!(h.bin(0), 1.0);
        assert_eq!(h.bin(9), 1.0);
        assert_eq!(h.overflow(), 1.0);
        assert_eq!(h.underflow(), 1.0);
        assert_eq!(h.integral(), 2.0);
        assert_eq!(h.entries(), 4);
    }

    #[test]
    fn variable_binning_lookup() {
        let mut h = Histogram1D::variable(&[0.0, 1.0, 4.0, 10.0]).unwrap();
        h.fill(0.5);
        h.fill(1.0); // lands in [1, 4)
        h.fill(3.9);
        h.fill(7.0);
        assert_eq!(h.bin(0), 1.0);
        assert_eq!(h.bin(1), 2.0);
        assert_eq!(h.bin(2), 1.0);
    }

    #[test]
    fn weights_and_errors() {
        let mut h = Histogram1D::uniform(2, 0.0, 2.0).unwrap();
        h.fill_weighted(0.5, 2.0);
        h.fill_weighted(0.5, 3.0);
        assert_eq!(h.bin(0), 5.0);
        assert_eq!(h.sum_weights_squared()[0], 13.0);
        assert!(close(h.bin_error(0), 13f64.sqrt()));
    }

    #[test]
    fn nan_and_nonfinite_are_rejected() {
        let mut h = Histogram1D::uniform(2, 0.0, 2.0).unwrap();
        h.fill(f64::NAN);
        h.fill_weighted(0.5, f64::INFINITY);
        assert_eq!(h.nan_fills(), 2);
        assert_eq!(h.entries(), 0);
        assert_eq!(h.integral(), 0.0);
    }

    #[test]
    fn mean_std_and_quantile() {
        let mut h = Histogram1D::uniform(100, 0.0, 100.0).unwrap();
        for i in 0..100 {
            h.fill(f64::from(i) + 0.5);
        }
        assert!(close(h.mean().unwrap(), 50.0));
        // Uniform on [0,100): σ = 100/√12 ≈ 28.87 (bin-centre approximation).
        assert!((h.std_dev().unwrap() - 100.0 / 12f64.sqrt()).abs() < 0.1);
        assert!((h.quantile(0.5).unwrap() - 50.0).abs() < 1.0);
        assert!(h.quantile(1.5).is_none());
        assert!(Histogram1D::uniform(3, 0.0, 1.0).unwrap().quantile(0.5).is_none());
    }

    #[test]
    fn merge_matches_single_fill_sequence() {
        let data: Vec<f64> =
            (0_u32..1000).map(|i| f64::from((i * 7919) % 1000) / 10.0 - 5.0).collect();

        let mut whole = Histogram1D::uniform(20, 0.0, 100.0).unwrap();
        for &x in &data {
            whole.fill(x);
        }

        let mut parts = [
            Histogram1D::uniform(20, 0.0, 100.0).unwrap(),
            Histogram1D::uniform(20, 0.0, 100.0).unwrap(),
            Histogram1D::uniform(20, 0.0, 100.0).unwrap(),
        ];
        for (i, &x) in data.iter().enumerate() {
            parts[i % 3].fill(x);
        }
        // Merge in a different order than filled to exercise commutativity.
        let mut merged = parts[2].clone();
        merged.merge(&parts[0]).unwrap();
        merged.merge(&parts[1]).unwrap();

        assert_eq!(merged.entries(), whole.entries());
        assert_eq!(merged.underflow(), whole.underflow());
        assert_eq!(merged.overflow(), whole.overflow());
        for i in 0..whole.nbins() {
            assert!(close(merged.bin(i), whole.bin(i)));
        }
        assert!(close(merged.mean().unwrap(), whole.mean().unwrap()));
    }

    #[test]
    fn merge_rejects_different_binning() {
        let mut a = Histogram1D::uniform(10, 0.0, 1.0).unwrap();
        let b = Histogram1D::uniform(11, 0.0, 1.0).unwrap();
        assert_eq!(a.merge(&b), Err(HistogramError::BinningMismatch));
    }
}

#[cfg(test)]
mod from_bins_tests {
    use super::*;

    #[test]
    fn from_bins_round_trips_contents_and_variances() {
        let h = Histogram1D::from_bins(&[0.0, 1.0, 3.0], &[2.0, 5.0], &[4.0, 25.0]).unwrap();
        assert_eq!(h.contents(), &[2.0, 5.0]);
        assert_eq!(h.sum_weights_squared(), &[4.0, 25.0]);
        assert!((h.bin_error(0) - 2.0).abs() < 1e-12);
        assert_eq!(h.entries(), 0);
        assert_eq!(h.mean(), None);
    }

    #[test]
    fn from_bins_rejects_bad_input() {
        assert_eq!(
            Histogram1D::from_bins(&[0.0, 1.0], &[1.0, 2.0], &[1.0]),
            Err(HistogramError::InvalidContents)
        );
        assert_eq!(
            Histogram1D::from_bins(&[0.0, 1.0], &[1.0], &[-1.0]),
            Err(HistogramError::InvalidContents)
        );
        assert_eq!(
            Histogram1D::from_bins(&[0.0, 1.0], &[f64::NAN], &[1.0]),
            Err(HistogramError::InvalidContents)
        );
        assert_eq!(
            Histogram1D::from_bins(&[1.0, 0.0], &[1.0], &[1.0]),
            Err(HistogramError::InvalidEdges)
        );
    }

    #[test]
    fn imported_histogram_accepts_further_fills() {
        let mut h = Histogram1D::from_bins(&[0.0, 1.0], &[1.0], &[1.0]).unwrap();
        h.fill_weighted(0.5, 2.0);
        assert_eq!(h.bin(0), 3.0);
        assert_eq!(h.sum_weights_squared()[0], 5.0);
    }
}
