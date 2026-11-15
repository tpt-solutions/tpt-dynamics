//! Profile histograms: the mean of a second quantity in each bin.
//!
//! * [`Profile1D`] fills `(x, y)` pairs and reports the mean of `y` in each
//!   `x`-bin, the classic `TProfile`.
//! * [`Profile2D`] fills `(x, y, z)` triples and reports the mean of `z` in each
//!   `(x, y)` cell, the classic `TProfile2D`.
//!
//! Each cell keeps weighted sums `Σw`, `Σwy`, `Σwy²` and `Σw²` plus the raw
//! entry count. The mean is `Σwy / Σw`. Its error is
//! `sqrt(var / n_eff)`, where `var` is the weighted (population) variance of `y`
//! in the cell and `n_eff = (Σw)² / Σw²`. For unit weights this is the usual
//! standard error of the mean. Fills outside the axis range are counted as
//! outside weight and do not enter any cell.

use crate::hist1d::{locate, uniform_edges, validate_edges, HistogramError, Slot};

/// Accumulated sums for one profile cell.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ProfileCell {
    /// `Σw`.
    pub sumw: f64,
    /// `Σwy`.
    pub sumwy: f64,
    /// `Σwy²`.
    pub sumwy2: f64,
    /// `Σw²`.
    pub sumw2: f64,
    /// Number of fills accepted into this cell.
    pub entries: u64,
}

impl ProfileCell {
    #[inline]
    fn fill(&mut self, y: f64, w: f64) {
        self.sumw += w;
        self.sumwy += w * y;
        self.sumwy2 += w * y * y;
        self.sumw2 += w * w;
        self.entries += 1;
    }

    fn merge(&mut self, other: &Self) {
        self.sumw += other.sumw;
        self.sumwy += other.sumwy;
        self.sumwy2 += other.sumwy2;
        self.sumw2 += other.sumw2;
        self.entries += other.entries;
    }

    /// Mean of `y`, or `None` for an empty cell.
    #[inline]
    pub fn mean(&self) -> Option<f64> {
        (self.entries > 0 && self.sumw != 0.0).then(|| self.sumwy / self.sumw)
    }

    /// Error on the mean, or `None` with fewer than two entries.
    pub fn mean_error(&self) -> Option<f64> {
        if self.entries < 2 || self.sumw == 0.0 || self.sumw2 == 0.0 {
            return None;
        }
        let mean = self.sumwy / self.sumw;
        let var = (self.sumwy2 / self.sumw - mean * mean).max(0.0);
        let neff = self.sumw * self.sumw / self.sumw2;
        Some((var / neff).sqrt())
    }
}

/// One-dimensional profile histogram of `y` versus `x`.
#[derive(Clone, Debug, PartialEq)]
pub struct Profile1D {
    pub(crate) edges: Vec<f64>,
    pub(crate) cells: Vec<ProfileCell>,
    pub(crate) underflow: f64,
    pub(crate) overflow: f64,
    pub(crate) nan_fills: u64,
}

impl Profile1D {
    /// Create a profile with `nbins` uniform `x`-bins over `[lo, hi)`.
    pub fn uniform(nbins: usize, lo: f64, hi: f64) -> Result<Self, HistogramError> {
        Ok(Self::from_edges_unchecked(uniform_edges(nbins, lo, hi)?))
    }

    /// Create a profile with variable `x`-bin edges.
    pub fn variable(edges: &[f64]) -> Result<Self, HistogramError> {
        validate_edges(edges)?;
        Ok(Self::from_edges_unchecked(edges.to_vec()))
    }

    pub(crate) fn from_parts(
        edges: Vec<f64>,
        cells: Vec<ProfileCell>,
        underflow: f64,
        overflow: f64,
        nan_fills: u64,
    ) -> Result<Self, HistogramError> {
        validate_edges(&edges)?;
        if cells.len() != edges.len() - 1 {
            return Err(HistogramError::InvalidEdges);
        }
        Ok(Self { edges, cells, underflow, overflow, nan_fills })
    }

    fn from_edges_unchecked(edges: Vec<f64>) -> Self {
        let n = edges.len() - 1;
        Self {
            edges,
            cells: vec![ProfileCell::default(); n],
            underflow: 0.0,
            overflow: 0.0,
            nan_fills: 0,
        }
    }

    /// Number of `x`-bins.
    #[inline]
    pub fn nbins(&self) -> usize {
        self.cells.len()
    }

    /// Bin edges in `x`, length `nbins() + 1`.
    #[inline]
    pub fn edges(&self) -> &[f64] {
        &self.edges
    }

    /// Raw sums for bin `i`.
    ///
    /// # Panics
    /// Panics if `i >= nbins()`.
    #[inline]
    pub fn cell(&self, i: usize) -> &ProfileCell {
        &self.cells[i]
    }

    /// Fill with `x`, `y` and unit weight.
    #[inline]
    pub fn fill(&mut self, x: f64, y: f64) {
        self.fill_weighted(x, y, 1.0);
    }

    /// Fill with `x`, `y` and weight `w`. Non-finite `x`, `y` or `w` are counted in
    /// [`Self::nan_fills`].
    pub fn fill_weighted(&mut self, x: f64, y: f64, w: f64) {
        if x.is_nan() || !y.is_finite() || !w.is_finite() {
            self.nan_fills += 1;
            return;
        }
        match locate(&self.edges, x) {
            Slot::Under => self.underflow += w,
            Slot::Over => self.overflow += w,
            Slot::In(i) => self.cells[i].fill(y, w),
        }
    }

    /// Mean of `y` in bin `i`, or `None` if the bin is empty.
    #[inline]
    pub fn mean(&self, i: usize) -> Option<f64> {
        self.cells[i].mean()
    }

    /// Error on the mean in bin `i`, or `None` with fewer than two entries.
    #[inline]
    pub fn mean_error(&self, i: usize) -> Option<f64> {
        self.cells[i].mean_error()
    }

    /// Number of entries accepted into bin `i`.
    #[inline]
    pub fn bin_entries(&self, i: usize) -> u64 {
        self.cells[i].entries
    }

    /// Weight of fills with `x` below the first edge.
    #[inline]
    pub fn underflow(&self) -> f64 {
        self.underflow
    }

    /// Weight of fills with `x` at or above the last edge.
    #[inline]
    pub fn overflow(&self) -> f64 {
        self.overflow
    }

    /// Number of fills rejected because a value or weight was not finite.
    #[inline]
    pub fn nan_fills(&self) -> u64 {
        self.nan_fills
    }

    /// Add the contents of `other` into `self`. Requires identical `x` edges.
    pub fn merge(&mut self, other: &Self) -> Result<(), HistogramError> {
        if self.edges != other.edges {
            return Err(HistogramError::BinningMismatch);
        }
        for (a, b) in self.cells.iter_mut().zip(&other.cells) {
            a.merge(b);
        }
        self.underflow += other.underflow;
        self.overflow += other.overflow;
        self.nan_fills += other.nan_fills;
        Ok(())
    }
}

/// Two-dimensional profile histogram of `z` versus `(x, y)`.
#[derive(Clone, Debug, PartialEq)]
pub struct Profile2D {
    pub(crate) xedges: Vec<f64>,
    pub(crate) yedges: Vec<f64>,
    /// Row-major, index `ix * ny + iy`.
    pub(crate) cells: Vec<ProfileCell>,
    pub(crate) outside: f64,
    pub(crate) nan_fills: u64,
}

impl Profile2D {
    /// Create a profile with `nx × ny` uniform cells over `[xlo, xhi) × [ylo, yhi)`.
    pub fn uniform(
        nx: usize,
        xlo: f64,
        xhi: f64,
        ny: usize,
        ylo: f64,
        yhi: f64,
    ) -> Result<Self, HistogramError> {
        Ok(Self::from_edges_unchecked(uniform_edges(nx, xlo, xhi)?, uniform_edges(ny, ylo, yhi)?))
    }

    /// Create a profile with variable edges on each axis.
    pub fn variable(xedges: &[f64], yedges: &[f64]) -> Result<Self, HistogramError> {
        validate_edges(xedges)?;
        validate_edges(yedges)?;
        Ok(Self::from_edges_unchecked(xedges.to_vec(), yedges.to_vec()))
    }

    pub(crate) fn from_parts(
        xedges: Vec<f64>,
        yedges: Vec<f64>,
        cells: Vec<ProfileCell>,
        outside: f64,
        nan_fills: u64,
    ) -> Result<Self, HistogramError> {
        validate_edges(&xedges)?;
        validate_edges(&yedges)?;
        if cells.len() != (xedges.len() - 1) * (yedges.len() - 1) {
            return Err(HistogramError::InvalidEdges);
        }
        Ok(Self { xedges, yedges, cells, outside, nan_fills })
    }

    fn from_edges_unchecked(xedges: Vec<f64>, yedges: Vec<f64>) -> Self {
        let n = (xedges.len() - 1) * (yedges.len() - 1);
        Self { xedges, yedges, cells: vec![ProfileCell::default(); n], outside: 0.0, nan_fills: 0 }
    }

    /// Number of bins along x.
    #[inline]
    pub fn nbins_x(&self) -> usize {
        self.xedges.len() - 1
    }

    /// Number of bins along y.
    #[inline]
    pub fn nbins_y(&self) -> usize {
        self.yedges.len() - 1
    }

    /// Bin edges along x.
    #[inline]
    pub fn x_edges(&self) -> &[f64] {
        &self.xedges
    }

    /// Bin edges along y.
    #[inline]
    pub fn y_edges(&self) -> &[f64] {
        &self.yedges
    }

    #[inline]
    fn index(&self, ix: usize, iy: usize) -> usize {
        assert!(ix < self.nbins_x() && iy < self.nbins_y(), "bin index out of range");
        ix * self.nbins_y() + iy
    }

    /// Raw sums for cell `(ix, iy)`.
    ///
    /// # Panics
    /// Panics if the indices are out of range.
    #[inline]
    pub fn cell(&self, ix: usize, iy: usize) -> &ProfileCell {
        &self.cells[self.index(ix, iy)]
    }

    /// Fill with `(x, y, z)` and unit weight.
    #[inline]
    pub fn fill(&mut self, x: f64, y: f64, z: f64) {
        self.fill_weighted(x, y, z, 1.0);
    }

    /// Fill with `(x, y, z)` and weight `w`. Non-finite `z` or `w`, or a `NaN`
    /// coordinate, are counted in [`Self::nan_fills`]. Fills outside the grid add
    /// to [`Self::outside`].
    pub fn fill_weighted(&mut self, x: f64, y: f64, z: f64, w: f64) {
        if x.is_nan() || y.is_nan() || !z.is_finite() || !w.is_finite() {
            self.nan_fills += 1;
            return;
        }
        match (locate(&self.xedges, x), locate(&self.yedges, y)) {
            (Slot::In(ix), Slot::In(iy)) => {
                let idx = ix * self.nbins_y() + iy;
                self.cells[idx].fill(z, w);
            }
            _ => self.outside += w,
        }
    }

    /// Mean of `z` in cell `(ix, iy)`, or `None` if it is empty.
    #[inline]
    pub fn mean(&self, ix: usize, iy: usize) -> Option<f64> {
        self.cell(ix, iy).mean()
    }

    /// Error on the mean in cell `(ix, iy)`, or `None` with fewer than two entries.
    #[inline]
    pub fn mean_error(&self, ix: usize, iy: usize) -> Option<f64> {
        self.cell(ix, iy).mean_error()
    }

    /// Weight of fills that fell outside the grid.
    #[inline]
    pub fn outside(&self) -> f64 {
        self.outside
    }

    /// Number of fills rejected because a value or weight was not finite.
    #[inline]
    pub fn nan_fills(&self) -> u64 {
        self.nan_fills
    }

    /// Add the contents of `other` into `self`. Both axes must have identical edges.
    pub fn merge(&mut self, other: &Self) -> Result<(), HistogramError> {
        if self.xedges != other.xedges || self.yedges != other.yedges {
            return Err(HistogramError::BinningMismatch);
        }
        for (a, b) in self.cells.iter_mut().zip(&other.cells) {
            a.merge(b);
        }
        self.outside += other.outside;
        self.nan_fills += other.nan_fills;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mean_of_y_per_x_bin() {
        let mut p = Profile1D::uniform(2, 0.0, 2.0).unwrap();
        p.fill(0.5, 1.0);
        p.fill(0.5, 3.0);
        p.fill(1.5, 10.0);
        assert_eq!(p.mean(0), Some(2.0));
        assert_eq!(p.mean(1), Some(10.0));
        assert_eq!(p.bin_entries(0), 2);
        assert_eq!(p.mean_error(1), None);
    }

    #[test]
    fn mean_error_matches_standard_error() {
        // y = 0, 2, 4, 6: mean 3, population variance 5, n = 4 -> error sqrt(5/4).
        let mut p = Profile1D::uniform(1, 0.0, 1.0).unwrap();
        for y in [0.0, 2.0, 4.0, 6.0] {
            p.fill(0.5, y);
        }
        assert!((p.mean_error(0).unwrap() - (5.0_f64 / 4.0).sqrt()).abs() < 1e-12);
    }

    #[test]
    fn weights_and_flows() {
        let mut p = Profile1D::variable(&[0.0, 1.0, 2.0]).unwrap();
        p.fill_weighted(0.5, 1.0, 1.0);
        p.fill_weighted(0.5, 4.0, 3.0);
        p.fill(-1.0, 5.0);
        p.fill(2.0, 5.0);
        p.fill(f64::NAN, 5.0);
        p.fill(0.5, f64::INFINITY);
        assert_eq!(p.mean(0), Some(13.0 / 4.0));
        assert_eq!(p.underflow(), 1.0);
        assert_eq!(p.overflow(), 1.0);
        assert_eq!(p.nan_fills(), 2);
    }

    #[test]
    fn profile1d_merge() {
        let data: Vec<(f64, f64)> =
            (0..300).map(|i| (f64::from(i % 10) * 0.1, f64::from(i))).collect();
        let mut whole = Profile1D::uniform(5, 0.0, 1.0).unwrap();
        let mut a = whole.clone();
        let mut b = whole.clone();
        for (i, &(x, y)) in data.iter().enumerate() {
            whole.fill(x, y);
            if i < 100 {
                a.fill(x, y);
            } else {
                b.fill(x, y);
            }
        }
        b.merge(&a).unwrap();
        for i in 0..5 {
            assert_eq!(b.bin_entries(i), whole.bin_entries(i));
            assert!((b.mean(i).unwrap() - whole.mean(i).unwrap()).abs() < 1e-9);
        }
        let other = Profile1D::uniform(4, 0.0, 1.0).unwrap();
        assert_eq!(b.merge(&other), Err(HistogramError::BinningMismatch));
    }

    #[test]
    fn profile2d_cells_and_outside() {
        let mut p = Profile2D::uniform(2, 0.0, 2.0, 2, 0.0, 2.0).unwrap();
        p.fill(0.5, 0.5, 4.0);
        p.fill(0.5, 0.5, 8.0);
        p.fill(1.5, 1.5, -1.0);
        p.fill(2.5, 0.5, 1.0);
        p.fill(0.5, -0.5, 1.0);
        assert_eq!(p.mean(0, 0), Some(6.0));
        assert_eq!(p.mean(1, 1), Some(-1.0));
        assert_eq!(p.mean(0, 1), None);
        assert_eq!(p.outside(), 2.0);
    }

    #[test]
    fn profile2d_merge_rejects_mismatch() {
        let mut a = Profile2D::uniform(2, 0.0, 1.0, 2, 0.0, 1.0).unwrap();
        let mut b = Profile2D::uniform(2, 0.0, 1.0, 2, 0.0, 1.0).unwrap();
        a.fill(0.1, 0.1, 1.0);
        b.fill(0.1, 0.1, 3.0);
        a.merge(&b).unwrap();
        assert_eq!(a.mean(0, 0), Some(2.0));
        let c = Profile2D::uniform(2, 0.0, 1.0, 3, 0.0, 1.0).unwrap();
        assert_eq!(a.merge(&c), Err(HistogramError::BinningMismatch));
    }
}
