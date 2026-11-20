//! Two-dimensional histogram with per-axis under/overflow and weights.
//!
//! The histogram stores an `(nx + 2) × (ny + 2)` grid: the inner `nx × ny` cells
//! are the regular bins, and the outer ring holds the under/overflow of each
//! axis (corners included). Bins are half-open, as in
//! [`crate::hist1d::Histogram1D`].

use crate::hist1d::{locate, uniform_edges, validate_edges, Histogram1D, HistogramError, Slot};

/// A 2D histogram with weighted fills.
#[derive(Clone, Debug, PartialEq)]
pub struct Histogram2D {
    pub(crate) xedges: Vec<f64>,
    pub(crate) yedges: Vec<f64>,
    /// Row-major grid, index `ix * (ny + 2) + iy` with `ix`, `iy` in `0..n + 2`
    /// (0 = underflow, `n + 1` = overflow).
    pub(crate) sumw: Vec<f64>,
    pub(crate) sumw2: Vec<f64>,
    pub(crate) entries: u64,
    pub(crate) nan_fills: u64,
}

impl Histogram2D {
    /// Create a histogram with `nx × ny` uniform bins over `[xlo, xhi) × [ylo, yhi)`.
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

    /// Create a histogram with variable bin edges on each axis.
    ///
    /// Each axis needs at least two finite, strictly increasing edges.
    pub fn variable(xedges: &[f64], yedges: &[f64]) -> Result<Self, HistogramError> {
        validate_edges(xedges)?;
        validate_edges(yedges)?;
        Ok(Self::from_edges_unchecked(xedges.to_vec(), yedges.to_vec()))
    }

    /// Build a histogram from raw state. Used by serialisation.
    pub(crate) fn from_parts(
        xedges: Vec<f64>,
        yedges: Vec<f64>,
        sumw: Vec<f64>,
        sumw2: Vec<f64>,
        entries: u64,
        nan_fills: u64,
    ) -> Result<Self, HistogramError> {
        validate_edges(&xedges)?;
        validate_edges(&yedges)?;
        let cells = (xedges.len() + 1) * (yedges.len() + 1);
        if sumw.len() != cells || sumw2.len() != cells {
            return Err(HistogramError::InvalidEdges);
        }
        Ok(Self { xedges, yedges, sumw, sumw2, entries, nan_fills })
    }

    fn from_edges_unchecked(xedges: Vec<f64>, yedges: Vec<f64>) -> Self {
        let cells = (xedges.len() + 1) * (yedges.len() + 1);
        Self {
            xedges,
            yedges,
            sumw: vec![0.0; cells],
            sumw2: vec![0.0; cells],
            entries: 0,
            nan_fills: 0,
        }
    }

    /// Build a histogram from bin contents and per-bin variances (`sum of w^2`).
    ///
    /// `contents` and `variances` are row-major over regular bins, index
    /// `ix * ny + iy`, with `ny = yedges.len() - 1`. Flow cells start at zero and
    /// `entries` is zero for imported data.
    ///
    /// # Errors
    /// Returns [`HistogramError::InvalidEdges`] for bad edges, and
    /// [`HistogramError::InvalidContents`] for a length mismatch or a non-finite
    /// content, or a negative or non-finite variance.
    pub fn from_bins(
        xedges: &[f64],
        yedges: &[f64],
        contents: &[f64],
        variances: &[f64],
    ) -> Result<Self, HistogramError> {
        validate_edges(xedges)?;
        validate_edges(yedges)?;
        let (nx, ny) = (xedges.len() - 1, yedges.len() - 1);
        if contents.len() != nx * ny
            || variances.len() != nx * ny
            || contents.iter().any(|c| !c.is_finite())
            || variances.iter().any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(HistogramError::InvalidContents);
        }
        let mut h = Self::from_edges_unchecked(xedges.to_vec(), yedges.to_vec());
        for ix in 0..nx {
            for iy in 0..ny {
                let idx = h.flat(ix + 1, iy + 1);
                h.sumw[idx] = contents[ix * ny + iy];
                h.sumw2[idx] = variances[ix * ny + iy];
            }
        }
        Ok(h)
    }

    /// Statistical variance `Σw²` of regular bin `(ix, iy)`.
    ///
    /// # Panics
    /// Panics if `ix >= nbins_x()` or `iy >= nbins_y()`.
    #[inline]
    pub fn bin_variance(&self, ix: usize, iy: usize) -> f64 {
        assert!(ix < self.nbins_x() && iy < self.nbins_y(), "bin index out of range");
        self.sumw2[self.flat(ix + 1, iy + 1)]
    }

    /// Number of regular bins along x.
    #[inline]
    pub fn nbins_x(&self) -> usize {
        self.xedges.len() - 1
    }

    /// Number of regular bins along y.
    #[inline]
    pub fn nbins_y(&self) -> usize {
        self.yedges.len() - 1
    }

    /// Bin edges along x, length `nbins_x() + 1`.
    #[inline]
    pub fn x_edges(&self) -> &[f64] {
        &self.xedges
    }

    /// Bin edges along y, length `nbins_y() + 1`.
    #[inline]
    pub fn y_edges(&self) -> &[f64] {
        &self.yedges
    }

    #[inline]
    fn flat(&self, ix: usize, iy: usize) -> usize {
        ix * (self.yedges.len() + 1) + iy
    }

    /// Grid index along one axis: 0 = underflow, `1..=n` = bins, `n + 1` = overflow.
    #[inline]
    fn axis_index(slot: Slot, nbins: usize) -> usize {
        match slot {
            Slot::Under => 0,
            Slot::In(i) => i + 1,
            Slot::Over => nbins + 1,
        }
    }

    #[inline]
    fn cell(&self, x: Slot, y: Slot) -> usize {
        self.flat(Self::axis_index(x, self.nbins_x()), Self::axis_index(y, self.nbins_y()))
    }

    /// Fill with unit weight.
    #[inline]
    pub fn fill(&mut self, x: f64, y: f64) {
        self.fill_weighted(x, y, 1.0);
    }

    /// Fill with weight `w`. A `NaN` coordinate or non-finite weight is counted in
    /// [`Self::nan_fills`] and otherwise ignored.
    pub fn fill_weighted(&mut self, x: f64, y: f64, w: f64) {
        if x.is_nan() || y.is_nan() || !w.is_finite() {
            self.nan_fills += 1;
            return;
        }
        self.entries += 1;
        let idx = self.cell(locate(&self.xedges, x), locate(&self.yedges, y));
        self.sumw[idx] += w;
        self.sumw2[idx] += w * w;
    }

    /// Weighted content of the cell addressed by per-axis positions. Use
    /// [`Slot::Under`] and [`Slot::Over`] to read the flow cells.
    pub fn content(&self, x: Slot, y: Slot) -> f64 {
        self.sumw[self.cell(x, y)]
    }

    /// Statistical error `sqrt(Σw²)` of the cell addressed by per-axis positions.
    pub fn content_error(&self, x: Slot, y: Slot) -> f64 {
        self.sumw2[self.cell(x, y)].sqrt()
    }

    /// Weighted content of regular bin `(ix, iy)`.
    ///
    /// # Panics
    /// Panics if `ix >= nbins_x()` or `iy >= nbins_y()`.
    #[inline]
    pub fn bin(&self, ix: usize, iy: usize) -> f64 {
        assert!(ix < self.nbins_x() && iy < self.nbins_y(), "bin index out of range");
        self.sumw[self.flat(ix + 1, iy + 1)]
    }

    /// Statistical error of regular bin `(ix, iy)`.
    ///
    /// # Panics
    /// Panics if `ix >= nbins_x()` or `iy >= nbins_y()`.
    #[inline]
    pub fn bin_error(&self, ix: usize, iy: usize) -> f64 {
        assert!(ix < self.nbins_x() && iy < self.nbins_y(), "bin index out of range");
        self.sumw2[self.flat(ix + 1, iy + 1)].sqrt()
    }

    /// Sum of weights inside the histogram range (excludes all flow cells).
    pub fn integral(&self) -> f64 {
        let (nx, ny) = (self.nbins_x(), self.nbins_y());
        let mut total = 0.0;
        for ix in 1..=nx {
            for iy in 1..=ny {
                total += self.sumw[self.flat(ix, iy)];
            }
        }
        total
    }

    /// Number of accepted fill calls (excludes rejected NaN / non-finite inputs).
    #[inline]
    pub fn entries(&self) -> u64 {
        self.entries
    }

    /// Number of fills rejected because a coordinate or weight was not finite.
    #[inline]
    pub fn nan_fills(&self) -> u64 {
        self.nan_fills
    }

    /// Project onto the x axis, summing over the in-range y bins.
    ///
    /// The result's under/overflow are the x-flows of the in-range y bins. Its
    /// mean and standard deviation are computed from bin centres, because the 2D
    /// histogram does not keep per-axis moments.
    pub fn project_x(&self) -> Histogram1D {
        let (nx, ny) = (self.nbins_x(), self.nbins_y());
        let sum_row =
            |ix: usize, grid: &[f64]| -> f64 { (1..=ny).map(|iy| grid[self.flat(ix, iy)]).sum() };
        let sumw: Vec<f64> = (1..=nx).map(|ix| sum_row(ix, &self.sumw)).collect();
        let sumw2: Vec<f64> = (1..=nx).map(|ix| sum_row(ix, &self.sumw2)).collect();
        build_projection(
            &self.xedges,
            sumw,
            sumw2,
            [
                sum_row(0, &self.sumw),
                sum_row(nx + 1, &self.sumw),
                sum_row(0, &self.sumw2),
                sum_row(nx + 1, &self.sumw2),
            ],
        )
    }

    /// Project onto the y axis, summing over the in-range x bins.
    ///
    /// See [`Self::project_x`] for the flow and moment conventions.
    pub fn project_y(&self) -> Histogram1D {
        let (nx, ny) = (self.nbins_x(), self.nbins_y());
        let sum_col =
            |iy: usize, grid: &[f64]| -> f64 { (1..=nx).map(|ix| grid[self.flat(ix, iy)]).sum() };
        let sumw: Vec<f64> = (1..=ny).map(|iy| sum_col(iy, &self.sumw)).collect();
        let sumw2: Vec<f64> = (1..=ny).map(|iy| sum_col(iy, &self.sumw2)).collect();
        build_projection(
            &self.yedges,
            sumw,
            sumw2,
            [
                sum_col(0, &self.sumw),
                sum_col(ny + 1, &self.sumw),
                sum_col(0, &self.sumw2),
                sum_col(ny + 1, &self.sumw2),
            ],
        )
    }

    /// Add the contents of `other` into `self`. Both axes must have identical edges.
    ///
    /// Like [`Histogram1D::merge`](crate::hist1d::Histogram1D::merge), this is
    /// associative and commutative up to floating-point rounding.
    pub fn merge(&mut self, other: &Self) -> Result<(), HistogramError> {
        if self.xedges != other.xedges || self.yedges != other.yedges {
            return Err(HistogramError::BinningMismatch);
        }
        for (a, b) in self.sumw.iter_mut().zip(&other.sumw) {
            *a += b;
        }
        for (a, b) in self.sumw2.iter_mut().zip(&other.sumw2) {
            *a += b;
        }
        self.entries += other.entries;
        self.nan_fills += other.nan_fills;
        Ok(())
    }
}

/// Build a 1D histogram from one axis and its summed contents.
fn build_projection(edges: &[f64], sumw: Vec<f64>, sumw2: Vec<f64>, flow: [f64; 4]) -> Histogram1D {
    let mut moments = [0.0_f64; 3];
    for (i, &c) in sumw.iter().enumerate() {
        let x = f64::midpoint(edges[i], edges[i + 1]);
        moments[0] += c;
        moments[1] += c * x;
        moments[2] += c * x * x;
    }
    // The edges were validated when this histogram was built.
    Histogram1D::from_parts(edges.to_vec(), sumw, sumw2, flow, 0, 0, moments)
        .expect("axis of a valid Histogram2D is a valid 1D axis")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn uniform_cells_and_flows() {
        let mut h = Histogram2D::uniform(2, 0.0, 2.0, 2, 0.0, 2.0).unwrap();
        h.fill(0.5, 0.5);
        h.fill(1.5, 1.5);
        h.fill(1.5, 1.5);
        h.fill(-1.0, 0.5); // x underflow
        h.fill(0.5, 2.0); // y overflow
        h.fill(5.0, -5.0); // corner: x overflow, y underflow
        assert_eq!(h.bin(0, 0), 1.0);
        assert_eq!(h.bin(1, 1), 2.0);
        assert_eq!(h.content(Slot::Under, Slot::In(0)), 1.0);
        assert_eq!(h.content(Slot::In(0), Slot::Over), 1.0);
        assert_eq!(h.content(Slot::Over, Slot::Under), 1.0);
        assert_eq!(h.integral(), 3.0);
        assert_eq!(h.entries(), 6);
    }

    #[test]
    fn upper_edge_is_overflow() {
        let mut h = Histogram2D::uniform(1, 0.0, 1.0, 1, 0.0, 1.0).unwrap();
        h.fill(1.0, 0.5);
        h.fill(0.5, 1.0);
        assert_eq!(h.integral(), 0.0);
        assert_eq!(h.content(Slot::Over, Slot::In(0)), 1.0);
        assert_eq!(h.content(Slot::In(0), Slot::Over), 1.0);
    }

    #[test]
    fn variable_binning_and_weights() {
        let mut h = Histogram2D::variable(&[0.0, 1.0, 10.0], &[-1.0, 1.0]).unwrap();
        h.fill_weighted(0.5, 0.0, 2.0);
        h.fill_weighted(5.0, 0.0, 3.0);
        assert_eq!(h.bin(0, 0), 2.0);
        assert_eq!(h.bin(1, 0), 3.0);
        assert!((h.bin_error(0, 0) - 2.0).abs() < 1e-12);
        assert_eq!(h.nbins_x(), 2);
        assert_eq!(h.nbins_y(), 1);
    }

    #[test]
    fn rejects_bad_input() {
        assert_eq!(Histogram2D::uniform(0, 0.0, 1.0, 1, 0.0, 1.0), Err(HistogramError::NoBins));
        assert_eq!(
            Histogram2D::variable(&[0.0, 1.0], &[1.0, 0.0]),
            Err(HistogramError::InvalidEdges)
        );
        let mut h = Histogram2D::uniform(1, 0.0, 1.0, 1, 0.0, 1.0).unwrap();
        h.fill(f64::NAN, 0.5);
        h.fill_weighted(0.5, 0.5, f64::INFINITY);
        assert_eq!(h.nan_fills(), 2);
        assert_eq!(h.entries(), 0);
    }

    #[test]
    fn projections_match_sums() {
        let mut h = Histogram2D::uniform(4, 0.0, 4.0, 2, 0.0, 2.0).unwrap();
        for ix in 0..4_u32 {
            for iy in 0..2_u32 {
                h.fill_weighted(f64::from(ix) + 0.5, f64::from(iy) + 0.5, f64::from(ix + 1));
            }
        }
        h.fill(-1.0, 0.5); // x underflow at in-range y
        let px = h.project_x();
        assert_eq!(px.contents(), &[2.0, 4.0, 6.0, 8.0]);
        assert_eq!(px.underflow(), 1.0);
        assert_eq!(px.integral(), 20.0);
        let py = h.project_y();
        assert_eq!(py.contents(), &[10.0, 10.0]);
        assert_eq!(py.underflow(), 0.0);
        assert_eq!(py.sum_weights_squared()[0], 1.0 + 4.0 + 9.0 + 16.0);
    }

    #[test]
    fn merge_matches_single_fill_and_rejects_mismatch() {
        let pts: Vec<(f64, f64)> = (0..500)
            .map(|i| (f64::from(i % 17) * 0.7 - 1.0, f64::from((i * 31) % 13) * 0.2))
            .collect();
        let mut whole = Histogram2D::uniform(5, 0.0, 10.0, 4, 0.0, 2.0).unwrap();
        let mut a = whole.clone();
        let mut b = whole.clone();
        for (i, &(x, y)) in pts.iter().enumerate() {
            whole.fill(x, y);
            if i % 2 == 0 {
                a.fill(x, y);
            } else {
                b.fill(x, y);
            }
        }
        b.merge(&a).unwrap();
        assert_eq!(b.entries(), whole.entries());
        for (x, y) in b.sumw.iter().zip(&whole.sumw) {
            assert!((x - y).abs() < 1e-9);
        }
        let other = Histogram2D::uniform(6, 0.0, 10.0, 4, 0.0, 2.0).unwrap();
        assert_eq!(b.merge(&other), Err(HistogramError::BinningMismatch));
    }
}

#[cfg(test)]
mod from_bins_tests {
    use super::*;

    #[test]
    fn from_bins_is_row_major_and_matches_fills() {
        let h = Histogram2D::from_bins(
            &[0.0, 1.0, 2.0],
            &[0.0, 1.0, 2.0, 3.0],
            &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
            &[1.0, 4.0, 9.0, 16.0, 25.0, 36.0],
        )
        .unwrap();
        assert_eq!(h.bin(0, 2), 3.0);
        assert_eq!(h.bin(1, 0), 4.0);
        assert_eq!(h.bin_variance(1, 2), 36.0);
        assert_eq!(h.integral(), 21.0);
        assert_eq!(h.entries(), 0);
    }

    #[test]
    fn from_bins_rejects_length_mismatch_and_negative_variance() {
        assert_eq!(
            Histogram2D::from_bins(&[0.0, 1.0], &[0.0, 1.0], &[1.0, 2.0], &[1.0]),
            Err(HistogramError::InvalidContents)
        );
        assert_eq!(
            Histogram2D::from_bins(&[0.0, 1.0], &[0.0, 1.0], &[1.0], &[-1.0]),
            Err(HistogramError::InvalidContents)
        );
    }
}
