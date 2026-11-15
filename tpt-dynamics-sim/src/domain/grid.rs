//! Domain grid: an axis-aligned split of an orthorhombic periodic box into
//! `nx × ny × nz` sub-domains.
//!
//! Each axis has its own list of boundaries, so the grid is separable: the
//! boundaries along x do not depend on y or z. Uniform grids use equal
//! spacing. [`DomainGrid::rebalanced`] moves the boundaries to atom-count
//! quantiles along each axis, which equalises the marginal counts per axis.

use rayon::prelude::*;
use tpt_dynamics_core::Vec3;

/// A separable grid of sub-domains over `[0, L_x) × [0, L_y) × [0, L_z)`.
#[derive(Clone, Debug, PartialEq)]
pub struct DomainGrid {
    dims: [usize; 3],
    /// Per axis: `dims[a] + 1` non-decreasing boundaries from `0` to `L_a`.
    bounds: [Vec<f64>; 3],
}

impl DomainGrid {
    /// Equal-width grid with `dims` sub-domains over a box of edge lengths `lengths`.
    ///
    /// # Panics
    /// Panics if any dimension is zero.
    pub fn uniform(dims: [usize; 3], lengths: Vec3) -> Self {
        assert!(dims.iter().all(|&d| d >= 1), "domain grid dimensions must be >= 1");
        let l = [lengths.x, lengths.y, lengths.z];
        let bounds = [0, 1, 2]
            .map(|a| (0..=dims[a]).map(|k| l[a] * k as f64 / dims[a] as f64).collect::<Vec<f64>>());
        Self { dims, bounds }
    }

    /// Number of sub-domains along each axis.
    pub fn dims(&self) -> [usize; 3] {
        self.dims
    }

    /// Total number of sub-domains.
    pub fn n_domains(&self) -> usize {
        self.dims.iter().product()
    }

    /// Linear index of the sub-domain at grid coordinates `c` (x fastest).
    pub fn index(&self, c: [usize; 3]) -> usize {
        c[0] + self.dims[0] * (c[1] + self.dims[1] * c[2])
    }

    /// Grid coordinates of linear sub-domain index `d`.
    pub fn coords(&self, d: usize) -> [usize; 3] {
        let x = d % self.dims[0];
        let y = (d / self.dims[0]) % self.dims[1];
        let z = d / (self.dims[0] * self.dims[1]);
        [x, y, z]
    }

    /// Lower and upper corners `(lo, hi)` of sub-domain `d`.
    pub fn region(&self, d: usize) -> ([f64; 3], [f64; 3]) {
        let c = self.coords(d);
        let lo = [0, 1, 2].map(|a| self.bounds[a][c[a]]);
        let hi = [0, 1, 2].map(|a| self.bounds[a][c[a] + 1]);
        (lo, hi)
    }

    /// Boundaries along axis `a` (`dims[a] + 1` values).
    pub fn boundaries(&self, a: usize) -> &[f64] {
        &self.bounds[a]
    }

    /// Index of the sub-domain that owns position `p`. The position must already
    /// be wrapped into the box; out-of-range values are clamped to the edge domains.
    pub fn owner(&self, p: Vec3) -> usize {
        let mut c = [0_usize; 3];
        for (a, value) in [p.x, p.y, p.z].into_iter().enumerate() {
            let b = &self.bounds[a];
            // Number of boundaries <= value, minus one, is the domain index. With
            // zero-width domains (repeated boundaries) this picks the non-empty one.
            let k = b.partition_point(|&x| x <= value).saturating_sub(1);
            c[a] = k.min(self.dims[a] - 1);
        }
        self.index(c)
    }

    /// A grid with the same dimensions whose boundaries are atom-count quantiles
    /// along each axis. `lengths` is the box edge.
    pub fn rebalanced(&self, positions: &[Vec3], lengths: Vec3) -> Self {
        let l = [lengths.x, lengths.y, lengths.z];
        let n = positions.len();
        let bounds = [0, 1, 2].map(|a| {
            let mut coords: Vec<f64> = positions.iter().map(|p| [p.x, p.y, p.z][a]).collect();
            coords.par_sort_unstable_by(f64::total_cmp);
            let mut b = Vec::with_capacity(self.dims[a] + 1);
            b.push(0.0);
            for k in 1..self.dims[a] {
                let v = if n == 0 {
                    l[a] * k as f64 / self.dims[a] as f64
                } else {
                    coords[k * n / self.dims[a]]
                };
                // Quantiles of sorted data are non-decreasing; clamp to the box.
                let prev = *b.last().unwrap_or(&0.0);
                b.push(v.clamp(prev, l[a]));
            }
            b.push(l[a]);
            b
        });
        Self { dims: self.dims, bounds }
    }
}

/// Chooses a factorisation `n = nx·ny·nz` that minimises the total boundary
/// surface of the sub-domains for a box with edges `lengths` (less halo traffic).
///
/// # Panics
/// Panics if `n == 0`.
pub fn auto_dims(n: usize, lengths: Vec3) -> [usize; 3] {
    assert!(n >= 1, "need at least one domain");
    let l = [lengths.x, lengths.y, lengths.z];
    let mut best = [n, 1, 1];
    let mut best_cost = f64::INFINITY;
    for nx in 1..=n {
        if n % nx != 0 {
            continue;
        }
        let rest = n / nx;
        for ny in 1..=rest {
            if rest % ny != 0 {
                continue;
            }
            let dims = [nx, ny, rest / ny];
            let w = [0, 1, 2].map(|a| l[a] / dims[a] as f64);
            let cost = w[0] * w[1] + w[1] * w[2] + w[0] * w[2];
            if cost < best_cost - 1e-12 {
                best_cost = cost;
                best = dims;
            }
        }
    }
    best
}

/// Load imbalance `max / mean` of per-domain atom counts (1.0 is perfect).
pub fn imbalance(counts: &[usize]) -> f64 {
    let total: usize = counts.iter().sum();
    if total == 0 || counts.is_empty() {
        return 1.0;
    }
    let mean = total as f64 / counts.len() as f64;
    counts.iter().copied().max().unwrap_or(0) as f64 / mean
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn owner_matches_region_and_index_round_trips() {
        let g = DomainGrid::uniform([3, 2, 4], Vec3::new(3.0, 2.0, 4.0));
        for d in 0..g.n_domains() {
            assert_eq!(g.index(g.coords(d)), d);
            let (lo, hi) = g.region(d);
            let mid =
                Vec3::new(0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1]), 0.5 * (lo[2] + hi[2]));
            assert_eq!(g.owner(mid), d);
        }
    }

    #[test]
    fn rebalance_handles_clustered_atoms_and_zero_width_domains() {
        let lengths = Vec3::new(10.0, 10.0, 10.0);
        // All atoms in the first tenth of x: quantile splits must squeeze the
        // first domains into a narrow slab while keeping boundaries ordered.
        let pos: Vec<Vec3> = (0..100).map(|i| Vec3::new(0.001 * i as f64, 5.0, 5.0)).collect();
        let g = DomainGrid::uniform([4, 1, 1], lengths).rebalanced(&pos, lengths);
        let b = g.boundaries(0);
        assert!(b.windows(2).all(|w| w[0] <= w[1]));
        assert_eq!((b[0], *b.last().unwrap()), (0.0, 10.0));
        let counts = count_per_domain(&g, &pos);
        assert!(imbalance(&counts) < 1.7, "counts {counts:?}");
    }

    #[test]
    fn auto_dims_prefers_compact_factorisation() {
        assert_eq!(auto_dims(8, Vec3::new(1.0, 1.0, 1.0)), [2, 2, 2]);
        let d = auto_dims(20, Vec3::new(100.0, 100.0, 25.0));
        assert_eq!(d.iter().product::<usize>(), 20);
    }

    fn count_per_domain(g: &DomainGrid, pos: &[Vec3]) -> Vec<usize> {
        let mut c = vec![0; g.n_domains()];
        for &p in pos {
            c[g.owner(p)] += 1;
        }
        c
    }
}
