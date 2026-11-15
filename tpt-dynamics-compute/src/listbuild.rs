//! SIMD cell-list pair search for neighbour (Verlet) lists.
//!
//! The builders here are pure functions over the structure-of-arrays
//! [`Positions`] view. They return every unordered atom pair `(i, j)` with
//! `i < j` whose (minimum-image, when periodic) separation is strictly less than
//! a radius. A Verlet list is the same search with `radius = cutoff + skin`.
//!
//! Method: atoms are binned into a grid of cells whose edge is at least `radius`,
//! so every partner of an atom lies in its own cell or one of the 26 adjacent
//! cells. Each cell pair is visited once and its atoms are compared four at a time
//! with `tpt-simd` `F64x4` lanes. Cell neighbours are de-duplicated, so grids with
//! one or two cells per periodic axis remain correct.
//!
//! Periodic boxes require `radius <= min(box) / 2`, the same validity limit as the
//! minimum-image convention used by the force kernels.
//!
//! The output is sorted, so results from different builders compare directly.

use tpt_simd::{F64x4, SimdMask};

use crate::simd_kernel::{min_image_v, LANES};
use crate::types::Positions;

/// Returns all unordered pairs `(i, j)`, `i < j`, with separation `< radius`.
///
/// `box_len` enables orthorhombic periodic boundaries with minimum-image distances.
///
/// # Panics
///
/// Panics if `radius` is larger than half of the smallest box edge.
pub fn build_pair_list(
    pos: Positions<'_>,
    radius: f64,
    box_len: Option<[f64; 3]>,
) -> Vec<(u32, u32)> {
    let n = pos.len();
    if n < 2 || radius.is_nan() || radius <= 0.0 {
        return Vec::new();
    }
    if let Some(l) = box_len {
        let min_edge = l.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(
            radius <= 0.5 * min_edge,
            "radius {radius} exceeds half the smallest box edge {min_edge}"
        );
    }

    let grid = Grid::new(pos, radius, box_len);
    let rc2 = F64x4::splat(radius * radius);
    let mut pairs = Vec::new();

    for c in 0..grid.ncells() {
        for d in grid.neighbour_cells(c) {
            if d < c {
                continue;
            }
            let (cs, ce) = grid.range(c);
            let (ds, de) = grid.range(d);
            for i in cs..ce {
                // Within one cell, only partners after i are visited.
                let js = if c == d { i + 1 } else { ds };
                let je = de;
                if js >= je {
                    continue;
                }
                grid.scan(i, js, je, rc2, box_len, &mut pairs);
            }
        }
    }
    pairs.sort_unstable();
    pairs
}

/// Builds a Verlet list: all pairs within `cutoff + skin`.
///
/// Pairs stay valid for as long as no atom has moved more than `skin / 2`.
/// Rebuilding criteria belong to the caller.
///
/// # Panics
///
/// Panics under the same condition as [`build_pair_list`].
pub fn build_verlet_list(
    pos: Positions<'_>,
    cutoff: f64,
    skin: f64,
    box_len: Option<[f64; 3]>,
) -> Vec<(u32, u32)> {
    build_pair_list(pos, cutoff + skin, box_len)
}

/// Atoms binned into cells, with coordinates copied in cell order.
struct Grid {
    dims: [usize; 3],
    periodic: Option<[f64; 3]>,
    /// `start[c]..start[c + 1]` is the range of sorted slots in cell `c`.
    start: Vec<usize>,
    /// Original atom index per sorted slot.
    orig: Vec<u32>,
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
}

impl Grid {
    fn new(pos: Positions<'_>, radius: f64, box_len: Option<[f64; 3]>) -> Self {
        let n = pos.len();
        let (origin, extent) = if let Some(l) = box_len {
            ([0.0; 3], l)
        } else {
            let mut lo = [f64::INFINITY; 3];
            let mut hi = [f64::NEG_INFINITY; 3];
            for a in 0..n {
                let p = [pos.x[a], pos.y[a], pos.z[a]];
                for k in 0..3 {
                    lo[k] = lo[k].min(p[k]);
                    hi[k] = hi[k].max(p[k]);
                }
            }
            let mut ext = [0.0; 3];
            for k in 0..3 {
                // Pad by one radius so a degenerate axis still has positive length.
                ext[k] = (hi[k] - lo[k]).max(radius);
            }
            (lo, ext)
        };

        let mut dims = [1usize; 3];
        let mut cell_len = [0.0; 3];
        for k in 0..3 {
            // Floor keeps every cell at least `radius` wide.
            dims[k] = ((extent[k] / radius).floor() as usize).max(1);
            cell_len[k] = extent[k] / dims[k] as f64;
        }
        let ncells = dims[0] * dims[1] * dims[2];

        let cell_of = |a: usize| -> usize {
            let mut idx = [0usize; 3];
            for k in 0..3 {
                let p = [pos.x[a], pos.y[a], pos.z[a]][k];
                let rel = match box_len {
                    Some(l) => (p - origin[k]).rem_euclid(l[k]),
                    None => p - origin[k],
                };
                idx[k] = ((rel / cell_len[k]) as usize).min(dims[k] - 1);
            }
            (idx[2] * dims[1] + idx[1]) * dims[0] + idx[0]
        };

        let cells: Vec<usize> = (0..n).map(cell_of).collect();
        let mut start = vec![0usize; ncells + 1];
        for &c in &cells {
            start[c + 1] += 1;
        }
        for c in 0..ncells {
            start[c + 1] += start[c];
        }
        let mut fill = start.clone();
        let mut orig = vec![0u32; n];
        let mut x = vec![0.0; n];
        let mut y = vec![0.0; n];
        let mut z = vec![0.0; n];
        for (a, &c) in cells.iter().enumerate() {
            let s = fill[c];
            fill[c] += 1;
            orig[s] = a as u32;
            x[s] = pos.x[a];
            y[s] = pos.y[a];
            z[s] = pos.z[a];
        }

        Self { dims, periodic: box_len, start, orig, x, y, z }
    }

    fn ncells(&self) -> usize {
        self.dims[0] * self.dims[1] * self.dims[2]
    }

    fn range(&self, c: usize) -> (usize, usize) {
        (self.start[c], self.start[c + 1])
    }

    /// Unique, ascending indices of the cells adjacent to `c`, including `c`.
    fn neighbour_cells(&self, c: usize) -> Vec<usize> {
        let [nx, ny, nz] = self.dims;
        let cx = c % nx;
        let cy = (c / nx) % ny;
        let cz = c / (nx * ny);
        let mut out = Vec::with_capacity(27);
        for dz in -1i64..=1 {
            for dy in -1i64..=1 {
                for dx in -1i64..=1 {
                    let mut idx = [0usize; 3];
                    let mut ok = true;
                    for (k, (base, off, m)) in
                        [(cx, dx, nx), (cy, dy, ny), (cz, dz, nz)].into_iter().enumerate()
                    {
                        let v = base as i64 + off;
                        let m = m as i64;
                        if self.periodic.is_some() {
                            idx[k] = v.rem_euclid(m) as usize;
                        } else if (0..m).contains(&v) {
                            idx[k] = v as usize;
                        } else {
                            ok = false;
                        }
                    }
                    if ok {
                        out.push((idx[2] * ny + idx[1]) * nx + idx[0]);
                    }
                }
            }
        }
        out.sort_unstable();
        out.dedup();
        out
    }

    /// Compares sorted atom `i` with sorted slots `js..je` four lanes at a time.
    fn scan(
        &self,
        i: usize,
        js: usize,
        je: usize,
        rc2: F64x4,
        box_len: Option<[f64; 3]>,
        out: &mut Vec<(u32, u32)>,
    ) {
        let xi = F64x4::splat(self.x[i]);
        let yi = F64x4::splat(self.y[i]);
        let zi = F64x4::splat(self.z[i]);
        let mut j = js;
        while j < je {
            let take = (je - j).min(LANES);
            let lane_of = |src: &[f64]| -> F64x4 {
                F64x4::from(std::array::from_fn(|l| if l < take { src[j + l] } else { 0.0 }))
            };
            let mut dx = xi - lane_of(&self.x);
            let mut dy = yi - lane_of(&self.y);
            let mut dz = zi - lane_of(&self.z);
            if let Some(l) = box_len {
                dx = min_image_v(dx, l[0]);
                dy = min_image_v(dy, l[1]);
                dz = min_image_v(dz, l[2]);
            }
            let r2 = dx * dx + dy * dy + dz * dz;
            let valid = SimdMask::<f64, LANES>::from_array(std::array::from_fn(|l| l < take));
            let mask = r2.simd_lt(rc2) & valid;
            let bits = mask.to_bitmask();
            for l in 0..take {
                if bits & (1u64 << l) != 0 {
                    let a = self.orig[i];
                    let b = self.orig[j + l];
                    out.push((a.min(b), a.max(b)));
                }
            }
            j += take;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::min_image;

    /// Small deterministic generator (xorshift64*), avoids a dev-dependency.
    struct Rng(u64);
    impl Rng {
        fn next_f64(&mut self) -> f64 {
            self.0 ^= self.0 >> 12;
            self.0 ^= self.0 << 25;
            self.0 ^= self.0 >> 27;
            let v = self.0.wrapping_mul(0x2545_F491_4F6C_DD1D);
            (v >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    fn random_positions(n: usize, extent: [f64; 3], seed: u64) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let mut rng = Rng(seed);
        let mut gen = |e: f64| (0..n).map(|_| rng.next_f64() * e).collect::<Vec<f64>>();
        (gen(extent[0]), gen(extent[1]), gen(extent[2]))
    }

    /// O(N^2) reference with the same arithmetic as the kernel.
    fn brute_force(pos: Positions<'_>, radius: f64, box_len: Option<[f64; 3]>) -> Vec<(u32, u32)> {
        let n = pos.len();
        let mut out = Vec::new();
        for i in 0..n {
            for j in i + 1..n {
                let mut dx = pos.x[i] - pos.x[j];
                let mut dy = pos.y[i] - pos.y[j];
                let mut dz = pos.z[i] - pos.z[j];
                if let Some(l) = box_len {
                    dx = min_image(dx, l[0]);
                    dy = min_image(dy, l[1]);
                    dz = min_image(dz, l[2]);
                }
                if dx * dx + dy * dy + dz * dz < radius * radius {
                    out.push((i as u32, j as u32));
                }
            }
        }
        out
    }

    fn check(n: usize, extent: [f64; 3], radius: f64, periodic: bool, seed: u64) {
        let (x, y, z) = random_positions(n, extent, seed);
        let pos = Positions::new(&x, &y, &z);
        let box_len = periodic.then_some(extent);
        let got = build_pair_list(pos, radius, box_len);
        let want = brute_force(pos, radius, box_len);
        assert_eq!(got, want, "n={n} extent={extent:?} r={radius} pbc={periodic}");
    }

    #[test]
    fn open_boundaries_match_brute_force() {
        for (seed, n) in [(1, 2), (2, 17), (3, 200), (4, 513)] {
            check(n, [10.0, 7.0, 4.0], 1.3, false, seed);
        }
    }

    #[test]
    fn periodic_cubic_matches_brute_force() {
        for (seed, n) in [(11, 5), (12, 64), (13, 400)] {
            check(n, [12.0, 12.0, 12.0], 2.5, true, seed);
        }
    }

    #[test]
    fn periodic_triclinic_axes_match_brute_force() {
        // Non-cubic box with uneven cell counts per axis.
        check(300, [9.0, 14.0, 21.0], 3.1, true, 21);
    }

    #[test]
    fn tiny_periodic_boxes_with_one_or_two_cells_match() {
        // Radius equal to half the edge gives 2 cells per axis (wrap dedup).
        check(40, [6.0, 6.0, 6.0], 3.0, true, 31);
        // A single cell per axis: radius just over a third of the edge is not
        // enough for 2 cells, so the wrap-around cell is the cell itself.
        check(30, [6.0, 6.0, 6.0], 2.9, true, 33);
        // Radius equal to half the edge (the largest allowed) also gives 2 cells.
        check(25, [4.0, 4.0, 4.0], 2.0, true, 32);
    }

    #[test]
    fn coincident_atoms_are_paired() {
        let x = [1.0, 1.0, 5.0];
        let y = [2.0, 2.0, 2.0];
        let z = [3.0, 3.0, 3.0];
        let pos = Positions::new(&x, &y, &z);
        assert_eq!(build_pair_list(pos, 0.5, None), vec![(0, 1)]);
    }

    #[test]
    fn degenerate_inputs_give_empty_lists() {
        let empty = Positions::new(&[], &[], &[]);
        assert_eq!(build_pair_list(empty, 1.0, None), Vec::<(u32, u32)>::new());
        let one = Positions::new(&[0.0], &[0.0], &[0.0]);
        assert_eq!(build_pair_list(one, 1.0, None), Vec::<(u32, u32)>::new());
        let two = Positions::new(&[0.0, 0.5], &[0.0, 0.0], &[0.0, 0.0]);
        assert_eq!(build_pair_list(two, 0.0, None), Vec::<(u32, u32)>::new());
    }

    #[test]
    fn verlet_list_is_superset_of_cutoff_pairs() {
        let extent = [15.0, 15.0, 15.0];
        let (x, y, z) = random_positions(300, extent, 41);
        let pos = Positions::new(&x, &y, &z);
        let (rc, skin) = (2.0, 0.4);
        let verlet = build_verlet_list(pos, rc, skin, Some(extent));
        let cut = brute_force(pos, rc, Some(extent));
        assert!(cut.iter().all(|p| verlet.binary_search(p).is_ok()));
        assert_eq!(verlet, brute_force(pos, rc + skin, Some(extent)));
    }

    #[test]
    #[should_panic(expected = "half the smallest box edge")]
    fn radius_beyond_half_box_is_rejected() {
        let x = [0.0, 1.0];
        let pos = Positions::new(&x, &x, &x);
        let _ = build_pair_list(pos, 3.0, Some([4.0, 4.0, 4.0]));
    }
}
