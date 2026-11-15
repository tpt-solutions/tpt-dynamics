//! Neighbour search: cell lists, brute-force reference, and Verlet lists with a
//! skin and a displacement-based rebuild criterion.
//!
//! A [`VerletList`] stores every pair `(i, j)` with `i < j` whose minimum-image
//! distance at build time is below `cutoff + skin`. The list stays valid while
//! no atom has moved more than `skin / 2` since the build: then any pair closer
//! than `cutoff` now was within `cutoff + skin` at build time.
//!
//! The list build is parallel (rayon) over atoms. Cells are at least as large as
//! the list radius, so only the 27 neighbouring cells (deduplicated for small
//! boxes) need to be searched.

use rayon::prelude::*;
use tpt_dynamics_core::Vec3;

use crate::error::SimError;
use crate::pbc::SimBox;

/// Brute-force O(N²) list of pairs `(i, j)`, `i < j`, with minimum-image
/// distance strictly below `radius`. Reference implementation for tests.
pub fn brute_force_pairs(positions: &[Vec3], sim_box: &SimBox, radius: f64) -> Vec<(usize, usize)> {
    let r2 = radius * radius;
    let mut out = Vec::new();
    for i in 0..positions.len() {
        for j in (i + 1)..positions.len() {
            let dr = sim_box.min_image(positions[i] - positions[j]);
            if dr.norm_sq() < r2 {
                out.push((i, j));
            }
        }
    }
    out
}

/// A uniform grid of cells over the periodic box.
#[derive(Clone, Debug)]
pub struct CellList {
    dims: [usize; 3],
    cell_of: Vec<usize>,
    /// CSR offsets into `members`, length `ncells + 1`.
    start: Vec<usize>,
    members: Vec<usize>,
}

impl CellList {
    /// Builds a cell list whose cells are at least `min_cell` wide.
    ///
    /// # Errors
    /// Returns [`SimError::Invalid`] if `min_cell` is not finite and positive.
    pub fn build(positions: &[Vec3], sim_box: &SimBox, min_cell: f64) -> Result<Self, SimError> {
        if !(min_cell.is_finite() && min_cell > 0.0) {
            return Err(SimError::Invalid(format!("cell size must be > 0, got {min_cell}")));
        }
        let l = sim_box.lengths;
        let dims = [
            ((l.x / min_cell).floor() as usize).max(1),
            ((l.y / min_cell).floor() as usize).max(1),
            ((l.z / min_cell).floor() as usize).max(1),
        ];
        let ncells = dims[0] * dims[1] * dims[2];
        let cell_of: Vec<usize> = positions
            .iter()
            .map(|&p| {
                let w = sim_box.wrap(p);
                let ix = cell_index(w.x / l.x, dims[0]);
                let iy = cell_index(w.y / l.y, dims[1]);
                let iz = cell_index(w.z / l.z, dims[2]);
                (ix * dims[1] + iy) * dims[2] + iz
            })
            .collect();
        let mut counts = vec![0usize; ncells + 1];
        for &c in &cell_of {
            counts[c + 1] += 1;
        }
        for c in 0..ncells {
            counts[c + 1] += counts[c];
        }
        let start = counts.clone();
        let mut fill = counts;
        let mut members = vec![0usize; positions.len()];
        for (i, &c) in cell_of.iter().enumerate() {
            members[fill[c]] = i;
            fill[c] += 1;
        }
        Ok(Self { dims, cell_of, start, members })
    }

    /// Number of cells along each axis.
    pub fn dims(&self) -> [usize; 3] {
        self.dims
    }

    /// Atoms stored in cell `c`.
    fn members(&self, c: usize) -> &[usize] {
        &self.members[self.start[c]..self.start[c + 1]]
    }

    /// Distinct cell indices within one cell of cell `c` (including `c`).
    fn neighbour_cells(&self, c: usize) -> Vec<usize> {
        let [nx, ny, nz] = self.dims;
        let iz = c % nz;
        let iy = (c / nz) % ny;
        let ix = c / (nz * ny);
        let mut out = Vec::with_capacity(27);
        for dx in 0..3 {
            for dy in 0..3 {
                for dz in 0..3 {
                    let jx = wrap_index(ix, dx, nx);
                    let jy = wrap_index(iy, dy, ny);
                    let jz = wrap_index(iz, dz, nz);
                    let idx = (jx * ny + jy) * nz + jz;
                    if !out.contains(&idx) {
                        out.push(idx);
                    }
                }
            }
        }
        out
    }

    /// Cell index of atom `i`.
    pub fn cell_of(&self, i: usize) -> usize {
        self.cell_of[i]
    }

    /// Serial list of pairs `(i, j)`, `i < j`, with distance below `radius`.
    /// `radius` must not exceed the minimum cell size used at build time.
    pub fn pairs_within(
        &self,
        positions: &[Vec3],
        sim_box: &SimBox,
        radius: f64,
    ) -> Vec<(usize, usize)> {
        let r2 = radius * radius;
        let mut out = Vec::new();
        for (i, &p) in positions.iter().enumerate() {
            for c in self.neighbour_cells(self.cell_of[i]) {
                for &j in self.members(c) {
                    if j > i {
                        let dr = sim_box.min_image(p - positions[j]);
                        if dr.norm_sq() < r2 {
                            out.push((i, j));
                        }
                    }
                }
            }
        }
        out.sort_unstable();
        out
    }
}

fn cell_index(frac: f64, n: usize) -> usize {
    let idx = (frac * n as f64).floor();
    if idx <= 0.0 {
        0
    } else {
        (idx as usize).min(n - 1)
    }
}

fn wrap_index(i: usize, off: usize, n: usize) -> usize {
    // `off` is 0, 1 or 2, meaning a neighbour offset of -1, 0 or +1.
    (i + n + off - 1) % n
}

/// Verlet neighbour list with a skin and a displacement-based rebuild rule.
#[derive(Clone, Debug)]
pub struct VerletList {
    cutoff: f64,
    skin: f64,
    /// CSR offsets: neighbours of atom `i` are `neighbours[offsets[i]..offsets[i+1]]`
    /// and all of them are greater than `i`.
    offsets: Vec<usize>,
    neighbours: Vec<usize>,
    reference: Vec<Vec3>,
    builds: u64,
    /// Build per-atom lists with rayon (true) or serially (false).
    parallel: bool,
}

impl VerletList {
    /// Creates an empty list for interaction `cutoff` (nm) and buffer `skin` (nm).
    ///
    /// # Errors
    /// Returns [`SimError::Invalid`] if `cutoff` is not positive or `skin` is negative.
    pub fn new(cutoff: f64, skin: f64) -> Result<Self, SimError> {
        if !(cutoff.is_finite() && cutoff > 0.0) {
            return Err(SimError::Invalid(format!("cutoff must be > 0, got {cutoff}")));
        }
        if !(skin.is_finite() && skin >= 0.0) {
            return Err(SimError::Invalid(format!("skin must be >= 0, got {skin}")));
        }
        Ok(Self {
            cutoff,
            skin,
            offsets: Vec::new(),
            neighbours: Vec::new(),
            reference: Vec::new(),
            builds: 0,
            parallel: true,
        })
    }

    /// Selects a parallel (rayon) or serial list build. The serial build makes
    /// no rayon calls. Both produce the same list.
    pub fn set_parallel(&mut self, parallel: bool) {
        self.parallel = parallel;
    }

    /// True if the list build uses rayon.
    pub fn is_parallel(&self) -> bool {
        self.parallel
    }

    /// Interaction cutoff in nm.
    pub fn cutoff(&self) -> f64 {
        self.cutoff
    }

    /// Buffer width in nm.
    pub fn skin(&self) -> f64 {
        self.skin
    }

    /// Number of times the list has been (re)built.
    pub fn builds(&self) -> u64 {
        self.builds
    }

    /// Number of stored pairs.
    pub fn pair_count(&self) -> usize {
        self.neighbours.len()
    }

    /// Iterates over all stored pairs `(i, j)` with `i < j`.
    pub fn pairs(&self) -> impl Iterator<Item = (usize, usize)> + '_ {
        (0..self.offsets.len().saturating_sub(1)).flat_map(move |i| {
            self.neighbours[self.offsets[i]..self.offsets[i + 1]].iter().map(move |&j| (i, j))
        })
    }

    /// Forces a rebuild on the next [`VerletList::update`], for example after the
    /// box has changed size.
    pub fn invalidate(&mut self) {
        self.reference.clear();
        self.offsets.clear();
    }

    /// Returns true if the list must be rebuilt for the current positions:
    /// never built, atom count changed, or some atom moved more than `skin / 2`.
    pub fn needs_rebuild(&self, positions: &[Vec3], sim_box: &SimBox) -> bool {
        if self.reference.len() != positions.len() || self.offsets.len() != positions.len() + 1 {
            return true;
        }
        let limit2 = (0.5 * self.skin) * (0.5 * self.skin);
        positions
            .iter()
            .zip(&self.reference)
            .any(|(&p, &r)| sim_box.min_image(p - r).norm_sq() > limit2)
    }

    /// Rebuilds the list if required. Returns true if a rebuild happened.
    pub fn update(&mut self, positions: &[Vec3], sim_box: &SimBox) -> Result<bool, SimError> {
        if self.needs_rebuild(positions, sim_box) {
            self.rebuild(positions, sim_box)?;
            Ok(true)
        } else {
            Ok(false)
        }
    }

    /// Unconditionally rebuilds the list from the given positions (parallel over atoms).
    ///
    /// # Errors
    /// Returns [`SimError::Invalid`] if the list radius does not fit in the box
    /// (`cutoff + skin` must be less than half the smallest box edge so that the
    /// minimum-image convention is unambiguous).
    pub fn rebuild(&mut self, positions: &[Vec3], sim_box: &SimBox) -> Result<(), SimError> {
        let radius = self.cutoff + self.skin;
        if 2.0 * radius >= sim_box.min_width() {
            return Err(SimError::Invalid(format!(
                "cutoff + skin ({radius} nm) must be less than half the smallest box width ({} nm)",
                sim_box.min_width()
            )));
        }
        let r2 = radius * radius;
        let n = positions.len();
        let parallel = self.parallel;
        let per_atom: Vec<Vec<usize>> = if sim_box.is_orthorhombic() {
            let cells = CellList::build(positions, sim_box, radius)?;
            let one = |i: usize| -> Vec<usize> {
                let p = positions[i];
                let mut out = Vec::new();
                for c in cells.neighbour_cells(cells.cell_of(i)) {
                    for &j in cells.members(c) {
                        if j > i {
                            let dr = sim_box.min_image(p - positions[j]);
                            if dr.norm_sq() < r2 {
                                out.push(j);
                            }
                        }
                    }
                }
                out.sort_unstable();
                out
            };
            if parallel {
                (0..n).into_par_iter().map(one).collect()
            } else {
                (0..n).map(one).collect()
            }
        } else {
            // Cell grids are not defined for skewed cells, so triclinic boxes use the
            // O(N^2) search. This is documented in the crate docs.
            let one = |i: usize| -> Vec<usize> {
                (i + 1..n)
                    .filter(|&j| sim_box.min_image(positions[i] - positions[j]).norm_sq() < r2)
                    .collect()
            };
            if parallel {
                (0..n).into_par_iter().map(one).collect()
            } else {
                (0..n).map(one).collect()
            }
        };

        let total: usize = per_atom.iter().map(Vec::len).sum();
        let mut offsets = Vec::with_capacity(positions.len() + 1);
        let mut neighbours = Vec::with_capacity(total);
        offsets.push(0);
        for list in &per_atom {
            neighbours.extend_from_slice(list);
            offsets.push(neighbours.len());
        }
        self.offsets = offsets;
        self.neighbours = neighbours;
        self.reference = positions.to_vec();
        self.builds += 1;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Small deterministic LCG-based position generator (no RNG dependency needed).
    fn random_positions(n: usize, sim_box: &SimBox, seed: u64) -> Vec<Vec3> {
        let mut s =
            seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
        let mut next = move || {
            s = s.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1_442_695_040_888_963_407);
            (s >> 11) as f64 / (1u64 << 53) as f64
        };
        let l = sim_box.lengths;
        (0..n).map(|_| Vec3::new(next() * l.x, next() * l.y, next() * l.z)).collect()
    }

    fn sorted(mut v: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
        v.sort_unstable();
        v
    }

    #[test]
    fn cell_list_matches_brute_force_over_box_shapes() {
        let cases = [
            (SimBox::new(10.0, 10.0, 10.0).unwrap(), 300, 1.7, 2.0),
            (SimBox::new(6.0, 9.0, 7.5).unwrap(), 250, 1.2, 1.3),
            // Very small boxes: fewer than three cells along some axes.
            (SimBox::new(2.2, 2.5, 3.0).unwrap(), 60, 1.0, 1.0),
            (SimBox::new(4.0, 4.0, 4.0).unwrap(), 150, 1.9, 1.9),
        ];
        for (k, (sim_box, n, radius, min_cell)) in cases.iter().enumerate() {
            let pos = random_positions(*n, sim_box, k as u64 + 1);
            let cells = CellList::build(&pos, sim_box, *min_cell).unwrap();
            let got = cells.pairs_within(&pos, sim_box, *radius);
            let want = brute_force_pairs(&pos, sim_box, *radius);
            assert_eq!(got, sorted(want), "case {k}");
        }
    }

    #[test]
    fn verlet_list_matches_brute_force_within_cutoff() {
        let sim_box = SimBox::new(8.0, 8.5, 9.0).unwrap();
        let pos = random_positions(400, &sim_box, 42);
        let (cutoff, skin) = (1.5, 0.4);
        let mut vl = VerletList::new(cutoff, skin).unwrap();
        vl.rebuild(&pos, &sim_box).unwrap();
        let mut got: Vec<(usize, usize)> = vl
            .pairs()
            .filter(|&(i, j)| sim_box.min_image(pos[i] - pos[j]).norm() < cutoff)
            .collect();
        got.sort_unstable();
        let want = brute_force_pairs(&pos, &sim_box, cutoff);
        assert_eq!(got, want);
        // Every stored pair is inside cutoff + skin.
        for (i, j) in vl.pairs() {
            assert!(sim_box.min_image(pos[i] - pos[j]).norm() < cutoff + skin);
        }
    }

    #[test]
    fn rebuild_criterion_triggers_on_half_skin_displacement() {
        let sim_box = SimBox::new(10.0, 10.0, 10.0).unwrap();
        let mut pos = random_positions(50, &sim_box, 5);
        let mut vl = VerletList::new(1.0, 0.6).unwrap();
        assert!(vl.needs_rebuild(&pos, &sim_box));
        vl.rebuild(&pos, &sim_box).unwrap();
        assert!(!vl.needs_rebuild(&pos, &sim_box));
        pos[3].x += 0.29; // just under skin / 2
        assert!(!vl.needs_rebuild(&pos, &sim_box));
        pos[3].x += 0.02; // now over
        assert!(vl.needs_rebuild(&pos, &sim_box));
        assert!(vl.update(&pos, &sim_box).unwrap());
        assert_eq!(vl.builds(), 2);
    }

    #[test]
    fn rejects_list_radius_too_large_for_box() {
        let sim_box = SimBox::new(3.0, 3.0, 3.0).unwrap();
        let pos = random_positions(10, &sim_box, 1);
        let mut vl = VerletList::new(1.2, 0.4).unwrap();
        assert!(vl.rebuild(&pos, &sim_box).is_err());
        assert!(VerletList::new(-1.0, 0.1).is_err());
        assert!(VerletList::new(1.0, -0.1).is_err());
    }
}
