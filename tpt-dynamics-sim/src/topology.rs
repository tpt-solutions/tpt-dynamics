//! Molecular topology: bonded term index lists and the pair lists that are
//! derived from the bond graph.
//!
//! The topology stores only atom indices. Parameters are attached later by the
//! force-field parameter set (see [`crate::ffparams`]). Angles, proper dihedrals
//! and 1-4 pairs are generated automatically from the bonds; impropers are
//! listed explicitly because their central atom is chemistry-specific.
//!
//! Exclusion convention (as in GROMACS with `nrexcl = 3`): atoms whose shortest
//! bond-path separation is 1, 2 or 3 bonds are excluded from the ordinary
//! non-bonded sum. Pairs separated by exactly three bonds (1-4 pairs) get an
//! explicit scaled interaction instead. Electrostatic exclusions are still
//! corrected for in reciprocal-space (Ewald) runs.

use std::collections::{BTreeSet, VecDeque};

use crate::error::SimError;

/// Bonded and pair index lists for a system.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Topology {
    /// Covalent bonds `(i, j)`.
    pub bonds: Vec<(usize, usize)>,
    /// Bond angles `[i, j, k]` with `j` the central atom.
    pub angles: Vec<[usize; 3]>,
    /// Proper dihedrals `[i, j, k, l]` along the chain `i-j-k-l`.
    pub dihedrals: Vec<[usize; 4]>,
    /// Improper dihedrals `[i, j, k, l]`, with `k` the central atom.
    pub impropers: Vec<[usize; 4]>,
    /// 1-4 pairs `(i, l)`: atoms separated by exactly three bonds.
    pub pairs14: Vec<(usize, usize)>,
    /// Excluded pairs `(i, j)`, `i < j`, separated by at most three bonds.
    pub exclusions: Vec<(usize, usize)>,
}

impl Topology {
    /// Builds the angle, dihedral, 1-4 and exclusion lists from a bond list for
    /// `n_atoms` atoms. `impropers` starts empty.
    ///
    /// # Errors
    /// Returns [`SimError::Invalid`] if a bond refers to an atom out of range or
    /// connects an atom to itself.
    pub fn from_bonds(n_atoms: usize, bonds: &[(usize, usize)]) -> Result<Self, SimError> {
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n_atoms];
        for &(i, j) in bonds {
            if i >= n_atoms || j >= n_atoms {
                return Err(SimError::Invalid(format!("bond ({i}, {j}) out of range")));
            }
            if i == j {
                return Err(SimError::Invalid(format!(
                    "bond ({i}, {i}) connects an atom to itself"
                )));
            }
            adj[i].push(j);
            adj[j].push(i);
        }
        for a in &mut adj {
            a.sort_unstable();
            a.dedup();
        }
        let mut topo = Self { bonds: bonds.to_vec(), ..Self::default() };

        // Angles: every pair of neighbours of each central atom.
        for (j, nb) in adj.iter().enumerate() {
            for (a, &i) in nb.iter().enumerate() {
                for &k in &nb[a + 1..] {
                    topo.angles.push([i, j, k]);
                }
            }
        }

        // Proper dihedrals: every central bond (j, k) with j < k, then every
        // outer pair i (neighbour of j) and l (neighbour of k) with i != l.
        let mut seen: BTreeSet<[usize; 4]> = BTreeSet::new();
        for &(a, b) in bonds {
            let (j, k) = if a < b { (a, b) } else { (b, a) };
            for &i in &adj[j] {
                if i == k {
                    continue;
                }
                for &l in &adj[k] {
                    if l == j || l == i {
                        continue;
                    }
                    let path = [i, j, k, l];
                    let rev = [l, k, j, i];
                    let key = if path <= rev { path } else { rev };
                    if seen.insert(key) {
                        topo.dihedrals.push(key);
                    }
                }
            }
        }

        // Exclusions and 1-4 pairs from shortest bond-path distances up to 3.
        for i in 0..n_atoms {
            let dist = bounded_bfs(&adj, i, 3);
            for (j, &d) in dist.iter().enumerate() {
                if j <= i || d == usize::MAX {
                    continue;
                }
                topo.exclusions.push((i, j));
                if d == 3 {
                    topo.pairs14.push((i, j));
                }
            }
        }
        Ok(topo)
    }

    /// Appends the terms of `other` with every atom index shifted by `offset`.
    /// Used when several molecules are placed into one system.
    pub fn append_shifted(&mut self, other: &Self, offset: usize) {
        let s2 = |p: (usize, usize)| (p.0 + offset, p.1 + offset);
        self.bonds.extend(other.bonds.iter().map(|&p| s2(p)));
        self.angles.extend(other.angles.iter().map(|a| a.map(|x| x + offset)));
        self.dihedrals.extend(other.dihedrals.iter().map(|a| a.map(|x| x + offset)));
        self.impropers.extend(other.impropers.iter().map(|a| a.map(|x| x + offset)));
        self.pairs14.extend(other.pairs14.iter().map(|&p| s2(p)));
        self.exclusions.extend(other.exclusions.iter().map(|&p| s2(p)));
    }

    /// Sorted list of excluded partners for every atom, both directions included.
    pub fn exclusion_partners(&self, n_atoms: usize) -> Vec<Vec<usize>> {
        let mut out = vec![Vec::new(); n_atoms];
        for &(i, j) in &self.exclusions {
            out[i].push(j);
            out[j].push(i);
        }
        for v in &mut out {
            v.sort_unstable();
            v.dedup();
        }
        out
    }
}

/// Breadth-first distances from `start` in bond steps, not expanded beyond `max`.
/// Unreached atoms get `usize::MAX`.
fn bounded_bfs(adj: &[Vec<usize>], start: usize, max: usize) -> Vec<usize> {
    let mut dist = vec![usize::MAX; adj.len()];
    dist[start] = 0;
    let mut queue = VecDeque::from([start]);
    while let Some(u) = queue.pop_front() {
        if dist[u] >= max {
            continue;
        }
        for &v in &adj[u] {
            if dist[v] == usize::MAX {
                dist[v] = dist[u] + 1;
                queue.push_back(v);
            }
        }
    }
    dist
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_chain_of_five_atoms() {
        // 0-1-2-3-4
        let t = Topology::from_bonds(5, &[(0, 1), (1, 2), (2, 3), (3, 4)]).unwrap();
        assert_eq!(t.angles.len(), 3);
        assert_eq!(t.dihedrals.len(), 2);
        assert_eq!(t.pairs14, vec![(0, 3), (1, 4)]);
        // Excluded: all pairs within three bonds: total pairs 10 minus the (0,4) pair.
        assert_eq!(t.exclusions.len(), 9);
        assert!(!t.exclusions.contains(&(0, 4)));
    }

    #[test]
    fn ring_has_no_duplicate_dihedrals() {
        // Four-membered ring 0-1-2-3-0.
        let t = Topology::from_bonds(4, &[(0, 1), (1, 2), (2, 3), (3, 0)]).unwrap();
        assert_eq!(t.angles.len(), 4);
        // Each of the 4 central bonds gives a path; i != l rules out degenerate ones.
        let mut sorted = t.dihedrals.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), t.dihedrals.len());
        // Every pair is within three bonds, so 1-4 pairs are the diagonal pairs at
        // distance 2 only, i.e. none at distance exactly 3.
        assert_eq!(t.pairs14, Vec::<(usize, usize)>::new());
        assert_eq!(t.exclusions.len(), 6);
    }

    #[test]
    fn rejects_bad_bonds() {
        assert!(Topology::from_bonds(2, &[(0, 5)]).is_err());
        assert!(Topology::from_bonds(2, &[(1, 1)]).is_err());
    }

    #[test]
    fn append_shifts_indices() {
        let a = Topology::from_bonds(2, &[(0, 1)]).unwrap();
        let mut t = Topology::default();
        t.append_shifted(&a, 0);
        t.append_shifted(&a, 2);
        assert_eq!(t.bonds, vec![(0, 1), (2, 3)]);
        assert_eq!(t.exclusion_partners(4)[3], vec![2]);
    }
}
