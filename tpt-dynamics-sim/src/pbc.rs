//! Periodic boundary conditions and the minimum-image convention for
//! orthorhombic and triclinic boxes.
//!
//! The cell is stored as a lower-triangular matrix in GROMACS convention:
//!
//! ```text
//! a = (lx, 0, 0)
//! b = (xy, ly, 0)
//! c = (xz, yz, lz)
//! ```
//!
//! [`SimBox::lengths`] holds the diagonal `(lx, ly, lz)` and [`SimBox::tilt`]
//! holds `[xy, xz, yz]`. A box with zero tilt takes the orthorhombic fast path.
//! Triclinic boxes are restricted to the reduced form `|xy| <= lx/2`,
//! `|xz| <= lx/2`, `|yz| <= ly/2`, which makes the 27-image search in
//! [`SimBox::min_image`] exact.

use tpt_dynamics_core::Vec3;

use crate::error::SimError;

/// A periodic box with edge lengths in nm, optionally triclinic.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SimBox {
    /// Diagonal edge lengths `(lx, ly, lz)` in nm. All strictly positive.
    pub lengths: Vec3,
    /// Off-diagonal tilt factors `[xy, xz, yz]` in nm (all zero when orthorhombic).
    pub tilt: [f64; 3],
}

impl SimBox {
    /// Creates an orthorhombic box, rejecting non-finite or non-positive edge lengths.
    ///
    /// # Errors
    /// Returns [`SimError::Invalid`] for an invalid edge length.
    pub fn new(lx: f64, ly: f64, lz: f64) -> Result<Self, SimError> {
        Self::triclinic(lx, ly, lz, 0.0, 0.0, 0.0)
    }

    /// Creates a (possibly triclinic) box in GROMACS lower-triangular form.
    ///
    /// # Errors
    /// Returns [`SimError::Invalid`] for invalid edge lengths, non-finite tilts,
    /// or tilts outside the reduced range `|xy| <= lx/2`, `|xz| <= lx/2`,
    /// `|yz| <= ly/2`.
    pub fn triclinic(
        lx: f64,
        ly: f64,
        lz: f64,
        xy: f64,
        xz: f64,
        yz: f64,
    ) -> Result<Self, SimError> {
        for (name, l) in [("x", lx), ("y", ly), ("z", lz)] {
            if !(l.is_finite() && l > 0.0) {
                return Err(SimError::Invalid(format!(
                    "box length along {name} must be finite and > 0, got {l}"
                )));
            }
        }
        if ![xy, xz, yz].iter().all(|t| t.is_finite()) {
            return Err(SimError::Invalid("box tilt must be finite".into()));
        }
        if xy.abs() > 0.5 * lx || xz.abs() > 0.5 * lx || yz.abs() > 0.5 * ly {
            return Err(SimError::Invalid(
                "triclinic tilts must satisfy |xy|,|xz| <= lx/2 and |yz| <= ly/2".into(),
            ));
        }
        Ok(Self { lengths: Vec3::new(lx, ly, lz), tilt: [xy, xz, yz] })
    }

    /// True if the box has no tilt.
    pub fn is_orthorhombic(&self) -> bool {
        self.tilt == [0.0; 3]
    }

    /// Box volume in nm³.
    pub fn volume(&self) -> f64 {
        let l = self.lengths;
        l.x * l.y * l.z
    }

    /// Smallest diagonal edge length in nm.
    pub fn min_edge(&self) -> f64 {
        self.lengths.x.min(self.lengths.y).min(self.lengths.z)
    }

    /// Smallest perpendicular distance between opposite faces (nm). Equal to
    /// [`SimBox::min_edge`] for orthorhombic boxes. Neighbour radii must stay
    /// below half of this value.
    pub fn min_width(&self) -> f64 {
        if self.is_orthorhombic() {
            return self.min_edge();
        }
        let [a, b, c] = self.vectors();
        let v = self.volume();
        let w = |x: Vec3, y: Vec3| v / cross(x, y).norm();
        w(b, c).min(w(c, a)).min(w(a, b))
    }

    /// The three cell vectors `a`, `b`, `c`.
    pub fn vectors(&self) -> [Vec3; 3] {
        let l = self.lengths;
        let [xy, xz, yz] = self.tilt;
        [Vec3::new(l.x, 0.0, 0.0), Vec3::new(xy, l.y, 0.0), Vec3::new(xz, yz, l.z)]
    }

    /// Converts lattice coordinates `s` (multiples of the cell vectors) to Cartesian.
    pub fn cartesian(&self, s: Vec3) -> Vec3 {
        let [a, b, c] = self.vectors();
        a * s.x + b * s.y + c * s.z
    }

    /// Converts a Cartesian vector to lattice coordinates.
    pub fn fractional(&self, r: Vec3) -> Vec3 {
        let l = self.lengths;
        let [xy, xz, yz] = self.tilt;
        let s2 = r.z / l.z;
        let s1 = (r.y - s2 * yz) / l.y;
        let s0 = (r.x - s1 * xy - s2 * xz) / l.x;
        Vec3::new(s0, s1, s2)
    }

    /// Applies the minimum-image convention to a displacement vector.
    ///
    /// Orthorhombic boxes shift each component by a whole number of box lengths.
    /// Triclinic boxes take the shortest of the 27 lattice translations around
    /// the rounded lattice coordinates.
    pub fn min_image(&self, dr: Vec3) -> Vec3 {
        if self.is_orthorhombic() {
            let l = self.lengths;
            return Vec3::new(
                dr.x - l.x * (dr.x / l.x).round(),
                dr.y - l.y * (dr.y / l.y).round(),
                dr.z - l.z * (dr.z / l.z).round(),
            );
        }
        let s = self.fractional(dr);
        let n = Vec3::new(s.x.round(), s.y.round(), s.z.round());
        let mut best = dr - self.cartesian(n);
        let mut best2 = best.norm_sq();
        for i in -1..=1 {
            for j in -1..=1 {
                for k in -1..=1 {
                    let shift = n + Vec3::new(f64::from(i), f64::from(j), f64::from(k));
                    let cand = dr - self.cartesian(shift);
                    let c2 = cand.norm_sq();
                    if c2 < best2 {
                        best2 = c2;
                        best = cand;
                    }
                }
            }
        }
        best
    }

    /// Wraps a position into the primary cell: `[0, L)` per axis for orthorhombic
    /// boxes, lattice coordinates in `[0, 1)` for triclinic boxes.
    pub fn wrap(&self, p: Vec3) -> Vec3 {
        if self.is_orthorhombic() {
            let l = self.lengths;
            return Vec3::new(wrap1(p.x, l.x), wrap1(p.y, l.y), wrap1(p.z, l.z));
        }
        let s = self.fractional(p);
        self.cartesian(Vec3::new(wrap1(s.x, 1.0), wrap1(s.y, 1.0), wrap1(s.z, 1.0)))
    }
}

fn cross(a: Vec3, b: Vec3) -> Vec3 {
    Vec3::new(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x)
}

fn wrap1(x: f64, l: f64) -> f64 {
    let w = x - l * (x / l).floor();
    // Guard against round-off returning exactly `l`.
    if w >= l {
        0.0
    } else {
        w
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rng::Rng;

    #[test]
    fn rejects_bad_lengths() {
        assert!(SimBox::new(0.0, 1.0, 1.0).is_err());
        assert!(SimBox::new(1.0, f64::NAN, 1.0).is_err());
        assert!(SimBox::new(1.0, 1.0, -2.0).is_err());
        assert!(SimBox::triclinic(2.0, 2.0, 2.0, 1.5, 0.0, 0.0).is_err());
        assert!(SimBox::triclinic(2.0, 2.0, 2.0, 0.0, 0.0, 1.5).is_err());
    }

    #[test]
    fn min_image_picks_nearest_copy() {
        let b = SimBox::new(10.0, 10.0, 10.0).unwrap();
        let d = b.min_image(Vec3::new(9.0, -6.0, 4.0));
        assert!((d.x + 1.0).abs() < 1e-12);
        assert!((d.y - 4.0).abs() < 1e-12);
        assert!((d.z - 4.0).abs() < 1e-12);
    }

    #[test]
    fn wrap_into_primary_box() {
        let b = SimBox::new(2.0, 3.0, 4.0).unwrap();
        let w = b.wrap(Vec3::new(-0.5, 7.0, 4.0));
        assert!((w.x - 1.5).abs() < 1e-12);
        assert!((w.y - 1.0).abs() < 1e-12);
        assert!(w.z.abs() < 1e-12);
        let tiny = b.wrap(Vec3::new(-1e-300, 0.0, 0.0));
        assert!(tiny.x < 2.0);
    }

    /// Reference search over every lattice translation whose index is within
    /// `R` cells of the origin, with `R` large enough to cover the displacement.
    fn brute_min_image(b: &SimBox, dr: Vec3) -> Vec3 {
        let reach = (dr.norm() / b.min_width()).ceil() as i64 + 2;
        let mut best = dr;
        let mut best2 = f64::MAX;
        for i in -reach..=reach {
            for j in -reach..=reach {
                for k in -reach..=reach {
                    let shift = Vec3::new(i as f64, j as f64, k as f64);
                    let c = dr - b.cartesian(shift);
                    if c.norm_sq() < best2 {
                        best2 = c.norm_sq();
                        best = c;
                    }
                }
            }
        }
        best
    }

    #[test]
    fn triclinic_min_image_matches_brute_force() {
        let mut rng = Rng::seed_from_u64(2024);
        for _ in 0..20 {
            let lx = 2.0 + 2.0 * rng.uniform();
            let ly = 2.0 + 2.0 * rng.uniform();
            let lz = 2.0 + 2.0 * rng.uniform();
            let xy = (rng.uniform() - 0.5) * lx;
            let xz = (rng.uniform() - 0.5) * lx;
            let yz = (rng.uniform() - 0.5) * ly;
            let b = SimBox::triclinic(lx, ly, lz, xy, xz, yz).unwrap();
            for _ in 0..200 {
                let dr = Vec3::new(rng.gaussian(), rng.gaussian(), rng.gaussian()) * 4.0;
                let fast = b.min_image(dr);
                let slow = brute_min_image(&b, dr);
                assert!(
                    (fast.norm() - slow.norm()).abs() < 1e-9,
                    "fast {} slow {}",
                    fast.norm(),
                    slow.norm()
                );
            }
        }
    }

    #[test]
    fn triclinic_wrap_and_fractional_round_trip() {
        let b = SimBox::triclinic(3.0, 2.5, 4.0, 1.0, -0.7, 0.4).unwrap();
        let p = Vec3::new(5.1, -2.2, 9.7);
        let back = b.cartesian(b.fractional(p));
        assert!((back - p).norm() < 1e-12);
        let sw = b.fractional(b.wrap(p));
        for c in [sw.x, sw.y, sw.z] {
            assert!((-1e-12..1.0 + 1e-12).contains(&c), "lattice coordinate {c} not wrapped");
        }
        assert!(b.min_width() > 0.0 && b.min_width() <= b.min_edge() + 1e-12);
    }
}
