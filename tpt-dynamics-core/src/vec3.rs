//! A minimal 3-component `f64` vector for physics code.
//!
//! This is a deliberately small, `no_std`-friendly type. The workspace will
//! migrate to `glam` or `tpt-math` types behind this interface if benchmarks
//! or APIs call for it (see `docs/dependencies.md`).

use core::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Sub, SubAssign};

/// A 3-dimensional vector of `f64` components.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
#[repr(C)]
pub struct Vec3 {
    /// x component.
    pub x: f64,
    /// y component.
    pub y: f64,
    /// z component.
    pub z: f64,
}

impl Vec3 {
    /// The zero vector.
    pub const ZERO: Self = Self { x: 0.0, y: 0.0, z: 0.0 };

    /// Construct from components.
    #[inline]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Dot product.
    #[inline]
    pub fn dot(self, rhs: Self) -> f64 {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }

    /// Cross product.
    #[inline]
    #[must_use]
    pub fn cross(self, rhs: Self) -> Self {
        Self {
            x: self.y * rhs.z - self.z * rhs.y,
            y: self.z * rhs.x - self.x * rhs.z,
            z: self.x * rhs.y - self.y * rhs.x,
        }
    }

    /// Squared Euclidean norm.
    #[inline]
    pub fn norm_sq(self) -> f64 {
        self.dot(self)
    }

    /// Euclidean norm.
    #[inline]
    pub fn norm(self) -> f64 {
        libm::sqrt(self.norm_sq())
    }

    /// Unit vector in the same direction, or `None` if the norm is zero.
    #[inline]
    pub fn normalized(self) -> Option<Self> {
        let n = self.norm();
        if n > 0.0 {
            Some(self / n)
        } else {
            None
        }
    }

    /// Component-wise multiplication (Hadamard product).
    #[inline]
    #[must_use]
    pub fn hadamard(self, rhs: Self) -> Self {
        Self { x: self.x * rhs.x, y: self.y * rhs.y, z: self.z * rhs.z }
    }
}

impl Add for Vec3 {
    type Output = Self;
    #[inline]
    fn add(self, rhs: Self) -> Self {
        Self { x: self.x + rhs.x, y: self.y + rhs.y, z: self.z + rhs.z }
    }
}

impl Sub for Vec3 {
    type Output = Self;
    #[inline]
    fn sub(self, rhs: Self) -> Self {
        Self { x: self.x - rhs.x, y: self.y - rhs.y, z: self.z - rhs.z }
    }
}

impl Neg for Vec3 {
    type Output = Self;
    #[inline]
    fn neg(self) -> Self {
        Self { x: -self.x, y: -self.y, z: -self.z }
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;
    #[inline]
    fn mul(self, s: f64) -> Self {
        Self { x: self.x * s, y: self.y * s, z: self.z * s }
    }
}

impl Mul<Vec3> for f64 {
    type Output = Vec3;
    #[inline]
    fn mul(self, v: Vec3) -> Vec3 {
        v * self
    }
}

impl Div<f64> for Vec3 {
    type Output = Self;
    #[inline]
    fn div(self, s: f64) -> Self {
        Self { x: self.x / s, y: self.y / s, z: self.z / s }
    }
}

impl AddAssign for Vec3 {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}

impl SubAssign for Vec3 {
    #[inline]
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl MulAssign<f64> for Vec3 {
    #[inline]
    fn mul_assign(&mut self, s: f64) {
        *self = *self * s;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arithmetic() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(-1.0, 0.5, 2.0);
        assert_eq!(a + b, Vec3::new(0.0, 2.5, 5.0));
        assert_eq!(a - b, Vec3::new(2.0, 1.5, 1.0));
        assert_eq!(-a, Vec3::new(-1.0, -2.0, -3.0));
        assert_eq!(2.0 * a, a * 2.0);
        assert_eq!(a / 2.0, Vec3::new(0.5, 1.0, 1.5));
    }

    #[test]
    #[allow(clippy::many_single_char_names)] // x, y, a, b, c are standard vector notation
    fn dot_and_cross() {
        let x = Vec3::new(1.0, 0.0, 0.0);
        let y = Vec3::new(0.0, 1.0, 0.0);
        assert_eq!(x.dot(y), 0.0);
        assert_eq!(x.cross(y), Vec3::new(0.0, 0.0, 1.0));
        // a × b is orthogonal to both inputs.
        let a = Vec3::new(0.3, -1.2, 2.5);
        let b = Vec3::new(4.0, 0.1, -0.7);
        let c = a.cross(b);
        assert!(c.dot(a).abs() < 1e-12);
        assert!(c.dot(b).abs() < 1e-12);
    }

    #[test]
    fn norm_and_normalize() {
        let v = Vec3::new(3.0, 4.0, 0.0);
        assert_eq!(v.norm(), 5.0);
        assert_eq!(v.norm_sq(), 25.0);
        let u = v.normalized().unwrap();
        assert!((u.norm() - 1.0).abs() < 1e-15);
        assert_eq!(Vec3::ZERO.normalized(), None);
    }

    #[test]
    fn assign_ops() {
        let mut v = Vec3::new(1.0, 1.0, 1.0);
        v += Vec3::new(1.0, 2.0, 3.0);
        v -= Vec3::new(0.0, 0.0, 1.0);
        v *= 2.0;
        assert_eq!(v, Vec3::new(4.0, 6.0, 6.0));
    }
}
