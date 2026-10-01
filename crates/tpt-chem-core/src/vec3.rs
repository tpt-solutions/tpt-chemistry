//! Fixed-size 3-vector of `f64`, the workhorse geometry type.
//!
//! Deliberately dependency-free so it works in `no_std` and can be stored
//! inline in SoA system buffers. Units are whatever the caller says — this
//! type is dimensionless.

use core::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

/// A 3-component double-precision vector.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
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
    pub const ZERO: Vec3 = Vec3 { x: 0.0, y: 0.0, z: 0.0 };

    /// Construct from components.
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Vec3 { x, y, z }
    }

    /// Component access by index (0, 1, 2). Panics out of range.
    pub fn get(&self, i: usize) -> f64 {
        match i {
            0 => self.x,
            1 => self.y,
            2 => self.z,
            _ => panic!("Vec3 index out of range"),
        }
    }

    /// Mutable component access by index (0, 1, 2). Panics out of range.
    pub fn set(&mut self, i: usize, v: f64) {
        match i {
            0 => self.x = v,
            1 => self.y = v,
            2 => self.z = v,
            _ => panic!("Vec3 index out of range"),
        }
    }

    /// Dot (inner) product.
    pub fn dot(self, o: Vec3) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// Cross product `self × o`.
    pub fn cross(self, o: Vec3) -> Vec3 {
        Vec3::new(
            self.y * o.z - self.z * o.y,
            self.z * o.x - self.x * o.z,
            self.x * o.y - self.y * o.x,
        )
    }

    /// Squared Euclidean norm (cheaper than [`Vec3::norm`]).
    pub fn norm_sq(self) -> f64 {
        self.dot(self)
    }

    /// Euclidean norm.
    pub fn norm(self) -> f64 {
        crate::num::sqrt(self.norm_sq())
    }

    /// Euclidean distance to `o`.
    pub fn dist(self, o: Vec3) -> f64 {
        (self - o).norm()
    }

    /// Squared Euclidean distance to `o`.
    pub fn dist_sq(self, o: Vec3) -> f64 {
        (self - o).norm_sq()
    }

    /// Unit vector in the same direction. The zero vector maps to itself.
    pub fn normalize(self) -> Vec3 {
        let n = self.norm();
        if n == 0.0 {
            Vec3::ZERO
        } else {
            self / n
        }
    }

    /// Scale so that the norm equals `l` (if the vector is nonzero).
    pub fn scale_to(self, l: f64) -> Vec3 {
        self.normalize() * l
    }

    /// Component-wise product (Hadamard).
    pub fn mul_add(self, a: f64, b: Vec3) -> Vec3 {
        Vec3::new(
            crate::num::fma(self.x, a, b.x),
            crate::num::fma(self.y, a, b.y),
            crate::num::fma(self.z, a, b.z),
        )
    }

    /// Linear interpolation: `self + t * (o - self)`.
    pub fn lerp(self, o: Vec3, t: f64) -> Vec3 {
        self + (o - self) * t
    }
}

impl From<[f64; 3]> for Vec3 {
    fn from(a: [f64; 3]) -> Self {
        Vec3::new(a[0], a[1], a[2])
    }
}

impl From<Vec3> for [f64; 3] {
    fn from(v: Vec3) -> Self {
        [v.x, v.y, v.z]
    }
}

impl Add for Vec3 {
    type Output = Vec3;
    fn add(self, o: Vec3) -> Vec3 {
        Vec3::new(self.x + o.x, self.y + o.y, self.z + o.z)
    }
}

impl Sub for Vec3 {
    type Output = Vec3;
    fn sub(self, o: Vec3) -> Vec3 {
        Vec3::new(self.x - o.x, self.y - o.y, self.z - o.z)
    }
}

impl AddAssign for Vec3 {
    fn add_assign(&mut self, o: Vec3) {
        *self = *self + o;
    }
}

impl SubAssign for Vec3 {
    fn sub_assign(&mut self, o: Vec3) {
        *self = *self - o;
    }
}

impl Neg for Vec3 {
    type Output = Vec3;
    fn neg(self) -> Vec3 {
        Vec3::new(-self.x, -self.y, -self.z)
    }
}

impl Mul<f64> for Vec3 {
    type Output = Vec3;
    fn mul(self, s: f64) -> Vec3 {
        Vec3::new(self.x * s, self.y * s, self.z * s)
    }
}

impl Mul<Vec3> for f64 {
    type Output = Vec3;
    fn mul(self, v: Vec3) -> Vec3 {
        v * self
    }
}

impl Div<f64> for Vec3 {
    type Output = Vec3;
    fn div(self, s: f64) -> Vec3 {
        Vec3::new(self.x / s, self.y / s, self.z / s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ops() {
        let a = Vec3::new(1.0, 2.0, 3.0);
        let b = Vec3::new(4.0, 5.0, 6.0);
        assert_eq!(a + b, Vec3::new(5.0, 7.0, 9.0));
        assert_eq!(b - a, Vec3::new(3.0, 3.0, 3.0));
        assert_eq!(a * 2.0, Vec3::new(2.0, 4.0, 6.0));
        assert_eq!(-a, Vec3::new(-1.0, -2.0, -3.0));
        assert_eq!(a.get(2), 3.0);
    }

    #[test]
    fn dot_cross_norm() {
        let x = Vec3::new(1.0, 0.0, 0.0);
        let y = Vec3::new(0.0, 1.0, 0.0);
        assert_eq!(x.dot(y), 0.0);
        assert_eq!(x.cross(y), Vec3::new(0.0, 0.0, 1.0));
        assert!((Vec3::new(3.0, 4.0, 0.0).norm() - 5.0).abs() < 1e-14);
        assert!((Vec3::new(1.0, 2.0, 2.0).dist(Vec3::ZERO) - 3.0).abs() < 1e-14);
    }

    #[test]
    fn normalize_zero_is_zero() {
        assert_eq!(Vec3::ZERO.normalize(), Vec3::ZERO);
        let u = Vec3::new(0.0, 0.0, -7.0).normalize();
        assert!((u.norm() - 1.0).abs() < 1e-14);
    }
}
