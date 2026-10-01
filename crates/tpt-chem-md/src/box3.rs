//! Orthorhombic simulation box with the minimum-image convention.

use tpt_chem_core::vec3::Vec3;

/// An orthorhombic periodic box with edges along the axes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Box3 {
    /// Box edge lengths (Å).
    pub length: Vec3,
}

impl Box3 {
    /// A cubic box with the given edge length (Å).
    pub fn cubic(l: f64) -> Self {
        Box3 {
            length: Vec3::new(l, l, l),
        }
    }

    /// Smallest-image displacement: wraps each component of `d` into
    /// `[-L/2, L/2)`.
    pub fn min_image(&self, d: Vec3) -> Vec3 {
        Vec3::new(
            min_image_component(d.x, self.length.x),
            min_image_component(d.y, self.length.y),
            min_image_component(d.z, self.length.z),
        )
    }

    /// Wrap a position back into `[0, L)`.
    pub fn wrap(&self, r: Vec3) -> Vec3 {
        Vec3::new(
            wrap_component(r.x, self.length.x),
            wrap_component(r.y, self.length.y),
            wrap_component(r.z, self.length.z),
        )
    }

    /// Box volume (Å³).
    pub fn volume(&self) -> f64 {
        self.length.x * self.length.y * self.length.z
    }
}

#[inline]
fn min_image_component(d: f64, l: f64) -> f64 {
    if l <= 0.0 {
        return d;
    }
    let mut v = d - (d / l).round() * l;
    // Numerical safety at exactly ±L/2.
    if v >= l * 0.5 {
        v -= l;
    } else if v < -l * 0.5 {
        v += l;
    }
    v
}

#[inline]
fn wrap_component(r: f64, l: f64) -> f64 {
    if l <= 0.0 {
        return r;
    }
    let mut v = r - (r / l).floor() * l;
    if v >= l {
        v -= l;
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn min_image_wraps() {
        let b = Box3::cubic(10.0);
        assert!((b.min_image(Vec3::new(11.0, 0.0, 0.0)).x - 1.0).abs() < 1e-12);
        assert!((b.min_image(Vec3::new(-11.0, 0.0, 0.0)).x + 1.0).abs() < 1e-12);
        assert!((b.min_image(Vec3::new(4.9, 0.0, 0.0)).x - 4.9).abs() < 1e-12);
        assert!((b.min_image(Vec3::new(-4.9, 0.0, 0.0)).x + 4.9).abs() < 1e-12);
        assert!((b.min_image(Vec3::new(5.1, 0.0, 0.0)).x + 4.9).abs() < 1e-12);
    }

    #[test]
    fn wrap_into_box() {
        let b = Box3::cubic(10.0);
        let w = b.wrap(Vec3::new(23.7, -3.2, 10.0));
        assert!(w.x >= 0.0 && w.x < 10.0);
        assert!(w.y >= 0.0 && w.y < 10.0);
        assert!(w.z >= 0.0 && w.z < 10.0);
        assert!((w.x - 3.7).abs() < 1e-12);
        assert!((w.y - 6.8).abs() < 1e-12);
    }

    #[test]
    fn min_image_idempotent() {
        let b = Box3::cubic(8.0);
        let d = Vec3::new(30.3, -17.2, 64.9);
        let m1 = b.min_image(d);
        let m2 = b.min_image(m1);
        assert_eq!(m1, m2);
    }
}
