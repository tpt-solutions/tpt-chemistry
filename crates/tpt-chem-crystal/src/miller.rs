//! Miller indices `(h k l)` and plane geometry.
//!
//! A [`Miller`] index identifies a lattice plane family; its d-spacing
//! follows from the reciprocal-lattice vector `g = h a* + k b* + l c*` as
//! `d = 1/|g*|` (crystallographic, no 2π).

use crate::bravais::Lattice;
use tpt_chem_core::num;
use tpt_chem_core::vec3::Vec3;

/// Miller index `(h k l)`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Miller {
    /// h index.
    pub h: i32,
    /// k index.
    pub k: i32,
    /// l index.
    pub l: i32,
}

impl Miller {
    /// Construct from integers.
    pub const fn new(h: i32, k: i32, l: i32) -> Self {
        Miller { h, k, l }
    }

    /// Parse "(1 0 -2)" / "(1 0 −2)" style strings.
    pub fn parse(s: &str) -> Option<Miller> {
        let cleaned = s.trim().trim_start_matches('(').trim_end_matches(')');
        let parts: Vec<&str> = cleaned.split_whitespace().collect();
        if parts.len() != 3 {
            return None;
        }
        Some(Miller::new(
            parse_signed(parts[0])?,
            parse_signed(parts[1])?,
            parse_signed(parts[2])?,
        ))
    }

    /// The reciprocal-lattice vector `g = h a* + k b* + l c*` (1/Å,
    /// crystallographic 2π-free convention).
    pub fn reciprocal_vector(&self, lat: &Lattice) -> Vec3 {
        let r = lat.reciprocal_no_2pi();
        r.a * self.h as f64 + r.b * self.k as f64 + r.c * self.l as f64
    }

    /// d-spacing `d = 1/|g*|` in Å.
    pub fn d_spacing(&self, lat: &Lattice) -> f64 {
        1.0 / self.reciprocal_vector(lat).norm()
    }

    /// Bragg angle `θ` (radians) for X-rays of wavelength `λ` (Å):
    /// `nλ = 2d sinθ` with n = 1.
    pub fn bragg_angle(&self, lat: &Lattice, lambda: f64) -> Option<f64> {
        let d = self.d_spacing(lat);
        let s = lambda / (2.0 * d);
        if s <= 0.0 || s > 1.0 {
            return None;
        }
        Some(num::asin(s))
    }

    /// Interplanar d-spacing for a CUBIC lattice: `d = a/√(h²+k²+l²)`.
    pub fn d_cubic(&self, a: f64) -> f64 {
        a / num::sqrt((self.h * self.h + self.k * self.k + self.l * self.l) as f64)
    }

    /// The greatest common divisor of (h, k, l); (2 4 6) → 2.
    pub fn gcd(&self) -> i32 {
        gcd(gcd(self.h.abs(), self.k.abs()), self.l.abs())
    }

    /// Reduce to the coprime form (dividing by the gcd).
    pub fn reduced(&self) -> Miller {
        let g = self.gcd();
        if g <= 1 {
            *self
        } else {
            Miller::new(self.h / g, self.k / g, self.l / g)
        }
    }
}

/// Parse one signed Miller index, accepting ASCII '-' and the Unicode minus.
fn parse_signed(token: &str) -> Option<i32> {
    let t = token.trim();
    if let Some(rest) = t.strip_prefix('-') {
        rest.trim().parse::<i32>().ok().map(|v| -v)
    } else if let Some(rest) = t.strip_prefix('\u{2212}') {
        rest.trim().parse::<i32>().ok().map(|v| -v)
    } else {
        t.parse::<i32>().ok()
    }
}

fn gcd(a: i32, b: i32) -> i32 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cubic_d_spacing() {
        let m = Miller::new(1, 1, 1);
        assert!((m.d_cubic(5.0) - 5.0 / 3.0f64.sqrt()).abs() < 1e-14);
        let m = Miller::new(2, 0, 0);
        assert!((m.d_cubic(5.0) - 2.5).abs() < 1e-14);
    }

    #[test]
    fn general_d_spacing_matches_cubic() {
        let lat = Lattice::cubic(5.431);
        let m = Miller::new(1, 1, 1);
        assert!((m.d_spacing(&lat) - m.d_cubic(5.431)).abs() < 1e-12);
        let m = Miller::new(3, 1, 1);
        assert!((m.d_spacing(&lat) - m.d_cubic(5.431)).abs() < 1e-12);
    }

    #[test]
    fn bragg_angle_si_111() {
        // Si 111 with Cu Kα (1.5406 Å): 2θ = 28.44° → θ ≈ 14.22°.
        let lat = Lattice::cubic(5.431);
        let m = Miller::new(1, 1, 1);
        let theta = m.bragg_angle(&lat, 1.5406).unwrap();
        assert!(
            (theta.to_degrees() - 14.22).abs() < 0.05,
            "{}",
            theta.to_degrees()
        );
    }

    #[test]
    fn bragg_angle_beyond_limit_is_none() {
        let lat = Lattice::cubic(5.0);
        // d = 5/√27 ≈ 0.96 < λ/2 = 0.77? λ = 2 → λ/2 = 1 > d → none.
        let m = Miller::new(3, 3, 3);
        assert!(m.bragg_angle(&lat, 2.0).is_none());
    }

    #[test]
    fn reduction_and_gcd() {
        let m = Miller::new(2, 4, 6);
        assert_eq!(m.gcd(), 2);
        assert_eq!(m.reduced(), Miller::new(1, 2, 3));
        let m = Miller::new(1, 1, 1);
        assert_eq!(m.gcd(), 1);
        assert_eq!(m.reduced(), m);
    }

    #[test]
    fn parse_formats() {
        assert_eq!(Miller::parse("(1 0 -2)"), Some(Miller::new(1, 0, -2)));
        assert_eq!(Miller::parse("(1 0 −2)"), Some(Miller::new(1, 0, -2)));
        assert_eq!(Miller::parse("1 0 2"), Some(Miller::new(1, 0, 2)));
        assert_eq!(Miller::parse("(1 0)"), None);
    }
}
