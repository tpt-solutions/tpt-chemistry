//! Bravais lattices and the direct/reciprocal lattice geometry.
//!
//! A [`Lattice`] is defined by the three conventional cell vectors (row
//! vectors, Å). [`Bravais`] names the 14 Bravais lattice types by their
//! metric symmetries; [`Lattice::classify`] infers the type from the
//! metric within a tolerance.

use tpt_chem_core::num;
use tpt_chem_core::vec3::Vec3;

/// The 14 Bravais lattice types (lattice system + centering).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Bravais {
    /// Cubic primitive.
    CubicP,
    /// Cubic body-centered.
    CubicI,
    /// Cubic face-centered.
    CubicF,
    /// Tetragonal primitive.
    TetragonalP,
    /// Tetragonal body-centered.
    TetragonalI,
    /// Orthorhombic primitive.
    OrthorhombicP,
    /// Orthorhombic body-centered.
    OrthorhombicI,
    /// Orthorhombic face-centered.
    OrthorhombicF,
    /// Orthorhombic base-centered.
    OrthorhombicC,
    /// Monoclinic primitive.
    MonoclinicP,
    /// Monoclinic base-centered.
    MonoclinicC,
    /// Triclinic.
    Triclinic,
    /// Trigonal (rhombohedral).
    TrigonalR,
    /// Hexagonal.
    Hexagonal,
}

impl core::fmt::Display for Bravais {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let s = match self {
            Bravais::CubicP => "cP",
            Bravais::CubicI => "cI",
            Bravais::CubicF => "cF",
            Bravais::TetragonalP => "tP",
            Bravais::TetragonalI => "tI",
            Bravais::OrthorhombicP => "oP",
            Bravais::OrthorhombicI => "oI",
            Bravais::OrthorhombicF => "oF",
            Bravais::OrthorhombicC => "oC",
            Bravais::MonoclinicP => "mP",
            Bravais::MonoclinicC => "mC",
            Bravais::Triclinic => "aP",
            Bravais::TrigonalR => "hR",
            Bravais::Hexagonal => "hP",
        };
        write!(f, "{s}")
    }
}

/// A lattice from three conventional cell vectors (row vectors, Å).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Lattice {
    /// First lattice vector.
    pub a: Vec3,
    /// Second lattice vector.
    pub b: Vec3,
    /// Third lattice vector.
    pub c: Vec3,
}

impl Lattice {
    /// Cubic lattice with edge `a`.
    pub fn cubic(a: f64) -> Self {
        Lattice {
            a: Vec3::new(a, 0.0, 0.0),
            b: Vec3::new(0.0, a, 0.0),
            c: Vec3::new(0.0, 0.0, a),
        }
    }

    /// Orthorhombic lattice.
    pub fn orthorhombic(a: f64, b: f64, c: f64) -> Self {
        Lattice {
            a: Vec3::new(a, 0.0, 0.0),
            b: Vec3::new(0.0, b, 0.0),
            c: Vec3::new(0.0, 0.0, c),
        }
    }

    /// Hexagonal lattice (a = b, γ = 120°).
    pub fn hexagonal(a: f64, c: f64) -> Self {
        Lattice {
            a: Vec3::new(a, 0.0, 0.0),
            b: Vec3::new(-a / 2.0, a * num::sqrt(3.0) / 2.0, 0.0),
            c: Vec3::new(0.0, 0.0, c),
        }
    }

    /// General triclinic lattice from edge lengths and angles (degrees).
    pub fn triclinic(a: f64, b: f64, c: f64, alpha: f64, beta: f64, gamma: f64) -> Self {
        // Convention: c along z; b in the xz-plane; a general.
        let (al, be, ga) = (alpha.to_radians(), beta.to_radians(), gamma.to_radians());
        let va = Vec3::new(a, 0.0, 0.0);
        let vb = Vec3::new(b * ga.cos(), b * ga.sin(), 0.0);
        let cx = c * be.cos();
        let cy = c * (al.cos() - be.cos() * ga.cos()) / ga.sin();
        let cz2 = c * c - cx * cx - cy * cy;
        let vc = Vec3::new(cx, cy, cz2.max(0.0).sqrt());
        Lattice {
            a: va,
            b: vb,
            c: vc,
        }
    }

    /// Cell volume (Å³), signed positive.
    pub fn volume(&self) -> f64 {
        self.a.dot(self.b.cross(self.c)).abs()
    }

    /// Fractional → Cartesian (Å).
    pub fn to_cartesian(&self, frac: Vec3) -> Vec3 {
        self.a * frac.x + self.b * frac.y + self.c * frac.z
    }

    /// Cartesian (Å) → fractional.
    pub fn to_fractional(&self, cart: Vec3) -> Vec3 {
        let r = self.reciprocal_no_2pi();
        Vec3::new(cart.dot(r.a), cart.dot(r.b), cart.dot(r.c))
    }

    /// Reciprocal lattice `bᵢ = 2π (a_j × a_k)/V` (without the 2π in
    /// [`ReciprocalLattice`]; this one carries the 2π, crystallographic
    /// convention uses 1/d — see [`Lattice::reciprocal_no_2pi`]).
    ///
    /// # Panics
    /// Panics on a degenerate (zero-volume) lattice, where the reciprocal
    /// is undefined. The lattice vectors are caller input.
    pub fn reciprocal(&self) -> ReciprocalLattice {
        let v = self.volume_with_sign();
        ReciprocalLattice {
            a: self.b.cross(self.c) * (2.0 * core::f64::consts::PI / v),
            b: self.c.cross(self.a) * (2.0 * core::f64::consts::PI / v),
            c: self.a.cross(self.b) * (2.0 * core::f64::consts::PI / v),
        }
    }

    /// Reciprocal vectors *without* the 2π factor (crystallographic
    /// convention: `aᵢ·aⱼ* = δᵢⱼ`), for fractional↔Cartesian conversion.
    ///
    /// # Panics
    /// Panics on a degenerate (zero-volume) lattice, where the reciprocal
    /// is undefined. The lattice vectors are caller input.
    pub fn reciprocal_no_2pi(&self) -> ReciprocalLattice {
        let v = self.volume_with_sign();
        ReciprocalLattice {
            a: self.b.cross(self.c) / v,
            b: self.c.cross(self.a) / v,
            c: self.a.cross(self.b) / v,
        }
    }

    fn volume_with_sign(&self) -> f64 {
        let v = self.a.dot(self.b.cross(self.c));
        if v == 0.0 {
            panic!("degenerate lattice: zero volume");
        }
        v
    }

    /// Classify the Bravais type from the metric (lengths + angles) within
    /// `tol` (relative for lengths, degrees for angles).
    pub fn classify(&self, tol: f64) -> Bravais {
        let (la, lb, lc) = (self.a.norm(), self.b.norm(), self.c.norm());
        let eq = |x: f64, y: f64| (x - y).abs() <= tol * x.max(y).max(1e-12);
        let angle = |u: Vec3, v: Vec3| num::acos(u.dot(v) / (u.norm() * v.norm())).to_degrees();
        let (al, be, ga) = (
            angle(self.b, self.c),
            angle(self.a, self.c),
            angle(self.a, self.b),
        );
        let ang_eq = |x: f64, y: f64| (x - y).abs() <= tol * 10.0;
        let right = |x: f64| ang_eq(x, 90.0);

        let abc = eq(la, lb) && eq(lb, lc);
        let all_right = right(al) && right(be) && right(ga);

        if abc && all_right {
            return Bravais::CubicP;
        }
        if all_right {
            // Right angles with a = b (but c free) → tetragonal; the cubic
            // case was handled above.
            return if eq(la, lb) {
                Bravais::TetragonalP
            } else {
                Bravais::OrthorhombicP
            };
        }
        if all_right {
            return Bravais::OrthorhombicP;
        }
        if right(al) && right(be) && !right(ga) && eq(la, lb) {
            return Bravais::Hexagonal;
        }
        if right(be) && right(ga) {
            return Bravais::MonoclinicP;
        }
        Bravais::Triclinic
    }
}

/// Reciprocal-lattice vectors (arbitrary 2π convention).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReciprocalLattice {
    /// Reciprocal of a.
    pub a: Vec3,
    /// Reciprocal of b.
    pub b: Vec3,
    /// Reciprocal of c.
    pub c: Vec3,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cubic_volume_and_cartesian() {
        let lat = Lattice::cubic(5.0);
        assert!((lat.volume() - 125.0).abs() < 1e-12);
        let cart = lat.to_cartesian(Vec3::new(0.25, 0.25, 0.25));
        assert!((cart.x - 1.25).abs() < 1e-12);
        assert!((cart.z - 1.25).abs() < 1e-12);
    }

    #[test]
    fn reciprocal_orthogonality() {
        let lat = Lattice::cubic(4.0);
        let r = lat.reciprocal_no_2pi();
        // a·a* = 1, a·b* = 0.
        assert!((lat.a.dot(r.a) - 1.0).abs() < 1e-12);
        assert!(lat.a.dot(r.b).abs() < 1e-12);
        assert!(lat.c.dot(r.b).abs() < 1e-12);
        // Reciprocal lattice constant for cubic: 1/a.
        assert!((r.a.norm() - 1.0 / 4.0).abs() < 1e-12);
    }

    #[test]
    fn roundtrip_fractional() {
        let lat = Lattice::triclinic(4.0, 5.0, 6.0, 75.0, 85.0, 95.0);
        let frac = Vec3::new(0.13, 0.47, 0.89);
        let cart = lat.to_cartesian(frac);
        let back = {
            let r = lat.reciprocal_no_2pi();
            Vec3::new(cart.dot(r.a), cart.dot(r.b), cart.dot(r.c))
        };
        assert!((back.x - frac.x).abs() < 1e-12);
        assert!((back.y - frac.y).abs() < 1e-12);
        assert!((back.z - frac.z).abs() < 1e-12);
    }

    #[test]
    fn classification() {
        assert_eq!(Lattice::cubic(5.0).classify(1e-6), Bravais::CubicP);
        assert_eq!(
            Lattice::orthorhombic(3.0, 4.0, 5.0).classify(1e-6),
            Bravais::OrthorhombicP
        );
        assert_eq!(
            Lattice::hexagonal(3.0, 7.0).classify(1e-6),
            Bravais::Hexagonal
        );
        assert_eq!(
            Lattice::triclinic(4.0, 5.0, 6.0, 70.0, 80.0, 100.0).classify(1e-6),
            Bravais::Triclinic
        );
    }

    #[test]
    fn hexagonal_volume() {
        // V = (√3/2) a² c.
        let lat = Lattice::hexagonal(3.0, 8.0);
        let want = 3.0f64.sqrt() / 2.0 * 9.0 * 8.0;
        assert!((lat.volume() - want).abs() < 1e-12);
    }
}
