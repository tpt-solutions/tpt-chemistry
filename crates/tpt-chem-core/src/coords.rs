//! Compile-time coordinate-frame safety and internal (Z-matrix) coordinates.
//!
//! The [`Frame`] trait's sealed marker types [`Cartesian`] and [`Internal`]
//! parameterize [`Position`]: a `Position<Cartesian>` (Å) and a
//! `Position<Internal>` (r, θ, φ triple) are distinct types, so mixing them
//! in an expression is a compile error — spec.txt §2's type-state principle.
//!
//! Converting a Z-matrix to Cartesian coordinates uses the natural extension
//! reference frame (NeRF) construction — see [`zmatrix_to_cartesian`].
//!
//! ```
//! use tpt_chem_core::coords::{Cartesian, Internal, Position};
//!
//! // Frame-tagged positions never mix:
//! let cart = Position::<Cartesian>::new([1.0, 2.0, 3.0].into());
//! let int = Position::<Internal>::spherical(1.5, 0.5, 0.25);
//! let sum = cart + Position::<Cartesian>::new([0.1, 0.0, 0.0].into());
//! assert!((sum.into_raw().x - 1.1).abs() < 1e-12);
//! assert!((int.r() - 1.5).abs() < 1e-12);
//! ```

use alloc::vec::Vec;
use core::marker::PhantomData;

use crate::num;
use crate::vec3::Vec3;

/// Sealed marker trait for coordinate frames.
pub trait Frame: private::Sealed {}

mod private {
    pub trait Sealed {}
    impl Sealed for super::Cartesian {}
    impl Sealed for super::Internal {}
}

/// Cartesian coordinate frame (Å).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Cartesian {}

/// Internal coordinate frame: (r, θ, φ) — distance in Å, polar and azimuthal
/// angles in radians, packed into the vector's x/y/z fields.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Internal {}

impl Frame for Cartesian {}
impl Frame for Internal {}

/// A position tagged with its coordinate frame.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Position<F: Frame> {
    v: Vec3,
    _frame: PhantomData<fn() -> F>,
}

impl<F: Frame> Position<F> {
    /// Wrap a raw vector (already expressed in frame `F`).
    pub fn from_raw(v: Vec3) -> Self {
        Position {
            v,
            _frame: PhantomData,
        }
    }

    /// Borrow the raw vector.
    pub fn raw(&self) -> &Vec3 {
        &self.v
    }

    /// Unwrap the raw vector.
    pub fn into_raw(self) -> Vec3 {
        self.v
    }
}

impl Position<Cartesian> {
    /// Construct from Cartesian components (Å).
    pub fn new(v: Vec3) -> Self {
        Position {
            v,
            _frame: PhantomData,
        }
    }
}

impl Position<Internal> {
    /// Construct from spherical coordinates: `r` in Å, `theta`/`phi` in rad.
    pub fn spherical(r: f64, theta: f64, phi: f64) -> Self {
        Position {
            v: Vec3::new(r, theta, phi),
            _frame: PhantomData,
        }
    }

    /// Distance component `r` (Å).
    pub fn r(&self) -> f64 {
        self.v.x
    }

    /// Polar angle `θ` (rad).
    pub fn theta(&self) -> f64 {
        self.v.y
    }

    /// Azimuthal angle `φ` (rad).
    pub fn phi(&self) -> f64 {
        self.v.z
    }
}

impl<F: Frame> core::ops::Add for Position<F> {
    type Output = Position<F>;
    fn add(self, o: Position<F>) -> Position<F> {
        Position::from_raw(self.v + o.v)
    }
}

impl<F: Frame> core::ops::Sub for Position<F> {
    type Output = Position<F>;
    fn sub(self, o: Position<F>) -> Position<F> {
        Position::from_raw(self.v - o.v)
    }
}

/// Bond angle at the middle point `b`, in radians (0..π).
///
/// # Panics
/// If `a == b` or `c == b` (zero-length arm).
pub fn angle(a: Vec3, b: Vec3, c: Vec3) -> f64 {
    let ab = a - b;
    let cb = c - b;
    let cos_theta = ab.dot(cb) / (ab.norm() * cb.norm());
    num::acos(cos_theta)
}

/// Dihedral (torsion) angle of the chain `a–b–c–d`, in radians (−π..π),
/// IUPAC sign convention.
///
/// # Panics
/// If any consecutive pair coincides (degenerate chain).
pub fn dihedral(a: Vec3, b: Vec3, c: Vec3, d: Vec3) -> f64 {
    // Praxeolitic's formulation.
    let b0 = a - b;
    let b1 = c - b;
    let b2 = d - c;
    let b1 = b1 / b1.norm();
    let v = b0 - b1 * b0.dot(b1);
    let w = b2 - b1 * b2.dot(b1);
    let x = v.dot(w);
    let y = b1.cross(v).dot(w);
    -num::atan2(y, x)
}

/// One row of a Z-matrix: an element and up to three internal-coordinate
/// references.
///
/// * Row 0 needs no references (placed at the origin).
/// * Row 1 needs `dist = (ref, r)`.
/// * Row 2 needs `dist` and `angle = (dist_ref, other, θ)`; θ is measured at
///   the distance reference.
/// * Rows ≥ 3 need `dist`, `angle`, and
///   `dihedral = (a, b, dist_ref, φ)`; φ is the dihedral a–b–dist_ref–new.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ZMatrixRow {
    /// Atomic number of the atom this row places.
    pub z: u8,
    /// Distance reference: `(atom index, distance in Å)`.
    pub dist: Option<(usize, f64)>,
    /// Angle reference: `(distance ref, second ref, angle in rad)`; the
    /// angle is measured at the distance reference.
    pub angle: Option<(usize, usize, f64)>,
    /// Dihedral reference: `(first, second, distance ref, angle in rad)` —
    /// the torsion `first–second–dist_ref–new`.
    pub dihedral: Option<(usize, usize, usize, f64)>,
}

/// Errors from [`zmatrix_to_cartesian`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoordError {
    /// A row references an atom that is not placed yet.
    ForwardReference,
    /// A row is missing the references its position requires.
    MissingReference,
    /// Degenerate geometry: three reference atoms are collinear where a
    /// dihedral is needed, or a reference distance is zero.
    Degenerate,
}

impl core::fmt::Display for CoordError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            CoordError::ForwardReference => write!(f, "reference to a not-yet-placed atom"),
            CoordError::MissingReference => write!(f, "row lacks required references"),
            CoordError::Degenerate => write!(f, "degenerate reference geometry"),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for CoordError {}

/// Convert a Z-matrix to Cartesian coordinates (Å) via the natural extension
/// reference frame algorithm.
///
/// The first atom lands at the origin, the second along +z, the third in the
/// xz-plane.
///
/// # Errors
/// [`CoordError`] if any row is missing references or refers to a
/// not-yet-placed atom, or the reference geometry is degenerate.
pub fn zmatrix_to_cartesian(rows: &[ZMatrixRow]) -> Result<Vec<Vec3>, CoordError> {
    let mut out: Vec<Vec3> = Vec::with_capacity(rows.len());
    for (i, row) in rows.iter().enumerate() {
        let p = match i {
            0 => Vec3::ZERO,
            1 => {
                let (ir, r) = row.dist.ok_or(CoordError::MissingReference)?;
                let base = *out.get(ir).ok_or(CoordError::ForwardReference)?;
                base + Vec3::new(0.0, 0.0, r)
            }
            2 => {
                let (ir, r) = row.dist.ok_or(CoordError::MissingReference)?;
                let (vertex, ja, theta) = row.angle.ok_or(CoordError::MissingReference)?;
                if vertex != ir {
                    // Convention: the angle vertex is the distance reference.
                    return Err(CoordError::MissingReference);
                }
                let base = *out.get(ir).ok_or(CoordError::ForwardReference)?;
                let other = *out.get(ja).ok_or(CoordError::ForwardReference)?;
                let dir = other - base;
                let len = dir.norm();
                if len < 1e-12 {
                    return Err(CoordError::Degenerate);
                }
                let dir = dir / len;
                // Choose a deterministic perpendicular axis to fix azimuth.
                let axis = if dir.x.abs() < 0.9 {
                    Vec3::new(1.0, 0.0, 0.0)
                } else {
                    Vec3::new(0.0, 1.0, 0.0)
                };
                let perp = axis.cross(dir).normalize();
                base + dir * (r * num::cos(theta)) + perp * (r * num::sin(theta))
            }
            _ => {
                let (ir, r) = row.dist.ok_or(CoordError::MissingReference)?;
                let (vertex, ja, theta) = row.angle.ok_or(CoordError::MissingReference)?;
                let (da, db, dc, phi) = row.dihedral.ok_or(CoordError::MissingReference)?;
                // The angle's second reference must be the dihedral's
                // middle atom (standard z-matrix convention K–J–I–new).
                if vertex != ir || dc != ir || ja != db {
                    return Err(CoordError::MissingReference);
                }
                let a = *out.get(da).ok_or(CoordError::ForwardReference)?;
                let b = *out.get(db).ok_or(CoordError::ForwardReference)?;
                let c = *out.get(dc).ok_or(CoordError::ForwardReference)?;
                // NeRF: angle θ at `c` between `b` and the new atom;
                // dihedral φ = a–b–c–new.
                let bc = c - b;
                let bc_len = bc.norm();
                if bc_len < 1e-12 {
                    return Err(CoordError::Degenerate);
                }
                let bc = bc / bc_len;
                let ab = b - a;
                let n = ab.cross(bc);
                if n.norm() < 1e-8 {
                    return Err(CoordError::Degenerate);
                }
                let n = n.normalize();
                // Coefficients place the new atom at torsion
                // `a–b–c–new = φ` (IUPAC convention, numerically calibrated
                // against [`dihedral`]).
                let d0 = -r * num::cos(theta);
                let d1 = -r * num::sin(theta) * num::sin(phi);
                let d2 = -r * num::sin(theta) * num::cos(phi);
                let bcn = bc.cross(n);
                c + bc * d0 + n * d1 + bcn * d2
            }
        };
        out.push(p);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Water z-matrix: O at origin, H 0.9572 Å along z from O,
    /// H at 0.9572 Å / 104.52°.
    fn water_rows() -> Vec<ZMatrixRow> {
        let r = 0.9572;
        let th = 104.52f64.to_radians();
        vec![
            ZMatrixRow {
                z: 8,
                dist: None,
                angle: None,
                dihedral: None,
            },
            ZMatrixRow {
                z: 1,
                dist: Some((0, r)),
                angle: None,
                dihedral: None,
            },
            ZMatrixRow {
                z: 1,
                dist: Some((0, r)),
                angle: Some((0, 1, th)),
                dihedral: None,
            },
        ]
    }

    #[test]
    fn water_geometry_from_zmatrix() {
        let xyz = zmatrix_to_cartesian(&water_rows()).unwrap();
        assert_eq!(xyz.len(), 3);
        assert_eq!(xyz[0], Vec3::ZERO);
        assert!((xyz[1].dist(xyz[0]) - 0.9572).abs() < 1e-12);
        let th = angle(xyz[2], xyz[0], xyz[1]);
        assert!((th - 104.52f64.to_radians()).abs() < 1e-12);
    }

    #[test]
    fn nerf_roundtrip_dihedral_chain() {
        // Build a 5-atom chain with fixed internal coordinates and verify
        // every internal coordinate round-trips through the Cartesian frame.
        let r = 1.54;
        let th = 109.47f64.to_radians();
        let phis = [60f64.to_radians(), -60f64.to_radians(), 180f64.to_radians()];
        let mut rows = vec![
            ZMatrixRow {
                z: 6,
                dist: None,
                angle: None,
                dihedral: None,
            },
            ZMatrixRow {
                z: 6,
                dist: Some((0, r)),
                angle: None,
                dihedral: None,
            },
            ZMatrixRow {
                z: 6,
                dist: Some((1, r)),
                angle: Some((1, 0, th)),
                dihedral: None,
            },
        ];
        for (k, &phi) in phis.iter().enumerate() {
            let base = k + 2;
            rows.push(ZMatrixRow {
                z: 1,
                dist: Some((base, r)),
                angle: Some((base, base - 1, th)),
                dihedral: Some((base - 2, base - 1, base, phi)),
            });
        }
        let xyz = zmatrix_to_cartesian(&rows).unwrap();
        for i in 1..xyz.len() {
            assert!((xyz[i].dist(xyz[i - 1]) - r).abs() < 1e-10, "distance {i}");
        }
        for i in 2..xyz.len() {
            let a = angle(xyz[i], xyz[i - 1], xyz[i - 2]);
            assert!((a - th).abs() < 1e-9, "angle {i}");
        }
        for (k, &phi) in phis.iter().enumerate() {
            let i = k + 3;
            let d = dihedral(xyz[i - 3], xyz[i - 2], xyz[i - 1], xyz[i]);
            // Wrap-aware angular comparison (−π and +π are the same angle).
            let diff = (d - phi + core::f64::consts::PI).rem_euclid(2.0 * core::f64::consts::PI)
                - core::f64::consts::PI;
            assert!(diff.abs() < 1e-9, "dihedral {i}: {d} vs {phi}");
        }
    }

    #[test]
    fn dihedral_trans_is_pi() {
        // Planar all-trans zigzag chain: torsion magnitude is 180°.
        let pts = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.5, 0.0, 2.0),
            Vec3::new(0.0, 0.0, 3.0),
        ];
        let d = dihedral(pts[0], pts[1], pts[2], pts[3]);
        assert!(
            ((d.abs() - core::f64::consts::PI).abs() < 1e-9),
            "dihedral was {d}"
        );
    }

    #[test]
    fn errors() {
        // Forward reference: row 1 points at atom 3, which does not exist.
        let rows = vec![
            ZMatrixRow {
                z: 1,
                dist: None,
                angle: None,
                dihedral: None,
            },
            ZMatrixRow {
                z: 1,
                dist: Some((3, 1.0)),
                angle: None,
                dihedral: None,
            },
        ];
        assert!(matches!(
            zmatrix_to_cartesian(&rows),
            Err(CoordError::ForwardReference)
        ));
        // Missing distance on row 1 (row 0 itself needs none).
        let rows = vec![
            ZMatrixRow {
                z: 1,
                dist: None,
                angle: None,
                dihedral: None,
            },
            ZMatrixRow {
                z: 1,
                dist: None,
                angle: None,
                dihedral: None,
            },
        ];
        assert!(matches!(
            zmatrix_to_cartesian(&rows),
            Err(CoordError::MissingReference)
        ));
    }

    #[test]
    fn internal_positions_stay_internal() {
        let p = Position::<Internal>::spherical(2.0, 1.0, 0.5);
        assert!((p.r() - 2.0).abs() < 1e-15);
        assert!((p.theta() - 1.0).abs() < 1e-15);
        assert!((p.phi() - 0.5).abs() < 1e-15);
    }
}
