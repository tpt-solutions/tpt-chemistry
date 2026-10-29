//! Space-group symmetry operations and Wyckoff orbit generation.
//!
//! A [`SymmetryOp`] is a 3×3 rotation (integer entries for the point-group
//! part) plus a fractional translation. A [`SpaceGroup`] carries the
//! operations of one of the included groups; Wyckoff orbits follow by
//! applying every operation to a seed position and reducing modulo 1.

use std::string::String;
use std::vec::Vec;

/// Fractional translation vector.
pub type FracTranslation = [f64; 3];

/// A symmetry operation: integer rotation `r` plus fractional translation
/// `t`, acting on fractional coordinates as `x' = r·x + t`.
#[derive(Clone, Debug, PartialEq)]
pub struct SymmetryOp {
    /// 3×3 integer rotation matrix (row-major).
    pub r: [[i32; 3]; 3],
    /// Fractional translation.
    pub t: FracTranslation,
}

impl SymmetryOp {
    /// The identity operation.
    pub fn identity() -> Self {
        SymmetryOp {
            r: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            t: [0.0; 3],
        }
    }

    /// Inversion through the origin with optional translation.
    pub fn inversion(t: FracTranslation) -> Self {
        SymmetryOp {
            r: [[-1, 0, 0], [0, -1, 0], [0, 0, -1]],
            t,
        }
    }

    /// Apply to a fractional point, wrapping into [0, 1).
    pub fn apply(&self, x: [f64; 3]) -> [f64; 3] {
        let mut out = [0.0; 3];
        for (row, o) in out.iter_mut().enumerate() {
            *o = (self.r[row][0] as f64) * x[0]
                + (self.r[row][1] as f64) * x[1]
                + (self.r[row][2] as f64) * x[2]
                + self.t[row];
            *o -= o.floor();
            // Snap values within fp noise of 1 back to 0.
            if (*o - 1.0).abs() < 1e-10 {
                *o = 0.0;
            }
        }
        out
    }
}

/// Operations of all 230 space groups in their standard settings, generated
/// by `tools/gen_spacegroups.py` from the spglib database.
const SPACEGROUP_DATA: &str = include_str!("../data/spacegroups.txt");

/// A named space group: identity plus its operations.
#[derive(Clone, Debug)]
pub struct SpaceGroup {
    /// International Tables number (1–230), when known.
    pub number: Option<u16>,
    /// Hermann-Mauguin symbol.
    pub symbol: String,
    /// Full operation set (includes the identity).
    pub operations: Vec<SymmetryOp>,
}

impl SpaceGroup {
    /// The standard-setting space group with International Tables number
    /// `number` (1–230), with all of its general-position operations.
    ///
    /// Covers all 230 groups. The operations come from the generated table
    /// in `data/spacegroups.txt`. `None` for numbers outside 1–230.
    pub fn from_number(number: u16) -> Option<SpaceGroup> {
        let mut symbol = None;
        let mut operations = Vec::new();
        let mut in_group = false;
        for line in SPACEGROUP_DATA.lines() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            if let Some(header) = line.strip_prefix("@ ") {
                if in_group {
                    break;
                }
                // "@ <number> <hall> <Hermann-Mauguin>"; underscores mark subscripts.
                let mut fields = header.splitn(3, ' ');
                let n: u16 = fields.next()?.parse().ok()?;
                let _hall = fields.next()?;
                if n == number {
                    in_group = true;
                    symbol = Some(fields.next()?.replace('_', ""));
                }
                continue;
            }
            if !in_group {
                continue;
            }
            let v: Vec<i64> = line
                .split_whitespace()
                .map(|x| x.parse().ok())
                .collect::<Option<_>>()?;
            if v.len() != 12 {
                return None;
            }
            let r = [
                [v[0] as i32, v[1] as i32, v[2] as i32],
                [v[3] as i32, v[4] as i32, v[5] as i32],
                [v[6] as i32, v[7] as i32, v[8] as i32],
            ];
            let t = [v[9] as f64 / 12.0, v[10] as f64 / 12.0, v[11] as f64 / 12.0];
            operations.push(SymmetryOp { r, t });
        }
        Some(SpaceGroup {
            number: Some(number),
            symbol: symbol?,
            operations,
        })
    }

    /// The operations of a space group by Hermann-Mauguin symbol, for the
    /// standard settings of the groups that have a common short symbol
    /// (`P1`, `P-1`, `P21`, `P21/c`, `P212121`, `Pna21`, `P4`, `P-3`,
    /// `P6/mmm`, `Fm-3m`). These come from the generated table, so they match
    /// [`SpaceGroup::from_number`]. `None` for other symbols; use
    /// `from_number` for all 230 groups.
    pub fn by_symbol(symbol: &str) -> Option<SpaceGroup> {
        let number = standard_number(symbol)?;
        let mut group = SpaceGroup::from_number(number)?;
        group.symbol = symbol.into();
        Some(group)
    }

    /// The orbit of a fractional position: every symmetry image reduced
    /// into [0,1), deduplicated within `eps`.
    pub fn orbit(&self, seed: [f64; 3], eps: f64) -> Vec<[f64; 3]> {
        let mut orbit: Vec<[f64; 3]> = Vec::new();
        for op in &self.operations {
            let image = op.apply(seed);
            if !orbit.iter().any(|p| {
                p.iter().zip(&image).all(|(a, b)| {
                    let d = (a - b).abs();
                    d.min((1.0 - d).abs()) < eps
                })
            }) {
                orbit.push(image);
            }
        }
        orbit
    }
}

/// International Tables number of one of the groups in [`SpaceGroup::by_symbol`].
fn standard_number(symbol: &str) -> Option<u16> {
    Some(match symbol {
        "P1" => 1,
        "P-1" => 2,
        "P21" => 4,
        "P21/c" => 14,
        "P212121" => 19,
        "Pna21" => 33,
        "P4" => 75,
        "P-3" => 147,
        "P6/mmm" => 191,
        "Fm-3m" => 225,
        _ => return None,
    })
}

/// Compose two symmetry operations: `(op2 ∘ op1)(x) = op2(op1(x))`.
pub fn compose(op2: &SymmetryOp, op1: &SymmetryOp) -> SymmetryOp {
    let mut r = [[0i32; 3]; 3];
    #[allow(clippy::needless_range_loop)]
    for i in 0..3 {
        for j in 0..3 {
            r[i][j] = (0..3).map(|k| op2.r[i][k] * op1.r[k][j]).sum();
        }
    }
    let mut t = [0.0; 3];
    for (i, ti) in t.iter_mut().enumerate() {
        #[allow(clippy::needless_range_loop)]
        for j in 0..3 {
            *ti += (op2.r[i][j] as f64) * op1.t[j];
        }
        *ti += op2.t[i];
        *ti -= ti.floor();
    }
    SymmetryOp { r, t }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn p1_orbit_is_single() {
        let g = SpaceGroup::by_symbol("P1").unwrap();
        assert_eq!(g.operations.len(), 1);
        let orbit = g.orbit([0.13, 0.27, 0.41], 1e-6);
        assert_eq!(orbit.len(), 1);
    }

    #[test]
    fn pminus1_doubles_inversion_center() {
        let g = SpaceGroup::by_symbol("P-1").unwrap();
        let orbit = g.orbit([0.1, 0.2, 0.3], 1e-6);
        assert_eq!(orbit.len(), 2);
        let inv = orbit.iter().find(|p| **p != [0.1, 0.2, 0.3]).unwrap();
        // Inversion through origin: (−x, −y, −z) mod 1.
        assert!((inv[0] - 0.9).abs() < 1e-10);
        assert!((inv[1] - 0.8).abs() < 1e-10);
        assert!((inv[2] - 0.7).abs() < 1e-10);
    }

    #[test]
    fn p21_screw_shifts_half() {
        let g = SpaceGroup::by_symbol("P21").unwrap();
        let orbit = g.orbit([0.1, 0.2, 0.3], 1e-6);
        assert_eq!(orbit.len(), 2);
        let other = orbit.iter().find(|p| **p != [0.1, 0.2, 0.3]).unwrap();
        // (−x, y+½, −z): x → −0.1 ≡ 0.9, y = 0.7, z = 0.7.
        assert!((other[0] - 0.9).abs() < 1e-10);
        assert!((other[1] - 0.7).abs() < 1e-10);
        assert!((other[2] - 0.7).abs() < 1e-10);
    }

    #[test]
    fn p21c_general_position_has_4() {
        let g = SpaceGroup::by_symbol("P21/c").unwrap();
        assert_eq!(g.operations.len(), 4);
        let orbit = g.orbit([0.11, 0.22, 0.33], 1e-6);
        assert_eq!(orbit.len(), 4);
    }

    #[test]
    fn wyckoff_special_position_collapses() {
        // On the inversion center of P-1 the orbit collapses to 1.
        let g = SpaceGroup::by_symbol("P-1").unwrap();
        let orbit = g.orbit([0.0, 0.0, 0.0], 1e-6);
        assert_eq!(orbit.len(), 1);
    }

    #[test]
    fn cubic_point_group_from_table() {
        // Pm-3m (221) is primitive, so its operations are the m-3m point group
        // itself: 48 operations with no translations.
        let ops = SpaceGroup::from_number(221).unwrap().operations;
        assert_eq!(ops.len(), 48);
        assert!(ops.iter().all(|op| op.t == [0.0; 3]));
        // All are orthogonal: each row has exactly one ±1.
        for op in &ops {
            for row in op.r {
                assert_eq!(row.iter().filter(|&&v| v != 0).count(), 1);
                assert!(row.iter().all(|&v| v == 0 || v.abs() == 1));
            }
        }
    }

    #[test]
    fn fm3m_general_orbit_is_192() {
        // F m -3 m general Wyckoff 192 positions (48 ops × 4 centerings).
        let g = SpaceGroup::by_symbol("Fm-3m").unwrap();
        assert_eq!(g.operations.len(), 192);
        let orbit = g.orbit([0.123, 0.157, 0.211], 1e-6);
        assert_eq!(orbit.len(), 192);
    }

    #[test]
    fn p6mmm_orbit_24() {
        // 6/mmm has point-group order 24 and P6/mmm is primitive.
        let g = SpaceGroup::by_symbol("P6/mmm").unwrap();
        assert_eq!(g.operations.len(), 24);
        let orbit = g.orbit([0.13, 0.29, 0.41], 1e-6);
        assert_eq!(orbit.len(), 24);
    }

    #[test]
    fn compose_identity_is_idempotent() {
        let g = SymmetryOp {
            r: [[0, -1, 0], [1, 0, 0], [0, 0, 1]], // 90° rotation
            t: [0.25, 0.0, 0.0],
        };
        let composed = compose(&g, &g); // 180° + 0.5 translation
        let x = [0.1, 0.2, 0.3];
        let directly = composed.apply(x);
        let once: [f64; 3] = g.apply(x);
        let twice = g.apply(once);
        for (a, b) in directly.iter().zip(&twice) {
            let d = (a - b).abs();
            assert!(d.min((1.0 - d).abs()) < 1e-10);
        }
    }
}
