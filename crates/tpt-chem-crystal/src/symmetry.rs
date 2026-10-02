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

/// A named space group: identity plus its operations.
#[derive(Clone, Debug)]
pub struct SpaceGroup {
    /// Hermann-Mauguin symbol.
    pub symbol: String,
    /// Full operation set (includes the identity).
    pub operations: Vec<SymmetryOp>,
}

impl SpaceGroup {
    /// The operations of one of the 11 included space groups by
    /// Hermann-Mauguin symbol. `None` for groups outside the table.
    pub fn by_symbol(symbol: &str) -> Option<SpaceGroup> {
        let ops = match symbol {
            "P1" => vec![SymmetryOp::identity()],
            "P-1" => vec![SymmetryOp::identity(), SymmetryOp::inversion([0.0; 3])],
            "P21" => {
                // 2₁ screw along b: (x, y+½, −z).
                vec![
                    SymmetryOp::identity(),
                    SymmetryOp {
                        r: [[1, 0, 0], [0, 1, 0], [0, 0, -1]],
                        t: [0.0, 0.5, 0.0],
                    },
                ]
            }
            "P21/c" => {
                // Standard setting: (x,y,z), (−x, y+½, −z+½), (−x,−y,−z), (x, −y+½, z+½).
                vec![
                    SymmetryOp::identity(),
                    SymmetryOp {
                        r: [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
                        t: [0.0, 0.5, 0.5],
                    },
                    SymmetryOp::inversion([0.0; 3]),
                    SymmetryOp {
                        r: [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
                        t: [0.0, 0.5, 0.5],
                    },
                ]
            }
            "P212121" => {
                // 2₁ screws along a, b, c.
                vec![
                    SymmetryOp::identity(),
                    SymmetryOp {
                        r: [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
                        t: [0.5, 0.5, 0.0],
                    },
                    SymmetryOp {
                        r: [[-1, 0, 0], [0, 1, 0], [0, 0, -1]],
                        t: [0.5, 0.0, 0.5],
                    },
                    SymmetryOp {
                        r: [[1, 0, 0], [0, -1, 0], [0, 0, -1]],
                        t: [0.0, 0.5, 0.5],
                    },
                ]
            }
            "Pna21" => {
                // n-glide ⊥ a, a-glide ⊥ b, 2₁ along c (origin choice 1 simplified):
                // standard operations for Pna2₁.
                vec![
                    SymmetryOp::identity(),
                    SymmetryOp {
                        r: [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
                        t: [0.0, 0.0, 0.5],
                    },
                    SymmetryOp {
                        r: [[1, 0, 0], [0, -1, 0], [0, 0, 1]],
                        t: [0.5, 0.5, 0.5],
                    },
                    SymmetryOp {
                        r: [[-1, 0, 0], [0, 1, 0], [0, 0, 1]],
                        t: [0.5, 0.5, 0.0],
                    },
                ]
            }
            "P4" => {
                // 4-fold rotation about c: (x,y,z), (−y,x,z), (−x,−y,z), (y,−x,z).
                let r90 = |s: i32| SymmetryOp {
                    r: [[0, -s, 0], [s, 0, 0], [0, 0, 1]],
                    t: [0.0; 3],
                };
                vec![
                    SymmetryOp::identity(),
                    r90(1),
                    SymmetryOp {
                        r: [[-1, 0, 0], [0, -1, 0], [0, 0, 1]],
                        t: [0.0; 3],
                    },
                    r90(-1),
                ]
            }
            "P-3" => vec![
                SymmetryOp::identity(),
                SymmetryOp {
                    r: [[0, -1, 0], [1, -1, 0], [0, 0, 1]],
                    t: [0.0; 3],
                },
                SymmetryOp {
                    r: [[-1, 1, 0], [-1, 0, 0], [0, 0, 1]],
                    t: [0.0; 3],
                },
                SymmetryOp::inversion([0.0; 3]),
                SymmetryOp {
                    r: [[0, 1, 0], [1, 1, 0], [0, 0, -1]],
                    t: [0.0; 3],
                },
                SymmetryOp {
                    r: [[1, -1, 0], [1, 0, 0], [0, 0, -1]],
                    t: [0.0; 3],
                },
            ],
            "P6/mmm" => {
                // 6/m m m: 6-fold rotations + inversion × 12 operations.
                let mut ops = vec![SymmetryOp::identity()];
                // Hexagonal rotations about c (60° steps, hexagonal axes).
                let rot60 = SymmetryOp {
                    r: [[1, -1, 0], [1, 0, 0], [0, 0, 1]],
                    t: [0.0; 3],
                };
                let mut cur = SymmetryOp::identity();
                for _ in 0..5 {
                    cur = compose(&rot60, &cur);
                    ops.push(cur.clone());
                }
                // Inversion × all.
                let inv_ops = ops.clone();
                for op in inv_ops {
                    ops.push(compose(&SymmetryOp::inversion([0.0; 3]), &op));
                }
                ops
            }
            "Fm-3m" => {
                // F m -3 m: 48 point operations × F-centering translations.
                let point_ops = cubic_point_ops();
                let mut ops = Vec::new();
                for op in &point_ops {
                    for &ct in &[
                        [0.0, 0.0, 0.0],
                        [0.0, 0.5, 0.5],
                        [0.5, 0.0, 0.5],
                        [0.5, 0.5, 0.0],
                    ] {
                        ops.push(SymmetryOp {
                            r: op.r,
                            t: [op.t[0] + ct[0], op.t[1] + ct[1], op.t[2] + ct[2]],
                        });
                    }
                }
                ops
            }
            _ => return None,
        };
        Some(SpaceGroup {
            symbol: symbol.into(),
            operations: ops,
        })
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

/// The 48 cubic point-group operations (m-3m point part, translations zero).
fn cubic_point_ops() -> Vec<SymmetryOp> {
    // All 3×3 signed permutation matrices with det = ±1: 48 operations.
    fn permute(perm: &mut [usize; 3], k: usize, signs: &[i32; 3], ops: &mut Vec<SymmetryOp>) {
        if k == 3 {
            let mut r = [[0i32; 3]; 3];
            for (row, &col) in perm.iter().enumerate() {
                r[row][col] = signs[row];
            }
            ops.push(SymmetryOp { r, t: [0.0; 3] });
            return;
        }
        for kk in 0..3 {
            if !perm[..k].contains(&kk) {
                perm[k] = kk;
                permute(perm, k + 1, signs, ops);
            }
        }
    }
    let mut ops = Vec::new();
    let mut signs = [1i32; 3];
    // Iterate all 8 sign combinations.
    for mask in 0..8usize {
        for (row, s) in signs.iter_mut().enumerate() {
            *s = if (mask >> row) & 1 == 0 { 1 } else { -1 };
        }
        permute(&mut [0usize; 3], 0, &signs, &mut ops);
    }
    ops
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
        // (x, y+½, −z): x stays 0.1, y = 0.7, z = 0.7.
        assert!((other[0] - 0.1).abs() < 1e-10);
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
    fn cubic_point_ops_count() {
        // m-3m point group: 48 operations.
        let ops = cubic_point_ops();
        assert_eq!(ops.len(), 48);
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
    fn p6mmm_orbit_12() {
        let g = SpaceGroup::by_symbol("P6/mmm").unwrap();
        assert_eq!(g.operations.len(), 12);
        let orbit = g.orbit([0.13, 0.29, 0.41], 1e-6);
        assert_eq!(orbit.len(), 12);
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
