//! Cartesian-to-spherical transformation for d shells.
//!
//! The integrals are evaluated over Cartesian components. A spherical d shell
//! uses the five real solid harmonics instead of the six Cartesian `d`
//! functions, and the extra function (an s-like contaminant) must not enter
//! the basis. Each spherical function is a linear combination of the
//! normalized Cartesian components, so every matrix element transforms as
//! `T_aᵀ M T_b` and every packed ERI quartet is contracted along each
//! spherical axis.

use crate::gaussian::{cartesian_components, double_factorial, Shell};
use crate::integrals::overlap_shell;

/// One monomial term of a harmonic: `((lx, ly, lz), coefficient)`.
type Term = ((u8, u8, u8), f64);

/// The five real degree-2 harmonics as monomial combinations. Their overall
/// scale is fixed numerically in [`transform`], so only the relative
/// coefficients matter here.
const D_HARMONICS: [&[Term]; 5] = [
    // 2z² − x² − y²
    &[((0, 0, 2), 2.0), ((2, 0, 0), -1.0), ((0, 2, 0), -1.0)],
    // xz
    &[((1, 0, 1), 1.0)],
    // yz
    &[((0, 1, 1), 1.0)],
    // x² − y²
    &[((2, 0, 0), 1.0), ((0, 2, 0), -1.0)],
    // xy
    &[((1, 1, 0), 1.0)],
];

/// `(2i−1)!!` for one Cartesian exponent `i`, with `(−1)!! = 1`.
fn odd_double_factorial(i: u8) -> f64 {
    if i == 0 {
        1.0
    } else {
        double_factorial(2 * i - 1)
    }
}

/// Matrix `t[c][m]` taking Cartesian components `c` (in
/// [`cartesian_components`] order) to the basis functions `m` of `shell`.
///
/// Cartesian shells map to the identity. Spherical d shells map to five
/// functions normalized to unit self-overlap, using the shell's own
/// Cartesian normalization (`shell.normalizations`).
///
/// Spherical shells with `l ≠ 2` are not supported. `BasisSet` only produces
/// spherical d shells, and for any other `l` a spherical shell is treated as
/// Cartesian.
pub(crate) fn transform(shell: &Shell) -> Vec<Vec<f64>> {
    let comps = cartesian_components(shell.l);
    let n_cart = comps.len();
    if !shell.spherical || shell.l != 2 {
        return (0..n_cart)
            .map(|c| (0..n_cart).map(|m| f64::from(u8::from(c == m))).collect())
            .collect();
    }

    let mut t = vec![vec![0.0f64; D_HARMONICS.len()]; n_cart];
    for (m, terms) in D_HARMONICS.iter().enumerate() {
        for &(mono, coef) in terms.iter() {
            let Some(c) = comps.iter().position(|&x| x == mono) else {
                continue;
            };
            // The Cartesian component is `mono · R(r)` up to the normalization
            // `norm_c / √D_c`, so `coef · mono · R` is `coef · √D_c / norm_c`
            // times the normalized component.
            let (lx, ly, lz) = mono;
            let d = odd_double_factorial(lx) * odd_double_factorial(ly) * odd_double_factorial(lz);
            t[c][m] = coef * d.sqrt() / shell.normalizations[c];
        }
    }

    // Fix the scale from the Cartesian overlap so each function has unit norm.
    let s = overlap_shell(shell, shell);
    for m in 0..D_HARMONICS.len() {
        let mut norm2 = 0.0;
        for c in 0..n_cart {
            for d in 0..n_cart {
                norm2 += t[c][m] * s[c][d] * t[d][m];
            }
        }
        let scale = 1.0 / norm2.sqrt();
        for row in t.iter_mut() {
            row[m] *= scale;
        }
    }
    t
}

/// `Tₐᵀ · block · T_b` for a Cartesian block (`block[c][d]`).
pub(crate) fn transform_block(
    block: &[Vec<f64>],
    ta: &[Vec<f64>],
    tb: &[Vec<f64>],
) -> Vec<Vec<f64>> {
    let na = ta.first().map_or(0, Vec::len);
    let nb = tb.first().map_or(0, Vec::len);
    let mut out = vec![vec![0.0f64; nb]; na];
    for (m, out_row) in out.iter_mut().enumerate() {
        for (n, value) in out_row.iter_mut().enumerate() {
            let mut acc = 0.0;
            for (c, block_row) in block.iter().enumerate() {
                for (d, &v) in block_row.iter().enumerate() {
                    acc += ta[c][m] * v * tb[d][n];
                }
            }
            *value = acc;
        }
    }
    out
}

/// Contract a flat `dims`-shaped tensor along `axis` with `t[c][m]`
/// (old index `c` → new index `m`). Returns the new tensor and its shape.
pub(crate) fn contract_axis(
    v: &[f64],
    dims: [usize; 4],
    axis: usize,
    t: &[Vec<f64>],
) -> (Vec<f64>, [usize; 4]) {
    let n_new = t.first().map_or(0, Vec::len);
    let mut new_dims = dims;
    new_dims[axis] = n_new;
    let stride = |d: [usize; 4]| [d[1] * d[2] * d[3], d[2] * d[3], d[3], 1];
    let s_old = stride(dims);
    let s_new = stride(new_dims);
    let total: usize = new_dims.iter().product();
    let mut out = vec![0.0f64; total];
    for (idx, slot) in out.iter_mut().enumerate() {
        let mut coord = [0usize; 4];
        let mut rem = idx;
        for k in 0..4 {
            coord[k] = rem / s_new[k];
            rem %= s_new[k];
        }
        let m = coord[axis];
        let mut acc = 0.0;
        for (c, row) in t.iter().enumerate().take(dims[axis]) {
            let mut old = 0;
            for k in 0..4 {
                let ck = if k == axis { c } else { coord[k] };
                old += ck * s_old[k];
            }
            acc += row[m] * v[old];
        }
        *slot = acc;
    }
    (out, new_dims)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::cc_pvdz_shells;
    use tpt_chem_core::vec3::Vec3;

    fn d_shell() -> Shell {
        // A single normalized spherical d shell at the origin.
        let mut sh = Shell::new(2, Vec3::ZERO, vec![(1.2, 1.0)]).with_spherical();
        sh.normalize();
        sh
    }

    #[test]
    fn spherical_d_functions_are_orthonormal() {
        let shell = d_shell();
        let t = transform(&shell);
        let s = overlap_shell(&shell, &shell);
        let ts = transform_block(&s, &t, &t);
        for (m, row) in ts.iter().enumerate() {
            for (n, &v) in row.iter().enumerate() {
                let expected = if m == n { 1.0 } else { 0.0 };
                assert!((v - expected).abs() < 1e-10, "S[{m}][{n}] = {v}");
            }
        }
    }

    #[test]
    fn cartesian_shell_is_identity() {
        let mut shell = Shell::new(2, Vec3::ZERO, vec![(1.2, 1.0)]);
        shell.normalize();
        let t = transform(&shell);
        assert_eq!(t.len(), 6);
        assert_eq!(t[0].len(), 6);
        assert_eq!(t[2][2], 1.0);
        assert_eq!(t[2][3], 0.0);
    }

    #[test]
    fn cc_pvdz_oxygen_has_five_d_functions() {
        let shells = cc_pvdz_shells(8, Vec3::ZERO).unwrap();
        let d: Vec<&Shell> = shells.iter().filter(|s| s.l == 2).collect();
        assert_eq!(d.len(), 1);
        assert!(d[0].spherical);
        assert_eq!(transform(d[0])[0].len(), 5);
    }

    #[test]
    fn spherical_d_density_is_rotation_invariant() {
        // The sum of squares over an orthonormal set of one l depends only on
        // |r| (addition theorem). It must not depend on direction, which a
        // wrong set of harmonics would break.
        use crate::cube::ao_values;
        let shells = [d_shell()];
        let radius = 1.7;
        let directions = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.3, -0.5, 0.8).normalize(),
            Vec3::new(-0.6, 0.2, -0.77).normalize(),
        ];
        let sums: Vec<f64> = directions
            .iter()
            .map(|u| {
                let v = ao_values(&shells, *u * radius);
                v.iter().map(|x| x * x).sum()
            })
            .collect();
        for s in &sums[1..] {
            assert!((s - sums[0]).abs() < 1e-12, "{sums:?}");
        }
    }

    #[test]
    fn contract_axis_matches_direct_sum() {
        // 2 × 3 tensor contracted on axis 1 with a 3 × 2 matrix.
        let v: Vec<f64> = (1..=6).map(f64::from).collect();
        let t = vec![vec![1.0, 0.0], vec![0.0, 1.0], vec![1.0, 1.0]];
        let (out, dims) = contract_axis(&v, [2, 3, 1, 1], 1, &t);
        assert_eq!(dims, [2, 2, 1, 1]);
        // Row 0: [1,2,3]→[1+3, 2+3]; row 1: [4,5,6]→[4+6, 5+6].
        assert_eq!(out, vec![4.0, 5.0, 10.0, 11.0]);
    }
}
