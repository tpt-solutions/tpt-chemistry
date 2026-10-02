//! Cartesian Gaussian-type orbitals (GTOs): primitives, shells, and
//! normalization.

use tpt_chem_core::num;
use tpt_chem_core::vec3::Vec3;

/// Angular-momentum components of a Cartesian shell: all `(ix, iy, iz)`
/// with `ix + iy + iz = l`.
pub fn cartesian_components(l: u8) -> Vec<(u8, u8, u8)> {
    let mut out = Vec::new();
    for ix in 0..=l {
        for iy in 0..=(l - ix) {
            out.push((ix, iy, l - ix - iy));
        }
    }
    out
}

/// Number of Cartesian functions in a shell of angular momentum `l`:
/// `(l+1)(l+2)/2`.
pub fn n_cartesian(l: u8) -> usize {
    ((l + 1) * (l + 2) / 2) as usize
}

/// Double factorial `n!!` for the small non-negative `n` used in
/// normalizations (n ≤ 9).
pub fn double_factorial(n: u8) -> f64 {
    match n {
        0 | 1 => 1.0,
        _ => {
            let mut acc = 1.0;
            let mut k = n;
            while k > 1 {
                acc *= f64::from(k);
                k -= 2;
            }
            acc
        }
    }
}

/// Normalization constant of a *primitive* Cartesian Gaussian with
/// exponents `alpha` and angular-momentum components `(lx, ly, lz)`:
///
/// `N = (2α/π)^{3/4} (4α)^{l/2} / √((2lx−1)!!(2ly−1)!!(2lz−1)!!)`
pub fn primitive_normalization(alpha: f64, (lx, ly, lz): (u8, u8, u8)) -> f64 {
    let l = f64::from(lx + ly + lz);
    let pref = num::powf(2.0 * alpha / core::f64::consts::PI, 0.75);
    let ang = num::powf(4.0 * alpha, l / 2.0);
    // (−1)!! = 1 by convention, so s components contribute 1.
    let df = |i: u8| {
        if i == 0 {
            1.0
        } else {
            double_factorial(2 * i - 1)
        }
    };
    let denom = num::sqrt(df(lx) * df(ly) * df(lz));
    pref * ang / denom
}

/// One contracted Gaussian shell at a center.
#[derive(Clone, Debug, PartialEq)]
pub struct Shell {
    /// Angular momentum `l` (0=s, 1=p, 2=d, …).
    pub l: u8,
    /// Center in Bohr.
    pub center: Vec3,
    /// Primitive data: `(exponent, contraction coefficient)` pairs.
    pub primitives: Vec<(f64, f64)>,
    /// Per-Cartesian-component renormalization factors, computed by
    /// [`Shell::normalize`]. Empty until normalized.
    pub normalizations: Vec<f64>,
}

impl Shell {
    /// A shell with the given primitives; contraction coefficients are
    /// *not* renormalized yet — call [`Shell::normalize`].
    pub fn new(l: u8, center: Vec3, primitives: Vec<(f64, f64)>) -> Self {
        let n = n_cartesian(l);
        Shell {
            l,
            center,
            primitives,
            normalizations: vec![1.0; n],
        }
    }

    /// Renormalize each contracted Cartesian component to unit self-overlap.
    pub fn normalize(&mut self) {
        let comps = cartesian_components(self.l);
        let mut norms = Vec::with_capacity(comps.len());
        for &comp in &comps {
            // Same-center self-overlap of the contracted Cartesian
            // component: the 1D factors are Gaussian moments
            // M(2k) = (2k-1)!! sqrt(pi) / (2^k p^(k+1/2)), NOT plain
            // overlaps (using (pi/p)^{3/2} for every component is only
            // correct for s shells and silently breaks all higher shells).
            let mut s = 0.0;
            for &(ai, ci) in &self.primitives {
                for &(aj, cj) in &self.primitives {
                    let ni = primitive_normalization(ai, comp);
                    let nj = primitive_normalization(aj, comp);
                    let p = ai + aj;
                    let mut m = 1.0;
                    for &lx in &[comp.0, comp.1, comp.2] {
                        let k = f64::from(u32::from(lx)); // moment order 2*lx
                                                          // (2lx-1)!! with the (-1)!! = 1 convention for lx = 0.
                        let df = if lx == 0 {
                            1.0
                        } else {
                            double_factorial(2 * lx - 1)
                        };
                        m *= df * num::sqrt(core::f64::consts::PI)
                            / (num::powf(2.0, k) * num::powf(p, k + 0.5));
                    }
                    s += ci * cj * ni * nj * m;
                }
            }
            norms.push(1.0 / num::sqrt(s));
        }
        self.normalizations = norms;
    }
}

/// All Cartesian components `(ix, iy, iz)` with `ix + iy + iz <= l`,
/// ordered by increasing total angular momentum (the intermediate
/// components the OS recurrences walk through).
pub fn all_components_up_to(l: u8) -> Vec<(u8, u8, u8)> {
    let mut out = Vec::new();
    for t in 0..=l {
        for ix in 0..=t {
            for iy in 0..=(t - ix) {
                out.push((ix, iy, t - ix - iy));
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn components_count() {
        assert_eq!(n_cartesian(0), 1);
        assert_eq!(n_cartesian(1), 3);
        assert_eq!(n_cartesian(2), 6);
        assert_eq!(
            cartesian_components(1),
            vec![(0, 0, 1), (0, 1, 0), (1, 0, 0)]
        );
    }

    #[test]
    fn s_normalization_overlap_is_one() {
        // Normalized s primitives at the same center have overlap 1.
        let alpha = 1.24;
        let n = primitive_normalization(alpha, (0, 0, 0));
        let s = n * n * num::powf(core::f64::consts::PI / (2.0 * alpha), 1.5);
        assert!((s - 1.0).abs() < 1e-14);
    }

    #[test]
    fn double_factorial_values() {
        assert_eq!(double_factorial(0), 1.0);
        assert_eq!(double_factorial(1), 1.0);
        assert_eq!(double_factorial(3), 3.0);
        assert_eq!(double_factorial(5), 15.0);
        assert_eq!(double_factorial(7), 105.0);
    }

    #[test]
    fn contracted_shell_is_normalized() {
        // A 3-term contracted s shell must have unit self-overlap.
        let mut shell = Shell::new(
            0,
            Vec3::ZERO,
            vec![
                (3.42525091, 0.15432897),
                (0.62391373, 0.53532814),
                (0.16885540, 0.44463454),
            ],
        );
        shell.normalize();
        let [n] = [shell.normalizations[0]];
        let mut s = 0.0;
        for &(ai, ci) in &shell.primitives {
            for &(aj, cj) in &shell.primitives {
                let ni = primitive_normalization(ai, (0, 0, 0));
                let nj = primitive_normalization(aj, (0, 0, 0));
                let p = ai + aj;
                s += ci * cj * ni * nj * num::powf(core::f64::consts::PI / p, 1.5);
            }
        }
        assert!((n * n * s - 1.0).abs() < 1e-12);
    }
}
