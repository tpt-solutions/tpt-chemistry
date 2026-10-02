//! From-scratch one- and two-electron integral evaluation over contracted
//! Cartesian Gaussians.
//!
//! All routines are primitive-first: the recurrences run per primitive
//! combination and the shell-level wrappers contract with the coefficients
//! and component normalizations.
//!
//! Recursion families:
//!
//! * **Overlap** — per-axis 1D Obara–Saika recurrence (Obara & Saika,
//!   J. Chem. Phys. 84, 3963 (1986)).
//! * **Kinetic** — derivatives of the overlap (`⟨a|−½∇²|b⟩`).
//! * **Nuclear attraction** — McMurchie–Davidson: the primitive pair is
//!   expanded in Hermite Gaussians about `P` (per-axis E-coefficients) and
//!   contracted against the Hermite Coulomb auxiliary `R^n_{tuv}`.
//! * **Electron repulsion** — McMurchie–Davidson (McMurchie & Davidson,
//!   J. Chem. Phys. 68, 3344 (1978); Helgaker, Jørgensen & Olsen,
//!   "Molecular Electronic-Structure Theory", ch. 9): E-coefficients of
//!   both primitive pairs contracted against the two-electron Hermite
//!   Coulomb auxiliary `R^n_{tuv}(ρ, P−Q)` at the reduced exponent
//!   `ρ = pq/(p+q)`, with the electron-2 Hermite indices carrying the
//!   `(−1)^{τ+ν+φ}` sign.
//!
//! Atomic units throughout.

use tpt_chem_core::num;
use tpt_chem_core::vec3::Vec3;

use crate::boys::boys_array;
use crate::gaussian::{cartesian_components, primitive_normalization, Shell};

/// One primitive of a contracted shell (the component-dependent Cartesian
/// normalization is applied at the shell wrappers).
#[derive(Clone, Copy)]
struct Prim {
    alpha: f64,
    coeff: f64,
}

/// Precomputed two-center pair quantities.
struct Pair {
    p: f64,
    /// The two primitive exponents.
    alpha: f64,
    beta: f64,
    /// Product-Gaussian center `P = (αA + βB)/p`.
    center: Vec3,
    /// `P − A`.
    pa: Vec3,
    /// `A − B` (horizontal transfer vector).
    ab: Vec3,
}

fn make_pair(alpha: f64, a: Vec3, beta: f64, b: Vec3) -> Pair {
    let p = alpha + beta;
    let ab = a - b;
    let center = (a * alpha + b * beta) / p;
    Pair {
        p,
        alpha,
        beta,
        center,
        pa: center - a,
        ab,
    }
}

fn prims(shell: &Shell) -> Vec<Prim> {
    shell
        .primitives
        .iter()
        .map(|&(alpha, coeff)| Prim { alpha, coeff })
        .collect()
}

/// Per-axis primitive overlap tables. Axis table `s[i][j]` covers
/// `i ≤ la+1` (the extra row feeds horizontal transfers) and `j ≤ jmax`,
/// with `s[0][0] = exp(−αβ/p·Δ²)` on that axis.
#[allow(clippy::needless_range_loop)]
fn overlap_tables(alpha: f64, beta: f64, pair: &Pair, la: u8, jmax: u8) -> [Vec<Vec<f64>>; 3] {
    let mut tables: [Vec<Vec<f64>>; 3] = Default::default();
    for axis in 0..3 {
        let n_i = la as usize + 2;
        let n_j = jmax as usize + 1;
        let mut s = vec![vec![0.0; n_j]; n_i];
        let pa = pair.pa.get(axis);
        let ab = pair.ab.get(axis);
        let delta = pair.ab.get(axis);
        s[0][0] = num::exp(-alpha * beta / pair.p * delta * delta);
        let inv2p = 0.5 / pair.p;
        // Iterate by increasing total i+j; within a total, raising cases
        // (i >= 1) first, then horizontal transfers (i = 0), because
        // (0,j) reads (1,j-1) at the same total.
        for total in 1..=(n_i + n_j) {
            for i in 1..n_i {
                for j in 0..n_j {
                    if i + j != total {
                        continue;
                    }
                    let mut v = pa * s[i - 1][j];
                    if i > 1 {
                        v += (i as f64 - 1.0) * inv2p * s[i - 2][j];
                    }
                    if j > 0 {
                        v += j as f64 * inv2p * s[i - 1][j - 1];
                    }
                    s[i][j] = v;
                }
            }
            for j in 1..n_j {
                if j != total {
                    continue;
                }
                // Horizontal: (0,j) = (1,j-1) - AB*(0,j-1).
                s[0][j] = s[1][j - 1] + ab * s[0][j - 1];
            }
        }
        tables[axis] = s;
    }
    tables
}

/// Read a 3D overlap component from the per-axis tables; 0 when any
/// angular index is out of range (a vanishing derivative term).
fn table_val(tables: &[Vec<Vec<f64>>; 3], a: (u8, u8, u8), b: (u8, u8, u8)) -> f64 {
    let (ix, iy, iz) = a;
    let (jx, jy, jz) = b;
    let t = |axis: usize, i: u8, j: u8| {
        tables[axis]
            .get(i as usize)
            .and_then(|row| row.get(j as usize))
            .copied()
            .unwrap_or(0.0)
    };
    t(0, ix, jx) * t(1, iy, jy) * t(2, iz, jz)
}

/// Shell-pair overlap matrix `[component of a][component of b]`,
/// contracted and normalized.
pub fn overlap_shell(a: &Shell, b: &Shell) -> Vec<Vec<f64>> {
    let ca = cartesian_components(a.l);
    let cb = cartesian_components(b.l);
    let mut out = vec![vec![0.0; cb.len()]; ca.len()];
    for pi in prims(a) {
        for pj in prims(b) {
            let pair = make_pair(pi.alpha, a.center, pj.alpha, b.center);
            // The per-axis table bases already carry the exp(-ab/p*d^2) factors
            // (product = pair.kab); the global factor is only (pi/p)^{3/2}.
            let global = num::powf(core::f64::consts::PI / pair.p, 1.5);
            let tables = overlap_tables(pi.alpha, pj.alpha, &pair, a.l, b.l);
            for (ia, &comp_a) in ca.iter().enumerate() {
                let na = primitive_normalization(pi.alpha, comp_a) * a.normalizations[ia];
                for (ib, &comp_b) in cb.iter().enumerate() {
                    let nb = primitive_normalization(pj.alpha, comp_b) * b.normalizations[ib];
                    out[ia][ib] +=
                        pi.coeff * pj.coeff * na * nb * global * table_val(&tables, comp_a, comp_b);
                }
            }
        }
    }
    out
}

/// Shell-pair kinetic-energy matrix `⟨a|−½∇²|b⟩` from derivatives of the
/// primitive overlap:
/// `T = −½ Σ_axis [ j(j−1)·S(a,b−2) − 2β(2j+1)·S(a,b) + 4β²·S(a,b+2) ]`.
pub fn kinetic_shell(a: &Shell, b: &Shell) -> Vec<Vec<f64>> {
    let ca = cartesian_components(a.l);
    let cb = cartesian_components(b.l);
    let mut out = vec![vec![0.0; cb.len()]; ca.len()];
    for pi in prims(a) {
        for pj in prims(b) {
            let beta = pj.alpha;
            let pair = make_pair(pi.alpha, a.center, beta, b.center);
            // The per-axis table bases already carry the exp(-ab/p*d^2) factors
            // (product = pair.kab); the global factor is only (pi/p)^{3/2}.
            let global = num::powf(core::f64::consts::PI / pair.p, 1.5);
            let tables = overlap_tables(pi.alpha, beta, &pair, a.l, b.l + 2);
            for (ia, &comp_a) in ca.iter().enumerate() {
                let na = primitive_normalization(pi.alpha, comp_a) * a.normalizations[ia];
                for (ib, &comp_b) in cb.iter().enumerate() {
                    let nb = primitive_normalization(pj.alpha, comp_b) * b.normalizations[ib];
                    let (jx, jy, jz) = comp_b;
                    let mut t = 0.0;
                    for axis in 0..3 {
                        let j = [jx, jy, jz][axis] as f64;
                        let shifted = |dj: i32| -> f64 {
                            let jj = [jx as i32, jy as i32, jz as i32][axis] + dj;
                            if jj < 0 {
                                return 0.0;
                            }
                            let mut comp = comp_b;
                            comp.0 = if axis == 0 { jj as u8 } else { jx };
                            comp.1 = if axis == 1 { jj as u8 } else { jy };
                            comp.2 = if axis == 2 { jj as u8 } else { jz };
                            global * table_val(&tables, comp_a, comp)
                        };
                        t += j * (j - 1.0) * shifted(-2)
                            - 2.0 * beta * (2.0 * j + 1.0) * shifted(0)
                            + 4.0 * beta * beta * shifted(2);
                    }
                    out[ia][ib] += pi.coeff * pj.coeff * na * nb * (-0.5) * t;
                }
            }
        }
    }
    out
}

/// Primitive ERI `(ab|cd)` for all component combinations, via the
/// McMurchie–Davidson Hermite expansion: per-axis E-coefficients for each
/// primitive pair (seeded with the per-axis `exp(−αβ/p·Δ²)` overlap
/// factor), contracted against the two-electron Hermite Coulomb auxiliary
/// `R^n_{tuv}(ρ, P−Q)` at the reduced exponent `ρ = pq/(p+q)`; the
/// electron-2 Hermite indices carry the `(−1)^{τ+ν+φ}` sign.
#[allow(clippy::too_many_lines)]
fn eri_primitive(
    pair1: &Pair,
    pair2: &Pair,
    la: u8,
    lb: u8,
    lc: u8,
    ld: u8,
) -> Vec<Vec<Vec<Vec<f64>>>> {
    let p = pair1.p;
    let q = pair2.p;
    let prefac = 2.0 * num::powf(core::f64::consts::PI, 2.5) / (p * q * num::sqrt(p + q));
    let mut e1: [Vec<Vec<Vec<f64>>>; 3] = Default::default();
    let mut e2: [Vec<Vec<Vec<f64>>>; 3] = Default::default();
    for axis in 0..3 {
        let d1 = pair1.ab.get(axis);
        let k1 = num::exp(-pair1.alpha * pair1.beta / p * d1 * d1);
        e1[axis] = hermite_coeffs_1d(
            pair1.pa.get(axis),
            (pair1.pa + pair1.ab).get(axis),
            p,
            la as usize,
            lb as usize,
            k1,
        );
        let d2 = pair2.ab.get(axis);
        let k2 = num::exp(-pair2.alpha * pair2.beta / q * d2 * d2);
        e2[axis] = hermite_coeffs_1d(
            pair2.pa.get(axis),
            (pair2.pa + pair2.ab).get(axis),
            q,
            lc as usize,
            ld as usize,
            k2,
        );
    }
    let rho = p * q / (p + q);
    let pq = pair1.center - pair2.center;
    let r_max = (la + lb + lc + ld) as usize;
    let r = hermite_coulomb_r(rho, pq, r_max, r_max, r_max, r_max);

    let a_comps = cartesian_components(la);
    let b_comps = cartesian_components(lb);
    let c_comps = cartesian_components(lc);
    let d_comps = cartesian_components(ld);
    let mut out =
        vec![vec![vec![vec![0.0f64; d_comps.len()]; c_comps.len()]; b_comps.len()]; a_comps.len()];
    for (ia, &a) in a_comps.iter().enumerate() {
        for (ib, &b) in b_comps.iter().enumerate() {
            for (ic, &c) in c_comps.iter().enumerate() {
                for (id, &d) in d_comps.iter().enumerate() {
                    let mut sum = 0.0;
                    for t in 0..=(a.0 + b.0) {
                        for u in 0..=(a.1 + b.1) {
                            for v in 0..=(a.2 + b.2) {
                                let ea = e1[0][a.0 as usize][b.0 as usize][t as usize]
                                    * e1[1][a.1 as usize][b.1 as usize][u as usize]
                                    * e1[2][a.2 as usize][b.2 as usize][v as usize];
                                for tau in 0..=(c.0 + d.0) {
                                    for nu in 0..=(c.1 + d.1) {
                                        for phi in 0..=(c.2 + d.2) {
                                            let sign =
                                                if (tau + nu + phi) % 2 == 0 { 1.0 } else { -1.0 };
                                            let ec = e2[0][c.0 as usize][d.0 as usize]
                                                [tau as usize]
                                                * e2[1][c.1 as usize][d.1 as usize][nu as usize]
                                                * e2[2][c.2 as usize][d.2 as usize][phi as usize];
                                            sum += sign
                                                * ea
                                                * ec
                                                * r[(t + tau) as usize][(u + nu) as usize]
                                                    [(v + phi) as usize][0];
                                        }
                                    }
                                }
                            }
                        }
                    }
                    out[ia][ib][ic][id] = prefac * sum;
                }
            }
        }
    }
    out
}

/// Shell-pair nuclear-attraction matrix `[a][b]` for
/// `⟨a|Σ_C −Z_C/r_C|b⟩`, contracted and normalized. `nuclei` is a list of
/// `(charge, position in Bohr)`.
pub fn nuclear_shell(a: &Shell, b: &Shell, nuclei: &[(f64, Vec3)]) -> Vec<Vec<f64>> {
    let ca = cartesian_components(a.l);
    let cb = cartesian_components(b.l);
    let mut out = vec![vec![0.0; cb.len()]; ca.len()];
    for pi in prims(a) {
        for pj in prims(b) {
            let pair = make_pair(pi.alpha, a.center, pj.alpha, b.center);
            for (ia, &comp_a) in ca.iter().enumerate() {
                let na = primitive_normalization(pi.alpha, comp_a) * a.normalizations[ia];
                let fa = pi.coeff * na;
                for (ib, &comp_b) in cb.iter().enumerate() {
                    let nb = primitive_normalization(pj.alpha, comp_b) * b.normalizations[ib];
                    let fb = pj.coeff * nb;
                    let mut sum = 0.0;
                    for &(z, c) in nuclei {
                        // nuclear_value already carries the -Z convention
                        // (the -(2pi/p) prefactor scales with charge).
                        sum += z * nuclear_value(&pair, c, a.l, b.l, comp_a, comp_b);
                    }
                    out[ia][ib] += fa * fb * sum;
                }
            }
        }
    }
    out
}

/// McMurchie-Davidson 1D Hermite expansion coefficients for the pair
/// `(x-A)^i (x-B)^j exp(-alpha (x-A)^2 - beta (x-B)^2)` about `P`:
/// returned as `E[i][j][t]`. Seeded with `E_0^{00} = kab_axis`.
fn hermite_coeffs_1d(
    pa: f64,
    pb: f64,
    p: f64,
    i_max: usize,
    j_max: usize,
    kab_axis: f64,
) -> Vec<Vec<Vec<f64>>> {
    let t_max = i_max + j_max;
    let mut e = vec![vec![vec![0.0f64; t_max + 2]; j_max + 1]; i_max + 1];
    e[0][0][0] = kab_axis;
    let inv2p = 0.5 / p;
    for i in 0..i_max {
        for t in 0..=t_max {
            let mut v = inv2p * if t >= 1 { e[i][0][t - 1] } else { 0.0 } + pa * e[i][0][t];
            if t + 1 < e[i][0].len() {
                v += (t + 1) as f64 * e[i][0][t + 1];
            }
            e[i + 1][0][t] = v;
        }
    }
    #[allow(clippy::needless_range_loop)]
    for i in 0..=i_max {
        for j in 0..j_max {
            for t in 0..=t_max {
                let mut v = inv2p * if t >= 1 { e[i][j][t - 1] } else { 0.0 } + pb * e[i][j][t];
                if t + 1 < e[i][j].len() {
                    v += (t + 1) as f64 * e[i][j][t + 1];
                }
                e[i][j + 1][t] = v;
            }
        }
    }
    e
}

/// Hermite Coulomb integrals `R[t][u][v][n]`:
/// `R_{t+1,u,v}^n = t R_{t-1,u,v}^{n+1} + (P-C)_x R_{tuv}^{n+1}` (and
/// symmetrically for u, v), seeded with
/// `R_{000}^n = (-2p)^n F_n(T)`, `T = p|P-C|^2`.
fn hermite_coulomb_r(
    p: f64,
    pc: Vec3,
    t_max: usize,
    u_max: usize,
    v_max: usize,
    max_m: usize,
) -> Vec<Vec<Vec<Vec<f64>>>> {
    let mut r = vec![vec![vec![vec![0.0f64; max_m + 2]; v_max + 1]; u_max + 1]; t_max + 1];
    let boys = boys_array(p * pc.norm_sq(), max_m);
    let mut sign = 1.0;
    for n in 0..=max_m {
        r[0][0][0][n] = sign * boys[n];
        sign *= -2.0 * p;
    }
    for total in 1..=(t_max + u_max + v_max) {
        for t in 0..=t_max.min(total) {
            for u in 0..=u_max.min(total - t) {
                for v in 0..=v_max.min(total - t - u) {
                    if t + u + v != total {
                        continue;
                    }
                    for n in 0..=max_m {
                        if n + 1 > max_m {
                            break;
                        }
                        let mut val = 0.0;
                        if t >= 1 {
                            val += pc.x * r[t - 1][u][v][n + 1];
                            if t >= 2 {
                                val += (t - 1) as f64 * r[t - 2][u][v][n + 1];
                            }
                        } else if u >= 1 {
                            val += pc.y * r[0][u - 1][v][n + 1];
                            if u >= 2 {
                                val += (u - 1) as f64 * r[0][u - 2][v][n + 1];
                            }
                        } else {
                            val += pc.z * r[0][0][v - 1][n + 1];
                            if v >= 2 {
                                val += (v - 1) as f64 * r[0][0][v - 2][n + 1];
                            }
                        }
                        r[t][u][v][n] = val;
                    }
                }
            }
        }
    }
    r
}

/// `m = 0` nuclear-attraction integral for one component pair via
/// McMurchie-Davidson. Carries the `-Z` attraction convention (negative for
/// positive `Z`); the caller multiplies by the charge.
#[allow(clippy::only_used_in_recursion)]
fn nuclear_value(
    pair: &Pair,
    c: Vec3,
    la: u8,
    lb: u8,
    comp_a: (u8, u8, u8),
    comp_b: (u8, u8, u8),
) -> f64 {
    // The i-raise path is the validated one; canonicalize so the
    // higher-momentum side is always the a-side. (The recursion is
    // depth-1: after the flip the higher side is always a-side.)
    let ta = u16::from(comp_a.0) + u16::from(comp_a.1) + u16::from(comp_a.2);
    let tb = u16::from(comp_b.0) + u16::from(comp_b.1) + u16::from(comp_b.2);
    if tb > ta {
        #[allow(clippy::only_used_in_recursion)]
        return {
            let a_center = pair.center - pair.pa;
            let b_center = a_center - pair.ab;
            let flipped = make_pair(pair.beta, b_center, pair.alpha, a_center);
            nuclear_value(&flipped, c, lb, la, comp_b, comp_a)
        };
    }
    let p = pair.p;
    let pa = pair.pa; // P - A
    let pb = pa + pair.ab; // P-B = (P-A) + (A-B)
    let i_ax = [comp_a.0 as usize, comp_a.1 as usize, comp_a.2 as usize];
    let j_ax = [comp_b.0 as usize, comp_b.1 as usize, comp_b.2 as usize];
    let t_max = i_ax[0] + j_ax[0] + i_ax[1] + j_ax[1] + i_ax[2] + j_ax[2];
    let mut axes_e: [Vec<Vec<Vec<f64>>>; 3] = Default::default();
    for axis in 0..3 {
        let delta = pair.ab.get(axis);
        let kab_axis = num::exp(-pair.alpha * pair.beta / p * delta * delta);
        axes_e[axis] = hermite_coeffs_1d(
            pa.get(axis),
            pb.get(axis),
            p,
            i_ax[axis],
            j_ax[axis],
            kab_axis,
        );
    }
    let pc = pair.center - c;
    let r = hermite_coulomb_r(p, pc, t_max, t_max, t_max, t_max);
    let tx = i_ax[0] + j_ax[0];
    let uy = i_ax[1] + j_ax[1];
    let vz = i_ax[2] + j_ax[2];
    let mut total = 0.0;
    for t in 0..=tx {
        let ex = axes_e[0][i_ax[0]][j_ax[0]][t];
        for u in 0..=uy {
            let ey = axes_e[1][i_ax[1]][j_ax[1]][u];
            for v in 0..=vz {
                let ezv = axes_e[2][i_ax[2]][j_ax[2]][v];
                total += ex * ey * ezv * r[t][u][v][0];
            }
        }
    }
    -(2.0 * core::f64::consts::PI / p) * total
}

/// Shell-quartet ERI tensor `[a][b][c][d]`, contracted and normalized.
pub fn eri_shell(a: &Shell, b: &Shell, c: &Shell, d: &Shell) -> Vec<Vec<Vec<Vec<f64>>>> {
    let ca = cartesian_components(a.l);
    let cb = cartesian_components(b.l);
    let cc = cartesian_components(c.l);
    let cd = cartesian_components(d.l);
    let mut out = vec![vec![vec![vec![0.0f64; cd.len()]; cc.len()]; cb.len()]; ca.len()];
    for pi in prims(a) {
        for pj in prims(b) {
            let pair1 = make_pair(pi.alpha, a.center, pj.alpha, b.center);
            for pk in prims(c) {
                for pl in prims(d) {
                    let pair2 = make_pair(pk.alpha, c.center, pl.alpha, d.center);
                    let v = eri_primitive(&pair1, &pair2, a.l, b.l, c.l, d.l);
                    for (ia, &comp_a) in ca.iter().enumerate() {
                        let na = primitive_normalization(pi.alpha, comp_a) * a.normalizations[ia];
                        for (ib, &comp_b) in cb.iter().enumerate() {
                            let nb =
                                primitive_normalization(pj.alpha, comp_b) * b.normalizations[ib];
                            for (ic, &comp_c) in cc.iter().enumerate() {
                                let nc = primitive_normalization(pk.alpha, comp_c)
                                    * c.normalizations[ic];
                                for (id, &comp_d) in cd.iter().enumerate() {
                                    let nd = primitive_normalization(pl.alpha, comp_d)
                                        * d.normalizations[id];
                                    out[ia][ib][ic][id] += pi.coeff
                                        * pj.coeff
                                        * pk.coeff
                                        * pl.coeff
                                        * na
                                        * nb
                                        * nc
                                        * nd
                                        * v[ia][ib][ic][id];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    out
}

/// Nuclear repulsion energy for `(charge, position in Bohr)` nuclei.
pub fn nuclear_repulsion(nuclei: &[(f64, Vec3)]) -> f64 {
    let mut e = 0.0;
    for i in 0..nuclei.len() {
        for j in (i + 1)..nuclei.len() {
            e += nuclei[i].0 * nuclei[j].0 / nuclei[i].1.dist(nuclei[j].1);
        }
    }
    e
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gaussian::Shell;

    fn s_shell(alpha: f64, center: Vec3) -> Shell {
        let mut sh = Shell::new(0, center, vec![(alpha, 1.0)]);
        sh.normalize();
        sh
    }

    #[test]
    fn overlap_same_center_is_one() {
        let a = s_shell(1.24, Vec3::ZERO);
        let s = overlap_shell(&a, &a);
        assert!((s[0][0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn overlap_ss_analytic() {
        // Normalized 1s pair: S = exp(−αβ/(α+β)·R²).
        let (alpha, r) = (1.0, 0.8);
        let a = s_shell(alpha, Vec3::new(-r / 2.0, 0.0, 0.0));
        let b = s_shell(alpha, Vec3::new(r / 2.0, 0.0, 0.0));
        let want = num::exp(-alpha / 2.0 * r * r);
        let got = overlap_shell(&a, &b)[0][0];
        assert!((got - want).abs() < 1e-12, "{got} vs {want}");
    }

    #[test]
    fn kinetic_same_center() {
        // ⟨a|−½∇²|a⟩ for a normalized 1s with exponent α is α/2.
        let alpha = 1.7;
        let a = s_shell(alpha, Vec3::new(0.3, -0.2, 0.5));
        let t = kinetic_shell(&a, &a)[0][0];
        assert!((t - 1.5 * alpha).abs() < 1e-12, "T = {t}");
    }

    #[test]
    fn nuclear_attraction_1s_analytic() {
        // ⟨g| r_C^{−1} |g⟩ with g a normalized 1s at C is 2√(α/π).
        let alpha = 1.24;
        let c = Vec3::new(0.4, 0.1, -0.3);
        let a = s_shell(alpha, c);
        let v = nuclear_shell(&a, &a, &[(1.0, c)])[0][0];
        let want = -2.0 * num::sqrt(2.0 * alpha / core::f64::consts::PI);
        assert!((v - want).abs() < 1e-10, "{v} vs {want}");
    }

    #[test]
    fn eri_ss_self_analytic() {
        // Four co-located normalized 1s Gaussians: (ss|ss) = 5α/8.
        let alpha = 0.9;
        let a = s_shell(alpha, Vec3::new(0.2, 0.1, 0.0));
        let v = eri_shell(&a, &a, &a, &a)[0][0][0][0];
        let want = 2.0 * num::sqrt(alpha / core::f64::consts::PI);
        assert!((v - want).abs() < 1e-10, "{v} vs {want}");
    }

    #[test]
    fn eri_symmetries() {
        // (12|34) = (21|34) = (12|43) = (34|12).
        let alpha = 1.3;
        let a = s_shell(alpha, Vec3::new(0.0, 0.0, 0.0));
        let b = s_shell(alpha, Vec3::new(0.5, 0.0, 0.0));
        let c = s_shell(alpha, Vec3::new(0.0, 0.7, 0.0));
        let d = s_shell(alpha, Vec3::new(0.2, 0.0, 0.9));
        let abcd = eri_shell(&a, &b, &c, &d)[0][0][0][0];
        let bacd = eri_shell(&b, &a, &c, &d)[0][0][0][0];
        let abdc = eri_shell(&a, &b, &d, &c)[0][0][0][0];
        let cdab = eri_shell(&c, &d, &a, &b)[0][0][0][0];
        assert!((abcd - bacd).abs() < 1e-12);
        assert!((abcd - abdc).abs() < 1e-12);
        assert!((abcd - cdab).abs() < 1e-12);
    }

    #[test]
    fn nuclear_repulsion_simple() {
        // Two protons 1 Bohr apart: 1 Hartree.
        let e = nuclear_repulsion(&[(1.0, Vec3::ZERO), (1.0, Vec3::new(1.0, 0.0, 0.0))]);
        assert!((e - 1.0).abs() < 1e-12);
    }
}

#[cfg(test)]
mod p_shell_tests {
    use super::*;
    use crate::gaussian::Shell;

    fn pz(alpha: f64, center: Vec3) -> Shell {
        // p_z: the (0, 0, 1) component is index 0 of l=1 components.
        let mut sh = Shell::new(1, center, vec![(alpha, 1.0)]);
        sh.normalize();
        sh
    }

    fn s_shell(alpha: f64, center: Vec3) -> Shell {
        let mut sh = Shell::new(0, center, vec![(alpha, 1.0)]);
        sh.normalize();
        sh
    }

    #[test]
    fn p_self_kinetic_and_nuclear() {
        // Closed forms (verified by numerical quadrature):
        // ⟨2p|T|2p⟩ = 5α/2;  ⟨2p|1/r|2p⟩ = (2/3)√(2α/π).
        let alpha = 1.3;
        let c = Vec3::new(0.2, -0.4, 0.6);
        let sh = pz(alpha, c);
        let t = kinetic_shell(&sh, &sh)[0][0];
        assert!((t - 2.5 * alpha).abs() < 1e-10, "T = {t}");
        let v = nuclear_shell(&sh, &sh, &[(1.0, c)])[0][0];
        let want = -(4.0 / 3.0) * num::sqrt(2.0 * alpha / core::f64::consts::PI);
        assert!((v - want).abs() < 1e-10, "V = {v} vs {want}");
    }

    #[test]
    fn p_s_overlap_analytic() {
        // ⟨p_z^A|s^B⟩ = n_p n_s Kab √(π/p) (P−A)_z.
        let (alpha, beta, r) = (1.1, 0.9, 0.9);
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(0.0, 0.0, r);
        let pa_sh = pz(alpha, a);
        let mut sb = Shell::new(0, b, vec![(beta, 1.0)]);
        sb.normalize();
        let got = overlap_shell(&pa_sh, &sb)[0][0];
        let p = alpha + beta;
        let center_z = (beta * r) / p; // P_z − A_z
        let kab = num::exp(-alpha * beta / p * r * r);
        let want = primitive_normalization(alpha, (0, 0, 1))
            * primitive_normalization(beta, (0, 0, 0))
            * kab
            * num::powf(core::f64::consts::PI / p, 1.5)
            * center_z;
        assert!((got - want).abs() < 1e-12, "{got} vs {want}");
    }

    #[test]
    fn p_p_overlap_sigma_analytic() {
        // ⟨p_z^A|p_z^B⟩ = n²Kab√(π/p)[1/(2p) + pa_z·pb_z].
        let (alpha, r) = (1.4, 1.0);
        let a = Vec3::new(0.0, 0.0, -r / 2.0);
        let b = Vec3::new(0.0, 0.0, r / 2.0);
        let sa = pz(alpha, a);
        let sb = pz(alpha, b);
        let got = overlap_shell(&sa, &sb)[0][0];
        let p = 2.0 * alpha;
        let pa_z = r / 2.0;
        let pb_z = -r / 2.0;
        let kab = num::exp(-alpha * alpha / p * r * r);
        let n2 = primitive_normalization(alpha, (0, 0, 1)).powi(2);
        let want =
            n2 * kab * num::powf(core::f64::consts::PI / p, 1.5) * (1.0 / (2.0 * p) + pa_z * pb_z);
        assert!((got - want).abs() < 1e-12, "{got} vs {want}");
    }

    #[test]
    fn eri_pp_ss_self_analytic() {
        // Co-located normalized 2p/1s primitives with exponent α:
        // (p_z p_z|ss) = (5/3)√(α/π). Derived from
        // (1/3)(−∂/∂p)·[2π^{5/2}/(pq√(p+q))] at p = q = 2α.
        let alpha = 0.9;
        let c = Vec3::new(0.2, -0.4, 0.6);
        let p = pz(alpha, c);
        let s = s_shell(alpha, c);
        let v = eri_shell(&p, &p, &s, &s)[0][0][0][0];
        let want = 5.0 / 3.0 * num::sqrt(alpha / core::f64::consts::PI);
        assert!((v - want).abs() < 1e-10, "{v} vs {want}");
    }

    #[test]
    fn eri_pp_pp_self_tensor() {
        // Co-located equal-exponent p shells: isotropy fixes the full
        // (p_μ p_ν|p_λ p_σ) tensor to
        //   T = a·δ_μν δ_λσ + b·(δ_μλ δ_νσ + δ_μσ δ_νλ),
        // with a, b cross-validated against 3D Gauss–Hermite quadrature of
        // ∫ρ_ab(r1) V_cd(r1) dr1 (α = 0.9):
        //   (p_z p_z|p_z p_z) = a + 2b,  (p_z p_z|p_x p_x) = a,
        //   (p_z p_x|p_z p_x) = (p_z p_x|p_x p_z) = b,
        // and every mismatched-index element vanishes.
        let alpha = 0.9;
        let a_val = 0.767_173_295_3;
        let b_val = 0.053_523_721_9;
        let c = Vec3::new(0.2, -0.4, 0.6);
        let p = pz(alpha, c);
        let v = eri_shell(&p, &p, &p, &p);
        for (i, row) in v.iter().enumerate() {
            for (j, r2) in row.iter().enumerate() {
                for (k, r3) in r2.iter().enumerate() {
                    for (l, &val) in r3.iter().enumerate() {
                        // Map component indices to exponent vectors.
                        let comp = |idx: usize| [(0u8, 0u8, 1u8), (0, 1, 0), (1, 0, 0)][idx.min(2)];
                        let (m, n, lam, sig) = (comp(i), comp(j), comp(k), comp(l));
                        let want = if m == n && lam == sig {
                            if m == lam {
                                a_val + 2.0 * b_val
                            } else {
                                a_val
                            }
                        } else if (m == lam && n == sig) || (m == sig && n == lam) {
                            b_val
                        } else {
                            0.0
                        };
                        assert!(
                            (val - want).abs() < 1e-7,
                            "v[{i}][{j}][{k}][{l}] = {val} vs {want}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn eri_ss_four_center_primitive_formula() {
        // Four normalized single-primitive s shells at distinct centers:
        // compare against the textbook closed form
        // (ss|ss) = 2π^{5/2}/(pq√(p+q))·K_ab·K_cd·F_0(ρ|P−Q|²),
        // ρ = pq/(p+q) — this pins the reduced exponent and the K factors
        // at nonzero Boys argument. Single-primitive s shells carry
        // per-function normalization (2α/π)^{3/4} in `eri_shell`.
        let alphas = [1.1, 0.7, 1.3, 0.5];
        let centers = [
            Vec3::new(0.1, 0.2, 0.3),
            Vec3::new(0.9, 0.0, 0.0),
            Vec3::new(0.0, 0.7, 0.4),
            Vec3::new(0.3, 0.3, 1.1),
        ];
        let shells: Vec<Shell> = centers
            .iter()
            .zip(alphas)
            .map(|(&c, a)| {
                let mut sh = Shell::new(0, c, vec![(a, 1.0)]);
                sh.normalize();
                sh
            })
            .collect();
        let got = eri_shell(&shells[0], &shells[1], &shells[2], &shells[3])[0][0][0][0];
        let (pa, pb, pc, pd) = (alphas[0], alphas[1], alphas[2], alphas[3]);
        let p = pa + pb;
        let q = pc + pd;
        let rho = p * q / (p + q);
        let ab = centers[0] - centers[1];
        let cd = centers[2] - centers[3];
        let pq = (centers[0] * pa + centers[1] * pb) / p - (centers[2] * pc + centers[3] * pd) / q;
        let boys = crate::boys::boys_array(rho * pq.norm_sq(), 0);
        let norm_prod: f64 = alphas
            .iter()
            .map(|&a| primitive_normalization(a, (0, 0, 0)))
            .product();
        let want = 2.0 * num::powf(core::f64::consts::PI, 2.5) / (p * q * num::sqrt(p + q))
            * num::exp(-pa * pb / p * ab.norm_sq())
            * num::exp(-pc * pd / q * cd.norm_sq())
            * boys[0]
            * norm_prod;
        assert!((got - want).abs() < 1e-12, "{got} vs {want}");
    }

    #[test]
    fn eri_p_shell_symmetries() {
        // Permutation identities in tensor-index form: swapping functions
        // 1↔2 transposes the first two index pair, 3↔4 the last pair, and
        // the electron exchange (12|34) = (34|12) swaps the index pairs.
        let alpha = 1.3;
        let pz_at = |c: Vec3| pz(alpha, c);
        let a = pz_at(Vec3::new(0.0, 0.0, 0.0));
        let b = pz_at(Vec3::new(0.5, 0.0, 0.0));
        let c = pz_at(Vec3::new(0.0, 0.7, 0.0));
        let d = pz_at(Vec3::new(0.2, 0.0, 0.9));
        let abcd = eri_shell(&a, &b, &c, &d);
        let bacd = eri_shell(&b, &a, &c, &d);
        let abdc = eri_shell(&a, &b, &d, &c);
        let cdab = eri_shell(&c, &d, &a, &b);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    for l in 0..3 {
                        let v = abcd[i][j][k][l];
                        assert!(
                            (v - bacd[j][i][k][l]).abs() < 1e-12,
                            "(12) swap at [{i}][{j}][{k}][{l}]"
                        );
                        assert!(
                            (v - abdc[i][j][l][k]).abs() < 1e-12,
                            "(34) swap at [{i}][{j}][{k}][{l}]"
                        );
                        assert!(
                            (v - cdab[k][l][i][j]).abs() < 1e-12,
                            "(13|24) swap at [{i}][{j}][{k}][{l}]"
                        );
                    }
                }
            }
        }
    }
}
