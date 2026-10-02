//! From-scratch one- and two-electron integral evaluation over contracted
//! Cartesian Gaussians, via the Obara–Saika (OS) recursion scheme
//! (Obara & Saika, J. Chem. Phys. 84, 3963 (1986); Thijssen,
//! "Computational Physics", ch. 15).
//!
//! All routines are primitive-first: the OS recurrences run per primitive
//! combination and the shell-level wrappers contract with the coefficients
//! and component normalizations.
//!
//! Recursion families:
//!
//! * **Overlap** — per-axis 1D OS recurrence.
//! * **Kinetic** — derivatives of the overlap (`⟨a|−½∇²|b⟩`).
//! * **Nuclear attraction** — OS with Boys base
//!   `(00)^m = (2π/p)·K_ab·F_m(T)` and the `q→∞` limit factors.
//! * **Electron repulsion** — the full two-electron OS recursion, built
//!   from the same raising machinery on both electrons.
//!
//! Atomic units throughout.

use tpt_chem_core::num;
use tpt_chem_core::vec3::Vec3;

use crate::boys::boys_array;
use crate::gaussian::{all_components_up_to, cartesian_components, primitive_normalization, Shell};

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
    /// `exp(−αβ/p·|A−B|²)`.
    kab: f64,
}

fn make_pair(alpha: f64, a: Vec3, beta: f64, b: Vec3) -> Pair {
    let p = alpha + beta;
    let ab = a - b;
    let kab = num::exp(-alpha * beta / p * ab.norm_sq());
    let center = (a * alpha + b * beta) / p;
    Pair {
        p,
        alpha,
        beta,
        center,
        pa: center - a,
        ab,
        kab,
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

/// The OS "raising family": all m-arrays `V(a, b)^m` for
/// `|a| ≤ la`, `|b| ≤ lb`, given the `(0,0)^m` base array.
///
/// * `pa`, `ab`, `wp` — the `P−A`, `A−B`, `W−P` vectors.
/// * `reduce` — the factor `q/(p+q)` multiplying the `^{m+1}` correction
///   of the lowering terms (ERI electron 1) or `1` (nuclear attraction,
///   the `q→∞` limit).
///
/// Recurrences (Obara–Saika):
/// * raise: `V(a+1,b)^m = PA·V(a,b)^m + WP·V(a+1,b)^{m+1}
///   + a/(2p)·[V(a−1,b)^m − R·V(a−1,b)^{m+1}]
///   + b/(2p)·[V(a,b−1)^m − R·V(a,b−1)^{m+1}]`
/// * horizontal: `(a, b+1) = (a+1, b) − AB·(a, b)`
fn raise_family(
    la: u8,
    lb: u8,
    p: f64,
    pa: Vec3,
    ab: Vec3,
    wp: Vec3,
    reduce: f64,
    base: &[f64],
) -> Vec<Vec<Vec<f64>>> {
    let max_m = base.len() - 1;
    // Intermediates are generated up to `la + 1`: the horizontal
    // transfer (0, b) = (1, b-1) - AB*(0, b-1) needs one raising step
    // beyond `la` even when `la == 0`. Consumers only read |a| <= la.
    let a_gen = all_components_up_to(la + 1);
    let b_comps = all_components_up_to(lb);
    let mut v = vec![vec![vec![0.0f64; max_m + 1]; b_comps.len()]; a_gen.len()];
    v[0][0].copy_from_slice(base);
    let inv2p = 0.5 / p;

    let total = |c: (u8, u8, u8)| u16::from(c.0) + u16::from(c.1) + u16::from(c.2);
    let idx = |comps: &[(u8, u8, u8)], c: (u8, u8, u8)| -> usize {
        comps
            .iter()
            .position(|&x| x == c)
            .expect("component generated earlier")
    };

    for t in 1..=(total((la, 0, 0)) + total((0, lb, 0))) as u16 {
        // Phase 1: raising cases (a ≠ 0).
        for (ia, &a) in a_gen.iter().enumerate() {
            for (ib, &b) in b_comps.iter().enumerate() {
                if total(a) + total(b) != t || total(a) == 0 {
                    continue;
                }
                let axis = if a.0 > 0 {
                    0
                } else if a.1 > 0 {
                    1
                } else {
                    2
                };
                let mut a_down = a;
                a_down.0 -= u8::from(axis == 0);
                a_down.1 -= u8::from(axis == 1);
                a_down.2 -= u8::from(axis == 2);
                let src = &v[idx(&a_gen, a_down)][ib];
                // Optional lowering of b on the same axis.
                let (b_down_idx, b_n): (Option<usize>, f64) = match axis {
                    0 if b.0 > 0 => (Some(idx(&b_comps, (b.0 - 1, b.1, b.2))), f64::from(b.0)),
                    1 if b.1 > 0 => (Some(idx(&b_comps, (b.0, b.1 - 1, b.2))), f64::from(b.1)),
                    2 if b.2 > 0 => (Some(idx(&b_comps, (b.0, b.1, b.2 - 1))), f64::from(b.2)),
                    _ => (None, 0.0),
                };
                let pa_axis = pa.get(axis);
                let wp_axis = wp.get(axis);
                // OS lowering coefficient: the recurrence's a_i is the
                // raised (target) index in the OS convention.
                let n_a = match axis {
                    0 => f64::from(a.0),
                    1 => f64::from(a.1),
                    _ => f64::from(a.2),
                };
                let mut arr = vec![0.0; max_m + 1];
                for m in (0..=max_m).rev() {
                    let hi = |arr_hi: &[f64]| arr_hi.get(m + 1).copied().unwrap_or(0.0);
                    let mut val = pa_axis * src[m] + wp_axis * hi(&arr);
                    if n_a > 0.0 {
                        val += n_a * inv2p * (src[m] - reduce * hi(src));
                    }
                    if let Some(bdi) = b_down_idx {
                        let bd = &v[ia][bdi];
                        val += b_n * inv2p * (bd[m] - reduce * hi(bd));
                    }
                    arr[m] = val;
                }
                v[ia][ib] = arr;
            }
        }
        // Phase 2: horizontal cases (a == 0, b ≠ 0):
        // (0, b) = (1, b−1) − AB·(0, b−1).
        for (ia, &a) in a_gen.iter().enumerate() {
            for (ib, &b) in b_comps.iter().enumerate() {
                if total(a) + total(b) != t || total(a) != 0 || total(b) == 0 {
                    continue;
                }
                let baxis = if b.0 > 0 {
                    0
                } else if b.1 > 0 {
                    1
                } else {
                    2
                };
                let mut a_up = a;
                a_up.0 += u8::from(baxis == 0);
                a_up.1 += u8::from(baxis == 1);
                a_up.2 += u8::from(baxis == 2);
                let mut b_down = b;
                b_down.0 -= u8::from(baxis == 0);
                b_down.1 -= u8::from(baxis == 1);
                b_down.2 -= u8::from(baxis == 2);
                let up = &v[idx(&a_gen, a_up)][idx(&b_comps, b_down)];
                let prev = &v[ia][idx(&b_comps, b_down)];
                let ab_axis = ab.get(baxis);
                let mut arr = vec![0.0; max_m + 1];
                for m in 0..=max_m {
                    arr[m] = up[m] + ab_axis * prev[m];
                }
                v[ia][ib] = arr;
            }
        }
    }
    v
}

/// Primitive ERI `(ab|cd)` for all component combinations, via the
/// Obara-Saika two-electron recursion: electron 1 is raised first
/// (`V1[a][b]` with `c = d = 0`), then electron 2 per `(a, b)` pair.
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
    let p_center = pair1.center;
    let q_center = pair2.center;
    let w_total = (p_center * p + q_center * q) / (p + q);
    let max_m = (la + lb + lc + ld) as usize;
    let t = p * q / (p + q) * (p_center - q_center).norm_sq();
    // The WP raise at the top of the (0,0,1) m-array reaches F_{max_m+1},
    // so the base carries one extra Boys term.
    let boys = boys_array(t, max_m + 1);
    let base: Vec<f64> = boys
        .iter()
        .map(|&f| {
            2.0 * num::powf(core::f64::consts::PI, 2.5) / (p * q * num::sqrt(p + q))
                * pair1.kab
                * pair2.kab
                * f
        })
        .collect();
    let reduce1 = q / (p + q);
    let v1 = raise_family(
        la,
        lb,
        p,
        pair1.pa,
        pair1.ab,
        w_total - p_center,
        reduce1,
        &base,
    );
    let a_comps = all_components_up_to(la);
    let b_comps = all_components_up_to(lb);
    let c_comps = all_components_up_to(lc);
    let d_comps = all_components_up_to(ld);
    let a_gen = all_components_up_to(la + 1);
    let mut out =
        vec![vec![vec![vec![0.0f64; d_comps.len()]; c_comps.len()]; b_comps.len()]; a_comps.len()];
    for (ia, &a) in a_comps.iter().enumerate() {
        for (ib, &b) in b_comps.iter().enumerate() {
            let used_m = (total(a) + total(b)) as usize;
            if used_m > max_m {
                continue;
            }
            let ia_gen = a_gen
                .iter()
                .position(|&x| x == a)
                .expect("component within la");
            let sub_base = &v1[ia_gen][ib][used_m..];
            let v2 = raise_family(
                lc,
                ld,
                q,
                pair2.pa,
                pair2.ab,
                w_total - q_center,
                p / (p + q),
                sub_base,
            );
            for (ic, _) in c_comps.iter().enumerate() {
                for (id, _) in d_comps.iter().enumerate() {
                    out[ia][ib][ic][id] = v2[ic][id][0];
                }
            }
        }
    }
    out
}

fn total(c: (u8, u8, u8)) -> u8 {
    c.0 + c.1 + c.2
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
fn nuclear_value(
    pair: &Pair,
    c: Vec3,
    la: u8,
    lb: u8,
    comp_a: (u8, u8, u8),
    comp_b: (u8, u8, u8),
) -> f64 {
    // The i-raise path is the validated one; canonicalize so the
    // higher-momentum side is always the a-side.
    let ta = u16::from(comp_a.0) + u16::from(comp_a.1) + u16::from(comp_a.2);
    let tb = u16::from(comp_b.0) + u16::from(comp_b.1) + u16::from(comp_b.2);
    if tb > ta {
        let a_center = pair.center - pair.pa;
        let b_center = a_center - pair.ab;
        let flipped = make_pair(pair.beta, b_center, pair.alpha, a_center);
        return nuclear_value(&flipped, c, lb, la, comp_b, comp_a);
    }
    let p = pair.p;
    let pa = pair.pa; // P - A
    let pb = pa - pair.ab; // P - B
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
mod dbg_tmp {
    use super::*;
    use crate::gaussian::Shell;

    #[test]
    fn dbg_kinetic_axes() {
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(0.0, 0.0, 0.8);
        let mut sa = Shell::new(0, a, vec![(1.0, 1.0)]);
        sa.normalize();
        let mut sb = Shell::new(0, b, vec![(1.0, 1.0)]);
        sb.normalize();
        let pi = Prim {
            alpha: 1.0,
            coeff: 1.0,
        };
        let pj = Prim {
            alpha: 1.0,
            coeff: 1.0,
        };
        let pair = make_pair(1.0, a, 1.0, b);
        let global = num::powf(core::f64::consts::PI / pair.p, 1.5);
        let tables = overlap_tables(1.0, 1.0, &pair, 0, 2);
        for axis in 0..3 {
            for dj in [0i32, 2] {
                let comp = match (axis, dj) {
                    (0, d) => (d as u8, 0, 0),
                    (1, d) => (0, d as u8, 0),
                    (2, d) => (0, 0, d as u8),
                    _ => unreachable!(),
                };
                println!(
                    "axis {axis} dj {dj}: table = {:.6}, with global = {:.6}",
                    table_val(&tables, (0, 0, 0), comp),
                    global * table_val(&tables, (0, 0, 0), comp)
                );
            }
        }
        println!("global = {global:.6}");
        let _ = (sa, sb, pi, pj);
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
}

#[cfg(test)]
mod dbg_p {
    use super::*;
    use crate::gaussian::Shell;

    #[test]
    fn dbg_p_kinetic_axes() {
        let alpha = 1.3f64;
        let mut sh = Shell::new(1, Vec3::ZERO, vec![(alpha, 1.0)]);
        sh.normalize();
        println!("shell_norm = {:?}", sh.normalizations);
        println!(
            "prim_norm(0,0,1) = {}",
            primitive_normalization(alpha, (0, 0, 1))
        );
        let pi = Prim { alpha, coeff: 1.0 };
        let pair = make_pair(alpha, Vec3::ZERO, alpha, Vec3::ZERO);
        let global = num::powf(core::f64::consts::PI / pair.p, 1.5);
        println!("p = {}, global = {:.6}", pair.p, global);
        let tables = overlap_tables(alpha, alpha, &pair, 1, 3);
        // z-axis table rows: s[0], s[1] for j = 0..=3.
        for (i, row) in tables[2].iter().take(2).enumerate() {
            println!(
                "z-table s[{i}] = {:?}",
                row.iter()
                    .map(|x| (x * 1e4).round() / 1e4)
                    .collect::<Vec<_>>()
            );
        }
        for (i, row) in tables[0].iter().take(2).enumerate() {
            println!(
                "x-table s[{i}] = {:?}",
                row.iter()
                    .map(|x| (x * 1e4).round() / 1e4)
                    .collect::<Vec<_>>()
            );
        }
        let t = kinetic_shell(&sh, &sh);
        println!("T(pz,pz) = {:.6}, T(px,px) = {:.6}", t[0][0], t[2][2]);
        let _ = (pi, sh);
    }
}

#[cfg(test)]
mod dbg_v {
    use super::*;
    use crate::gaussian::Shell;

    #[test]
    fn dbg_pz_self_v() {
        let alpha = 1.3f64;
        let c = Vec3::new(0.2, -0.4, 0.6);
        let mut sh = Shell::new(1, c, vec![(alpha, 1.0)]);
        sh.normalize();
        println!("shell_norm = {}", sh.normalizations[0]);
        let _pi = Prim { alpha, coeff: 1.0 };
        let _pj = Prim { alpha, coeff: 1.0 };
        let pair = make_pair(alpha, c, alpha, c);
        println!(
            "p = {}, kab = {}, pa = {:?}",
            pair.p,
            pair.kab,
            (pair.pa.x, pair.pa.y, pair.pa.z)
        );
        let max_m = 1usize;
        let t = pair.p * (pair.center - c).norm_sq();
        println!("T = {t}");
        let boys = boys_array(t, max_m);
        let base: Vec<f64> = boys
            .iter()
            .map(|&f| 2.0 * core::f64::consts::PI / pair.p * pair.kab * f)
            .collect();
        println!("base = {base:?}");
        let tables = raise_family(1, 0, pair.p, pair.pa, pair.ab, Vec3::ZERO, 1.0, &base);
        let a_gen = crate::gaussian::all_components_up_to(2);
        for (i, comp) in a_gen.iter().enumerate() {
            println!("v[{i}] ({comp:?}) m0 = {:.6}", tables[i][0][0]);
        }
        let nv = nuclear_value(
            &make_pair(alpha, c, alpha, c),
            c,
            1,
            0,
            (0, 0, 1),
            (0, 0, 0),
        );
        println!("nuclear_value raw = {nv:.6}");
        let na = crate::gaussian::primitive_normalization(alpha, (0, 0, 1));
        let nb = crate::gaussian::primitive_normalization(alpha, (0, 0, 1));
        println!("na = {na:.6}, nb = {nb:.6}, product = {:.6}", na * nb);
        println!("expected V = {:.6}", -nv * na * nb);
        let v = nuclear_shell(&sh, &sh, &[(1.0, c)]);
        println!("V(pz,pz) = {:.6}", v[0][0]);
        println!("shell normalizations = {:?}", sh.normalizations);
    }
}
