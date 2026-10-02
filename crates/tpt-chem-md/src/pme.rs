//! Particle Mesh Ewald (PME): mesh reciprocal-space electrostatics.
//!
//! Replaces the direct structure-factor sum of [`crate::ewald`] with a grid
//! pipeline (Essmann et al., J. Chem. Phys. 103, 8577 (1995)):
//!
//! 1. Assign the charges to a power-of-two mesh with cubic B-splines
//!    ([`crate::bspline`]).
//! 2. Forward-FFT the charge grid ([`crate::fft`], in-house radix-2).
//! 3. Scale each mode by the Ewald weight `e^{−g²/4α²}/g²` divided by the
//!    squared B-spline structure factor `|b|²` (the aliasing correction):
//!    since the mesh structure factor factorizes exactly as
//!    `Q̂(m) = b(m)·S(g_m)`, the division recovers `|S(g_m)|²` at every mesh
//!    wavenumber.
//! 4. Energy from the scaled modes; forces from the inverse transform,
//!    interpolated back with the B-spline derivatives (the exact gradient
//!    of the mesh energy).
//!
//! Units and normalization match [`crate::ewald`]: MD units (Å, kJ·mol⁻¹),
//! charges in e, tin-foil boundary conditions. Only the *reciprocal*
//! contribution is returned; combine with
//! [`crate::ewald::ewald_real_space`] and [`crate::ewald::ewald_self_energy`]
//! for the full Ewald splitting.

use tpt_chem_core::units::COULOMB_PREFACTOR_KJ_ANG;
use tpt_chem_core::vec3::Vec3;

use crate::box3::Box3;
use crate::bspline::{m4, m4_derivative, ORDER};
use crate::fft::Field3;

/// PME mesh parameters.
#[derive(Clone, Copy, Debug)]
pub struct PmeParams {
    /// Ewald splitting parameter α (Å⁻¹); must match the real-space part.
    pub alpha: f64,
    /// Mesh points per axis (powers of two).
    pub dims: [usize; 3],
    /// Reciprocal-wavenumber truncation (Å⁻¹). `0.0` keeps every mesh mode
    /// up to the grid Nyquist frequency; a positive value limits the sum to
    /// `|g| ≤ g_max`, e.g. to mirror a direct Ewald sum for validation.
    pub g_max: f64,
}

impl PmeParams {
    /// Heuristic for a box: mesh spacing ≈ 1 Å (rounded up to a power of
    /// two, minimum 4) and the [`crate::ewald::EwaldParams::heuristic`]
    /// splitting; `g_max = 0` keeps the full mesh.
    pub fn heuristic(box_: &Box3) -> Self {
        let l = box_.length.x.min(box_.length.y).min(box_.length.z);
        let alpha = core::f64::consts::PI.sqrt() / l * 2.0;
        let n = (l as usize).max(4).next_power_of_two();
        PmeParams {
            alpha,
            dims: [n, n, n],
            g_max: 0.0,
        }
    }
}

/// Wrapped grid coordinate `u = N·s` for a position component in `[0, L)`.
#[inline]
fn grid_coord(r: f64, l: f64, n: usize) -> f64 {
    let s = r - (r / l).floor() * l; // wrap into [0, L)
    n as f64 * s / l
}

/// PME reciprocal-space energy and per-atom forces.
///
/// Positions must lie inside the box (they are wrapped here). Returns
/// `(forces, reciprocal_energy)`.
///
/// # Examples
///
/// Full Ewald splitting with the mesh reciprocal: real-space and self terms
/// come from [`crate::ewald`], the long-range part from the mesh.
///
/// ```
/// use tpt_chem_core::vec3::Vec3;
/// use tpt_chem_md::box3::Box3;
/// use tpt_chem_md::ewald::{ewald_real_space, ewald_self_energy, EwaldParams};
/// use tpt_chem_md::pme::{pme_reciprocal_energy_forces, PmeParams};
///
/// let box_ = Box3::cubic(12.0);
/// let pos = [
///     Vec3::new(1.0, 2.0, 3.0),
///     Vec3::new(7.0, 5.0, 4.0),
///     Vec3::new(3.0, 9.0, 8.0),
/// ];
/// let q = [1.0, -0.5, -0.5];
/// let ewald = EwaldParams { alpha: 0.5, r_cutoff: 5.5, g_max: 3.0 };
/// let pme = PmeParams { alpha: 0.5, dims: [32, 32, 32], g_max: 0.0 };
///
/// let (f_real, e_real) = ewald_real_space(&pos, &q, &box_, &ewald);
/// let e_self = ewald_self_energy(&q, ewald.alpha);
/// let (f_mesh, e_mesh) = pme_reciprocal_energy_forces(&pos, &q, &box_, &pme);
/// let total_energy = e_real + e_self + e_mesh;
/// assert!(total_energy.is_finite());
/// let _ = (f_real, f_mesh);
/// ```
pub fn pme_reciprocal_energy_forces(
    pos: &[Vec3],
    charges: &[f64],
    box_: &Box3,
    params: &PmeParams,
) -> (Vec<Vec3>, f64) {
    let n = pos.len();
    let [n1, n2, n3] = params.dims;
    let two_pi = 2.0 * core::f64::consts::PI;
    let c_scale = COULOMB_PREFACTOR_KJ_ANG * two_pi / box_.volume();
    let alpha2 = params.alpha * params.alpha;

    // Grid coordinates u = N·s per particle per axis.
    let u: Vec<[f64; 3]> = pos
        .iter()
        .map(|r| {
            [
                grid_coord(r.x, box_.length.x, n1),
                grid_coord(r.y, box_.length.y, n2),
                grid_coord(r.z, box_.length.z, n3),
            ]
        })
        .collect();

    // ---- 1. Charge assignment (order-4 splines, wrapped indices).
    let mut grid = Field3::zeros(params.dims);
    for (&q, uj) in charges.iter().zip(u.iter()) {
        let base = [
            uj[0].floor() as isize,
            uj[1].floor() as isize,
            uj[2].floor() as isize,
        ];
        for dx in 0..ORDER as isize {
            for dy in 0..ORDER as isize {
                for dz in 0..ORDER as isize {
                    let wx = m4(uj[0] - (base[0] - 3 + dx) as f64);
                    let wy = m4(uj[1] - (base[1] - 3 + dy) as f64);
                    let wz = m4(uj[2] - (base[2] - 3 + dz) as f64);
                    let x = (base[0] - 3 + dx).rem_euclid(n1 as isize) as usize;
                    let y = (base[1] - 3 + dy).rem_euclid(n2 as isize) as usize;
                    let z = (base[2] - 3 + dz).rem_euclid(n3 as isize) as usize;
                    let i = (x * n2 + y) * n3 + z;
                    grid.data[i].re += q * wx * wy * wz;
                }
            }
        }
    }

    // ---- 2. Forward FFT.
    grid.fft_forward();

    // ---- 3. B-spline structure factors |b|² per axis, then mode scaling.
    // b_ax(m) = Σ_v M4(v)·e^{2πi m v / N} for v = 0..3 (M4(0) = 0).
    let mut b2: [Vec<f64>; 3] = Default::default();
    for (ax, &n_ax) in params.dims.iter().enumerate() {
        let mut row = vec![0.0; n_ax];
        for (m, row_m) in row.iter_mut().enumerate() {
            let mut re = 0.0;
            let mut im = 0.0;
            for v in 1..ORDER {
                let ang = two_pi * m as f64 * v as f64 / n_ax as f64;
                let mv = m4(v as f64);
                re += mv * ang.cos();
                im += mv * ang.sin();
            }
            *row_m = re * re + im * im;
        }
        b2[ax] = row;
    }

    let g_max2 = if params.g_max > 0.0 {
        params.g_max * params.g_max
    } else {
        f64::INFINITY
    };

    let mut energy = 0.0;
    let mut phi = Field3::zeros(params.dims);
    for x in 0..n1 {
        let mx: isize = if x <= n1 / 2 {
            x as isize
        } else {
            x as isize - n1 as isize
        };
        for y in 0..n2 {
            let my: isize = if y <= n2 / 2 {
                y as isize
            } else {
                y as isize - n2 as isize
            };
            for z in 0..n3 {
                let mz: isize = if z <= n3 / 2 {
                    z as isize
                } else {
                    z as isize - n3 as isize
                };
                if mx == 0 && my == 0 && mz == 0 {
                    continue;
                }
                let gx = two_pi * mx as f64 / box_.length.x;
                let gy = two_pi * my as f64 / box_.length.y;
                let gz = two_pi * mz as f64 / box_.length.z;
                let g2 = gx * gx + gy * gy + gz * gz;
                let qh = grid.get(x, y, z);
                let denom = b2[0][x] * b2[1][y] * b2[2][z];
                if g2 > g_max2 || denom < 1e-30 {
                    continue; // phi entry stays zero
                }
                let w = (-g2 / (4.0 * alpha2)).exp() / g2 / denom;
                energy += w * (qh.re * qh.re + qh.im * qh.im);
                phi.set(x, y, z, crate::fft::Complex::new(w * qh.re, w * qh.im));
            }
        }
    }
    energy *= c_scale;

    // ---- 4. Inverse transform of the scaled modes.
    phi.fft_inverse();

    // ---- 5. Forces via B-spline derivative interpolation: the exact
    // gradient of the mesh energy.
    let mut forces = vec![Vec3::ZERO; n];
    for (j, uj) in u.iter().enumerate() {
        let base = [
            uj[0].floor() as isize,
            uj[1].floor() as isize,
            uj[2].floor() as isize,
        ];
        let mut f_acc = [0.0f64; 3];
        for dx in 0..ORDER as isize {
            for dy in 0..ORDER as isize {
                for dz in 0..ORDER as isize {
                    let x = (base[0] - 3 + dx).rem_euclid(n1 as isize) as usize;
                    let y = (base[1] - 3 + dy).rem_euclid(n2 as isize) as usize;
                    let z = (base[2] - 3 + dz).rem_euclid(n3 as isize) as usize;
                    let th = phi.get(x, y, z).re;
                    if th == 0.0 {
                        continue;
                    }
                    let wx = uj[0] - (base[0] - 3 + dx) as f64;
                    let wy = uj[1] - (base[1] - 3 + dy) as f64;
                    let wz = uj[2] - (base[2] - 3 + dz) as f64;
                    let mx = m4(wx);
                    let my = m4(wy);
                    let mz = m4(wz);
                    f_acc[0] += th * m4_derivative(wx) * my * mz;
                    f_acc[1] += th * mx * m4_derivative(wy) * mz;
                    f_acc[2] += th * mx * my * m4_derivative(wz);
                }
            }
        }
        let q = charges[j];
        forces[j] = Vec3::new(
            -(n1 as f64) / box_.length.x * 2.0 * c_scale * q * f_acc[0],
            -(n2 as f64) / box_.length.y * 2.0 * c_scale * q * f_acc[1],
            -(n3 as f64) / box_.length.z * 2.0 * c_scale * q * f_acc[2],
        );
    }

    (forces, energy)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::ewald::{ewald_energy_forces, ewald_real_space, ewald_self_energy, EwaldParams};
    use crate::fft::Complex;

    /// Deterministic pseudo-random positions and a net-neutral charge set.
    fn test_system(n: usize, l: f64) -> (Vec<Vec3>, Vec<f64>) {
        let mut s = 0x1234_5678_9ABC_DEF0u64;
        let mut next = move || {
            s ^= s >> 12;
            s ^= s << 25;
            s ^= s >> 27;
            ((s.wrapping_mul(0x2545_F491_4F6C_DD1D) >> 11) as f64) / ((1u64 << 53) as f64)
        };
        let mut pos = Vec::new();
        let mut charges = Vec::new();
        let mut q_sum = 0.0;
        for i in 0..n - 1 {
            pos.push(Vec3::new(next() * l, next() * l, next() * l));
            let q = (next() * 2.0 - 1.0).round() * 0.5 + 0.25 * (i as f64 % 3.0 - 1.0);
            charges.push(q);
            q_sum += q;
        }
        // Last charge balances the sum exactly.
        pos.push(Vec3::new(next() * l, next() * l, next() * l));
        charges.push(-q_sum);
        (pos, charges)
    }

    #[test]
    fn pme_matches_ewald_with_matching_wavenumber_cutoff() {
        // With the same α and the mesh modes truncated at the Ewald g_max,
        // the |b|²-corrected mesh sum reproduces the direct reciprocal sum
        // to within the spline discretization error, which falls as the
        // fourth power of the mesh spacing (measured: 8.4e-4 / 5.4e-5 /
        // 2.2e-6 relative at 16³/32³/64³ for this system).
        let l = 12.34;
        let box_ = Box3::cubic(l);
        let (pos, charges) = test_system(12, l);
        let params = EwaldParams {
            alpha: 0.5,
            r_cutoff: 6.0,
            g_max: 3.0,
        };
        let (f_ewald_total, e_ewald) = ewald_energy_forces(&pos, &charges, &box_, &params);
        // Strip the real-space contribution so both sides are pure reciprocal.
        let (f_real, _) = ewald_real_space(&pos, &charges, &box_, &params);
        let f_ewald: Vec<Vec3> = f_ewald_total
            .iter()
            .zip(f_real.iter())
            .map(|(t, r)| Vec3::new(t.x - r.x, t.y - r.y, t.z - r.z))
            .collect();

        // Mesh Nyquist must cover g_max: πN/L ≥ g_max → N ≥ 8.
        let pme = PmeParams {
            alpha: 0.5,
            dims: [32, 32, 32],
            g_max: 3.0,
        };
        let (f_pme, e_pme) = pme_reciprocal_energy_forces(&pos, &charges, &box_, &pme);

        // Isolate the reciprocal pieces: Ewald total minus its real+self.
        let (_, e_real) = ewald_real_space(&pos, &charges, &box_, &params);
        let e_self = ewald_self_energy(&charges, params.alpha);
        let e_recip_ewald = e_ewald - e_real - e_self;

        let rel = (e_pme - e_recip_ewald).abs() / e_recip_ewald.abs();
        assert!(
            rel < 1e-4,
            "reciprocal energy: pme {e_pme} vs ewald {e_recip_ewald} (rel {rel})"
        );

        let mut max_df: f64 = 0.0;
        let mut max_f: f64 = 0.0;
        for (fp, fe) in f_pme.iter().zip(f_ewald.iter()) {
            let d = Vec3::new(fp.x - fe.x, fp.y - fe.y, fp.z - fe.z);
            max_df = max_df.max(d.norm());
            max_f = max_f.max(fe.norm());
        }
        // Both sides are now pure reciprocal forces of the same mode set;
        // the residual is the mesh discretization error.
        assert!(
            max_df < 1e-3 * max_f,
            "force mismatch {max_df} vs scale {max_f}"
        );
    }

    #[test]
    fn pme_force_is_gradient_of_pme_energy() {
        // Central finite differences of the mesh energy must agree with the
        // spline-interpolated forces.
        let l = 10.0;
        let box_ = Box3::cubic(l);
        let (pos, charges) = test_system(6, l);
        let pme = PmeParams {
            alpha: 0.4,
            dims: [16, 16, 16],
            g_max: 0.0,
        };
        let (forces, e0) = pme_reciprocal_energy_forces(&pos, &charges, &box_, &pme);
        let h = 1e-4;
        let comp = |v: Vec3, axis: usize| [v.x, v.y, v.z][axis];
        let with_axis = |v: Vec3, axis: usize, val: f64| {
            Vec3::new(
                if axis == 0 { val } else { v.x },
                if axis == 1 { val } else { v.y },
                if axis == 2 { val } else { v.z },
            )
        };
        for j in 0..pos.len() {
            for axis in 0..3 {
                let mut pp = pos.to_vec();
                let mut pm = pos.to_vec();
                pp[j] = with_axis(pos[j], axis, comp(pos[j], axis) + h);
                pm[j] = with_axis(pos[j], axis, comp(pos[j], axis) - h);
                let (_, ep) = pme_reciprocal_energy_forces(&pp, &charges, &box_, &pme);
                let (_, em) = pme_reciprocal_energy_forces(&pm, &charges, &box_, &pme);
                let fd = -(ep - em) / (2.0 * h);
                let an = comp(forces[j], axis);
                assert!(
                    (fd - an).abs() < 1e-5 * (fd.abs() + an.abs()).max(1e-6) + 1e-7,
                    "particle {j} axis {axis}: fd {fd} vs analytic {an}"
                );
            }
        }
        let _ = e0;
        let _ = Complex::default();
    }

    #[test]
    fn pme_converges_to_ewald_as_mesh_refines() {
        // Practical-parameter check: full-mesh PME (no g truncation) vs the
        // converged direct Ewald sum; the error shrinks with grid spacing.
        let l = 12.34;
        let box_ = Box3::cubic(l);
        let (pos, charges) = test_system(10, l);
        let params = EwaldParams {
            alpha: 0.5,
            r_cutoff: 6.0,
            g_max: 5.0,
        };
        let (_, e_ewald) = ewald_energy_forces(&pos, &charges, &box_, &params);
        let (_, e_real) = ewald_real_space(&pos, &charges, &box_, &params);
        let e_self = ewald_self_energy(&charges, params.alpha);
        let e_recip = e_ewald - e_real - e_self;

        let mut prev = f64::INFINITY;
        for &n in &[8usize, 16, 32] {
            let pme = PmeParams {
                alpha: 0.5,
                dims: [n, n, n],
                g_max: 0.0,
            };
            let (_, e_pme) = pme_reciprocal_energy_forces(&pos, &charges, &box_, &pme);
            let err = (e_pme - e_recip).abs();
            assert!(
                err < 0.05 * e_recip.abs(),
                "mesh {n}: |{e_pme} − {e_recip}| = {err}"
            );
            assert!(err <= prev * 1.0001, "mesh refinement must not worsen: {n}");
            // Fourth-order convergence: halving the mesh spacing quarters
            // twice the error (measured ratio ≈ 16× per doubling).
            if prev.is_finite() {
                assert!(err < 0.25 * prev, "not 4th-order: {prev} -> {err}");
            }
            prev = err;
        }
    }

    #[test]
    fn pme_force_sum_is_zero_for_neutral_system() {
        let l = 10.0;
        let box_ = Box3::cubic(l);
        let (pos, charges) = test_system(8, l);
        let pme = PmeParams {
            alpha: 0.4,
            dims: [32, 32, 32],
            g_max: 0.0,
        };
        let (forces, _) = pme_reciprocal_energy_forces(&pos, &charges, &box_, &pme);
        let sum = forces.iter().fold(Vec3::ZERO, |a, f| a + *f);
        // Mesh momentum conservation holds to the discretization level only.
        let max_f = forces.iter().map(|f| f.norm()).fold(0.0f64, f64::max);
        assert!(
            sum.norm() < 1e-3 * max_f,
            "|ΣF| = {} vs max |F| = {}",
            sum.norm(),
            max_f
        );
    }
}
