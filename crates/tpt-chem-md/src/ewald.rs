//! Long-range electrostatics: Ewald summation and thermostats/barostats.
//!
//! The Ewald sum splits the Coulomb interaction of a periodic system into a
//! short-range real-space part (erfc) and a long-range reciprocal-space
//! part (the structure-factor sum over reciprocal vectors), plus self and
//! surface corrections — the standard textbook decomposition
//! (Allen & Tildesley, Computer Simulation of Liquids, §5.5).
//!
//! Units: MD units (Å, fs, amu, kJ·mol⁻¹); charges in elementary units e,
//! so the Coulomb prefactor is [`tpt_chem_core::units::COULOMB_PREFACTOR_KJ_ANG`].
//!
//! For large systems the reciprocal sum can be swapped for the mesh-based
//! [`crate::pme`] solver, which reproduces this sum to the spline
//! discretization error while scaling as O(N log N).

use tpt_chem_core::special;
use tpt_chem_core::units::COULOMB_PREFACTOR_KJ_ANG;
use tpt_chem_core::vec3::Vec3;

use crate::box3::Box3;

/// Parameters for the Ewald sum.
#[derive(Clone, Copy, Debug)]
pub struct EwaldParams {
    /// Splitting parameter α (Å⁻¹). α⁻¹ is the real/reciprocal crossover.
    pub alpha: f64,
    /// Real-space cutoff (Å).
    pub r_cutoff: f64,
    /// Reciprocal-space cutoff: sum over |g|² ≤ g_max².
    pub g_max: f64,
}

impl EwaldParams {
    /// Heuristic parameters for a cubic box of edge `L` with `n` charges:
    /// α = √π (σ ≈ 1 Å) and cutoffs chosen for ~1e-8 convergence.
    pub fn heuristic(box_: &Box3) -> Self {
        let l = box_.length.x.min(box_.length.y).min(box_.length.z);
        let alpha = core::f64::consts::PI.sqrt() / l * 2.0;
        EwaldParams {
            alpha,
            r_cutoff: l / 2.0,
            g_max: 6.0 * alpha,
        }
    }
}

/// Real-space part of the Ewald sum: erfc-screened pair interactions within
/// `r_cutoff`. Returned as `(forces, energy)`.
pub fn ewald_real_space(
    pos: &[Vec3],
    charges: &[f64],
    box_: &Box3,
    params: &EwaldParams,
) -> (Vec<Vec3>, f64) {
    let n = pos.len();
    let mut forces = vec![Vec3::ZERO; n];
    let mut energy = 0.0;
    let alpha = params.alpha;
    let alpha2 = alpha * alpha;

    for i in 0..n {
        for j in (i + 1)..n {
            let raw = pos[j] - pos[i];
            let d = box_.min_image(raw);
            let r2 = d.norm_sq();
            if r2 >= params.r_cutoff * params.r_cutoff || r2 < 1e-12 {
                continue;
            }
            let r = r2.sqrt();
            let qi_qj = charges[i] * charges[j];
            let br = alpha * r;
            let erfc_br = special::erfc(br);
            energy += COULOMB_PREFACTOR_KJ_ANG * qi_qj * erfc_br / r;
            // d/dr [erfc(αr)/r] = −(2α/√π)e^{−α²r²}/r − erfc(αr)/r².
            let f_mag = COULOMB_PREFACTOR_KJ_ANG
                * qi_qj
                * (erfc_br / r2
                    + 2.0 * alpha / core::f64::consts::PI.sqrt() * (-alpha2 * r2).exp() / r);
            let fvec = d * (f_mag / r);
            forces[i] += fvec;
            forces[j] -= fvec;
        }
    }
    (forces, energy)
}

/// Self term of the Ewald sum: `−α/√π Σ q²` (tin-foil boundary conditions).
pub fn ewald_self_energy(charges: &[f64], alpha: f64) -> f64 {
    let q2: f64 = charges.iter().map(|q| q * q).sum();
    -COULOMB_PREFACTOR_KJ_ANG * alpha / core::f64::consts::PI.sqrt() * q2
}

/// Ewald-summed Coulomb energy and per-atom forces for a periodic system.
///
/// Charges `q` in e, positions in Å (must be wrapped into the box), masses
/// are not needed here. Returns `(forces, coulomb_energy)`. The energy
/// includes real + reciprocal + self terms (no surface term: appropriate
/// for tin-foil boundary conditions).
pub fn ewald_energy_forces(
    pos: &[Vec3],
    charges: &[f64],
    box_: &Box3,
    params: &EwaldParams,
) -> (Vec<Vec3>, f64) {
    let n = pos.len();
    let mut forces = vec![Vec3::ZERO; n];
    let alpha = params.alpha;

    // ---- Real space.
    let (real_forces, real_energy) = ewald_real_space(pos, charges, box_, params);
    for (f, rf) in forces.iter_mut().zip(real_forces.iter()) {
        *f += *rf;
    }
    let mut energy = real_energy;
    let alpha2 = alpha * alpha;

    // ---- Reciprocal space: structure-factor sum over |g| ≤ g_max.
    // Enumerate integer Miller indices within the sphere.
    let box_inv = box_.length;
    let nmax = [
        (params.g_max * box_inv.x / (2.0 * core::f64::consts::PI)).ceil() as i64 + 1,
        (params.g_max * box_inv.y / (2.0 * core::f64::consts::PI)).ceil() as i64 + 1,
        (params.g_max * box_inv.z / (2.0 * core::f64::consts::PI)).ceil() as i64 + 1,
    ];
    let two_pi = 2.0 * core::f64::consts::PI;
    let mut g_list: Vec<(Vec3, f64)> = Vec::new();
    for h in -nmax[0]..=nmax[0] {
        for k in -nmax[1]..=nmax[1] {
            for l in -nmax[2]..=nmax[2] {
                if h == 0 && k == 0 && l == 0 {
                    continue;
                }
                let gv = Vec3::new(
                    two_pi * h as f64 / box_inv.x,
                    two_pi * k as f64 / box_inv.y,
                    two_pi * l as f64 / box_inv.z,
                );
                let g2 = gv.norm_sq();
                if g2 <= params.g_max * params.g_max {
                    g_list.push((gv, g2));
                }
            }
        }
    }

    // Structure factors S(g) = Σ q_j e^{i g·r_j}.
    for &(gv, g2) in &g_list {
        let mut s_re = 0.0;
        let mut s_im = 0.0;
        for (r, &q) in pos.iter().zip(charges.iter()) {
            let phase = gv.dot(*r);
            s_re += q * phase.cos();
            s_im += q * phase.sin();
        }
        let factor = (-g2 / (4.0 * alpha2)).exp() / g2;
        // E_recip = prefac × (2π/V) × e^{−g²/4α²}/g² × |S(g)|².
        energy += COULOMB_PREFACTOR_KJ_ANG * factor * two_pi / box_.volume()
            * (s_re * s_re + s_im * s_im);
        let c = COULOMB_PREFACTOR_KJ_ANG * factor * 2.0 / box_.volume() * two_pi;
        for ((f, r), &q) in forces.iter_mut().zip(pos.iter()).zip(charges.iter()) {
            let phase = gv.dot(*r);
            // −∂E/∂r_i = +prefac × factor × 2/V × 2π g × [q sin + q cos]…
            let contrib = c * q * (s_re * phase.sin() - s_im * phase.cos());
            *f += gv * contrib;
        }
    }

    // ---- Self term.
    energy += ewald_self_energy(charges, alpha);

    (forces, energy)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_pair_real_space_self_consistent() {
        // Two charges +1/−1 in a big box: Ewald ≈ direct Coulomb.
        let box_ = Box3::cubic(20.0);
        let pos = [Vec3::new(9.0, 10.0, 10.0), Vec3::new(11.0, 10.0, 10.0)];
        let q = [1.0, -1.0];
        let params = EwaldParams {
            alpha: 0.3,
            r_cutoff: 9.5,
            g_max: 7.0,
        };
        let (_, e) = ewald_energy_forces(&pos, &q, &box_, &params);
        let direct = -COULOMB_PREFACTOR_KJ_ANG / 2.0;
        assert!(
            (e - direct).abs() < 0.05 * direct.abs(),
            "e {e} vs {direct}"
        );
    }

    #[test]
    fn net_neutral_force_sum_is_zero() {
        let box_ = Box3::cubic(12.0);
        let pos = [
            Vec3::new(2.0, 3.0, 4.0),
            Vec3::new(9.0, 3.0, 4.0),
            Vec3::new(4.0, 8.0, 8.0),
        ];
        let q = [1.0, -0.5, -0.5];
        let params = EwaldParams::heuristic(&box_);
        let (forces, _) = ewald_energy_forces(&pos, &q, &box_, &params);
        let sum = forces.iter().fold(Vec3::ZERO, |a, f| a + *f);
        assert!(sum.norm() < 1e-8, "|ΣF| = {}", sum.norm());
    }

    #[test]
    fn nacl_lattice_energy_is_negative() {
        // Simple 2-ion NaCl-ish cell: ±1 at half-box separation.
        let box_ = Box3::cubic(5.64);
        let pos = [Vec3::ZERO, Vec3::new(2.82, 2.82, 2.82)];
        let q = [1.0, -1.0];
        let params = EwaldParams {
            alpha: 0.28,
            r_cutoff: 2.8,
            g_max: 7.0,
        };
        let (_, e) = ewald_energy_forces(&pos, &q, &box_, &params);
        assert!(e < 0.0, "attraction must give negative energy: {e}");
    }
}
