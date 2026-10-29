//! Energy minimisation: steepest descent with an adaptive line step, and
//! limited-memory BFGS ([`minimise_lbfgs`]) for faster convergence on
//! stiff or ill-conditioned geometries.
//!
//! Used to relax input geometries before dynamics (jittered lattices,
//! parsed structures with clashes). Forces and energies come from a
//! [`ForceModel`], so the minimised state is a stationary point of exactly
//! the model the subsequent run uses.

use std::collections::VecDeque;

use crate::forces::ForceModel;
use crate::system::System;
use tpt_chem_core::vec3::Vec3;

/// Outcome of [`minimise`].
#[derive(Clone, Copy, Debug)]
pub struct MinimiseReport {
    /// Line-search steps taken.
    pub steps: usize,
    /// Energy before the first step (kJ·mol⁻¹).
    pub initial_energy: f64,
    /// Energy at the final state (kJ·mol⁻¹).
    pub final_energy: f64,
    /// Largest per-atom force magnitude at the final state
    /// (kJ·mol⁻¹·Å⁻¹).
    pub final_max_force: f64,
    /// Whether the force-norm criterion was met before `max_steps`.
    pub converged: bool,
}

/// Steepest descent: `r ← r + s·F̂` with an adaptive step `s` (backtracking
/// on energy increase, growth on decrease). Positions are rewrapped into
/// the periodic box after every move.
///
/// Stops when the largest per-atom force magnitude falls below `f_tol`
/// (kJ·mol⁻¹·Å⁻¹) or after `max_steps` steps.
pub fn minimise(
    sys: &mut System,
    model: &ForceModel,
    max_steps: usize,
    f_tol: f64,
    initial_step: f64,
) -> MinimiseReport {
    let mut step = initial_step;
    let (f0, e0) = model.evaluate(sys);
    sys.force = f0;
    sys.potential_energy = e0;
    let initial_energy = e0;
    let mut energy = e0;
    let mut converged = false;
    let mut steps = 0usize;
    for _ in 0..max_steps {
        let max_f = sys.force.iter().map(|f| f.norm()).fold(0.0f64, f64::max);
        if max_f < f_tol {
            converged = true;
            break;
        }
        // Trial move along the (normalized) force direction.
        let moved: Vec<tpt_chem_core::vec3::Vec3> = sys
            .pos
            .iter()
            .zip(sys.force.iter())
            .map(|(r, f)| *r + *f * (step / f.norm().max(1e-30)))
            .collect();
        let saved = sys.pos.clone();
        sys.pos = moved;
        if let Some(b) = sys.box_ {
            for p in &mut sys.pos {
                *p = b.wrap(*p);
            }
        }
        let (_, e_trial) = model.evaluate(sys);
        if e_trial < energy {
            // Accept; grow the step (capped) for the next line search.
            energy = e_trial;
            step = (step * 1.5).min(0.5);
            sys.force = model.evaluate(sys).0;
            steps += 1;
        } else {
            // Reject; backtrack and retry within the same iteration budget.
            sys.pos = saved;
            step *= 0.5;
            if step < 1e-10 {
                break;
            }
        }
    }
    let final_max_force = sys.force.iter().map(|f| f.norm()).fold(0.0f64, f64::max);
    MinimiseReport {
        steps,
        initial_energy,
        final_energy: energy,
        final_max_force,
        converged,
    }
}

/// Inner product of two per-atom vector fields.
fn field_dot(a: &[Vec3], b: &[Vec3]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x.dot(*y)).sum()
}

/// Limited-memory BFGS over the atomic coordinates.
///
/// Keeps the last `memory` curvature pairs `(s, y)` and builds the search
/// direction with the standard two-loop recursion. Each trial step is
/// accepted under the Armijo sufficient-decrease condition, halving the step
/// on failure. `max_disp` caps the largest per-atom displacement of one
/// iteration (Å), which keeps early steps from flinging atoms across the
/// box. Stops under the same force criterion as [`minimise`], returning the
/// same [`MinimiseReport`].
///
/// If the direction is not a descent direction, or the line search fails
/// with a non-empty history, the history is discarded and the search restarts
/// from steepest descent. If steepest descent itself cannot satisfy the
/// Armijo condition, the search stops and `converged` stays `false`.
pub fn minimise_lbfgs(
    sys: &mut System,
    model: &ForceModel,
    max_steps: usize,
    f_tol: f64,
    max_disp: f64,
    memory: usize,
) -> MinimiseReport {
    const ARMIJO_C1: f64 = 1e-4;
    const MAX_BACKTRACKS: usize = 30;

    let (f0, e0) = model.evaluate(sys);
    sys.force = f0;
    sys.potential_energy = e0;
    let initial_energy = e0;
    let mut energy = e0;
    let mut converged = false;
    let mut steps = 0usize;
    // Curvature history: (s, y) with s = Δr and y = Δg, where g = −F.
    let mut history: VecDeque<(Vec<Vec3>, Vec<Vec3>)> = VecDeque::with_capacity(memory);

    for _ in 0..max_steps {
        let max_f = sys.force.iter().map(|f| f.norm()).fold(0.0f64, f64::max);
        if max_f < f_tol {
            converged = true;
            break;
        }
        let grad: Vec<Vec3> = sys.force.iter().map(|f| -*f).collect();

        // Two-loop recursion: d = −H·g with H the implicit inverse Hessian.
        let mut q = grad.clone();
        let mut alphas = Vec::with_capacity(history.len());
        for (s, y) in history.iter().rev() {
            let rho = 1.0 / field_dot(y, s);
            let a = rho * field_dot(s, &q);
            for (qi, yi) in q.iter_mut().zip(y) {
                *qi -= *yi * a;
            }
            alphas.push(a);
        }
        // Initial Hessian scaling γ = sᵀy / yᵀy from the newest pair.
        let gamma = match history.back() {
            Some((s, y)) => field_dot(s, y) / field_dot(y, y),
            None => 1.0,
        };
        let mut r: Vec<Vec3> = q.iter().map(|qi| *qi * gamma).collect();
        for ((s, y), a) in history.iter().zip(alphas.iter().rev()) {
            let rho = 1.0 / field_dot(y, s);
            let b = rho * field_dot(y, &r);
            for (ri, si) in r.iter_mut().zip(s) {
                *ri += *si * (a - b);
            }
        }
        let mut dir: Vec<Vec3> = r.iter().map(|ri| -*ri).collect();

        let mut slope = field_dot(&grad, &dir);
        if slope >= 0.0 {
            // Not a descent direction: fall back to steepest descent.
            history.clear();
            dir = grad.iter().map(|g| -*g).collect();
            slope = -field_dot(&grad, &grad);
        }
        let max_len = dir.iter().map(|d| d.norm()).fold(0.0f64, f64::max);
        if max_len > max_disp {
            let scale = max_disp / max_len;
            for d in &mut dir {
                *d = *d * scale;
            }
            slope *= scale;
        }

        // Armijo backtracking along `dir`.
        let saved = sys.pos.clone();
        let mut alpha = 1.0;
        let mut accepted: Option<(f64, Vec<Vec3>)> = None;
        for _ in 0..MAX_BACKTRACKS {
            let mut trial: Vec<Vec3> = saved
                .iter()
                .zip(&dir)
                .map(|(r, d)| *r + *d * alpha)
                .collect();
            if let Some(b) = sys.box_ {
                for p in &mut trial {
                    *p = b.wrap(*p);
                }
            }
            sys.pos = trial;
            let (f_new, e_new) = model.evaluate(sys);
            if e_new <= energy + ARMIJO_C1 * alpha * slope {
                accepted = Some((e_new, f_new));
                break;
            }
            alpha *= 0.5;
        }
        let Some((e_new, f_new)) = accepted else {
            sys.pos = saved;
            if history.is_empty() {
                break;
            }
            // A quasi-Newton direction can be useless across a potential
            // discontinuity (e.g. an unshifted cutoff): restart from steepest descent.
            history.clear();
            continue;
        };

        // Record the curvature pair with the actual displacement and gradient change.
        let s: Vec<Vec3> = dir.iter().map(|d| *d * alpha).collect();
        let new_grad: Vec<Vec3> = f_new.iter().map(|f| -*f).collect();
        let y: Vec<Vec3> = new_grad.iter().zip(&grad).map(|(a, b)| *a - *b).collect();
        if field_dot(&s, &y) > 1e-12 {
            if history.len() == memory {
                history.pop_front();
            }
            history.push_back((s, y));
        }

        energy = e_new;
        sys.force = f_new;
        steps += 1;
    }

    let final_max_force = sys.force.iter().map(|f| f.norm()).fold(0.0f64, f64::max);
    MinimiseReport {
        steps,
        initial_energy,
        final_energy: energy,
        final_max_force,
        converged,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::box3::Box3;
    use tpt_chem_core::forcefield::LennardJones;

    #[test]
    fn relaxes_a_clashing_dimer() {
        // Two argon atoms at 2.5 Å (deep on the repulsive wall, U > 0).
        let mut sys = System::new();
        let lj = LennardJones {
            sigma: 3.4,
            epsilon: 0.997,
        };
        sys.add_atom(lj, 0.0, 39.95, Vec3::new(0.0, 0.0, 0.0));
        sys.add_atom(lj, 0.0, 39.95, Vec3::new(2.5, 0.0, 0.0));
        let model = ForceModel::AllPairs { cutoff: None };
        let report = minimise(&mut sys, &model, 500, 1e-4, 0.05);
        assert!(report.final_energy < 0.0, "must reach the well: {report:?}");
        // σ_min of the LB combination: equilibrium at r = 3.4·2^(1/6).
        let r_eq = 3.4 * f64::powf(2.0, 1.0 / 6.0);
        let r = sys.pos[1].x - sys.pos[0].x;
        assert!((r - r_eq).abs() < 1e-2, "r = {r} vs {r_eq}");
        assert!(report.converged, "{report:?}");
    }

    #[test]
    fn relaxes_jittered_periodic_lattice() {
        let mut sys = System::new();
        let lj = LennardJones {
            sigma: 3.4,
            epsilon: 0.997,
        };
        let mut rng = tpt_chem_core::rng::Rng::new(21);
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    sys.add_atom(
                        lj,
                        0.0,
                        39.95,
                        Vec3::new(
                            i as f64 * 4.5 + rng.range(-0.6, 0.6),
                            j as f64 * 4.5 + rng.range(-0.6, 0.6),
                            k as f64 * 4.5 + rng.range(-0.6, 0.6),
                        ),
                    );
                }
            }
        }
        sys.set_box(Box3::cubic(13.5));
        sys.wrap_in_box();
        let model = ForceModel::AllPairs { cutoff: Some(6.5) };
        let report = minimise(&mut sys, &model, 2000, 1e-3, 0.02);
        assert!(
            report.final_energy < report.initial_energy,
            "energy must decrease: {report:?}"
        );
        // The relaxed force scale must be far below the starting one
        // (clashing jitter has |F| of several kJ/mol/Å).
        assert!(report.final_max_force < 1.0, "{report:?}");
    }

    #[test]
    fn lbfgs_relaxes_a_clashing_dimer() {
        let mut sys = System::new();
        let lj = LennardJones {
            sigma: 3.4,
            epsilon: 0.997,
        };
        sys.add_atom(lj, 0.0, 39.95, Vec3::new(0.0, 0.0, 0.0));
        sys.add_atom(lj, 0.0, 39.95, Vec3::new(2.5, 0.0, 0.0));
        let model = ForceModel::AllPairs { cutoff: None };
        let report = minimise_lbfgs(&mut sys, &model, 200, 1e-4, 0.2, 5);
        assert!(report.converged, "{report:?}");
        let r_eq = 3.4 * f64::powf(2.0, 1.0 / 6.0);
        let r = sys.pos[1].x - sys.pos[0].x;
        assert!((r - r_eq).abs() < 1e-3, "r = {r} vs {r_eq}");
        assert!(report.final_energy < 0.0, "{report:?}");
    }

    #[test]
    fn lbfgs_matches_steepest_descent_on_jittered_cluster() {
        // Non-periodic and uncut-off, so the energy is smooth everywhere.
        // A periodic box with an unshifted LJ cutoff has energy jumps that no
        // gradient method can cross, so it is not a fair benchmark here.
        let build = || {
            let mut sys = System::new();
            let lj = LennardJones {
                sigma: 3.4,
                epsilon: 0.997,
            };
            let mut rng = tpt_chem_core::rng::Rng::new(21);
            for i in 0..3 {
                for j in 0..3 {
                    for k in 0..3 {
                        sys.add_atom(
                            lj,
                            0.0,
                            39.95,
                            Vec3::new(
                                i as f64 * 4.0 + rng.range(-0.6, 0.6),
                                j as f64 * 4.0 + rng.range(-0.6, 0.6),
                                k as f64 * 4.0 + rng.range(-0.6, 0.6),
                            ),
                        );
                    }
                }
            }
            sys
        };
        let model = ForceModel::AllPairs { cutoff: None };

        let mut sd = build();
        let sd_rep = minimise(&mut sd, &model, 5000, 1e-3, 0.02);
        let mut lb = build();
        let lb_rep = minimise_lbfgs(&mut lb, &model, 5000, 1e-3, 0.2, 8);

        assert!(lb_rep.final_energy < lb_rep.initial_energy, "{lb_rep:?}");
        assert!(lb_rep.converged, "{lb_rep:?}");
        assert!(lb_rep.final_max_force < 1e-3, "{lb_rep:?}");
        // Both settle into stationary points; the cluster has several local
        // minima, so the two methods need not land in the same one.
        assert!(sd_rep.converged, "{sd_rep:?}");
        // The quasi-Newton search should need far fewer iterations than steepest descent.
        assert!(
            lb_rep.steps < sd_rep.steps,
            "lbfgs {lb_rep:?} vs sd {sd_rep:?}"
        );
    }

    #[test]
    fn already_minimized_system_is_a_fixed_point() {
        let mut sys = System::new();
        let lj = LennardJones {
            sigma: 3.4,
            epsilon: 0.997,
        };
        let r_eq = 3.4 * f64::powf(2.0, 1.0 / 6.0);
        sys.add_atom(lj, 0.0, 39.95, Vec3::ZERO);
        sys.add_atom(lj, 0.0, 39.95, Vec3::new(r_eq, 0.0, 0.0));
        let model = ForceModel::AllPairs { cutoff: None };
        let report = minimise(&mut sys, &model, 50, 1e-6, 0.05);
        assert!(report.converged);
        assert_eq!(report.steps, 0);
        assert!((report.final_energy - report.initial_energy).abs() < 1e-12);
    }
}
