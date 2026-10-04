//! Energy minimisation: steepest descent with an adaptive line step.
//!
//! Used to relax input geometries before dynamics (jittered lattices,
//! parsed structures with clashes). Forces and energies come from a
//! [`ForceModel`], so the minimised state is a stationary point of exactly
//! the model the subsequent run uses.

use crate::forces::ForceModel;
use crate::system::System;

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::box3::Box3;
    use tpt_chem_core::forcefield::LennardJones;
    use tpt_chem_core::vec3::Vec3;

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
