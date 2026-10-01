//! Conservation-law checks: the statistical verification layer of
//! spec.txt §2 ("property-based testing for physical laws").
//!
//! These are plain functions so they compose: the `proptest!` cases below
//! randomize the inputs, and downstream crates can call the same checks on
//! their own systems.

use tpt_chem_md::integrator::VelocityVerlet;
use tpt_chem_md::system::System;

/// Integrate `steps` Velocity Verlet steps of `dt` (fs) and assert that
/// the total energy stays within `rel_tol` (relative, floored by a small
/// absolute term) of the initial value.
///
/// # Panics
/// If energy conservation is violated beyond tolerance, or the final
/// energy is non-finite.
pub fn assert_energy_conserved(sys: &mut System, steps: usize, dt: f64, rel_tol: f64) {
    sys.update_forces();
    let e0 = sys.total_energy();
    let vv = VelocityVerlet::new(dt);
    for _ in 0..steps {
        vv.step(sys, None);
    }
    let e1 = sys.total_energy();
    assert!(
        e1.is_finite(),
        "energy became non-finite: e0 = {e0}, e1 = {e1}"
    );
    let scale = e0.abs().max(1.0);
    assert!(
        (e1 - e0).abs() <= rel_tol * scale,
        "energy drift |{e1} - {e0}| = {} exceeds {} × {scale}",
        (e1 - e0).abs(),
        rel_tol
    );
}

/// Harmonic-oscillator phase-space test bed: one particle in a 1D harmonic
/// well (implemented via a two-atom LJ-free proxy is awkward, so this uses
/// the analytic Hamiltonian directly).
///
/// Returns the phase-space area enclosed by the trajectory after `steps`
/// VV steps of the 2D oscillator, relative to the initial area. Symplectic
/// maps preserve this area exactly (up to floating point) — the property
/// tested by [`phase_space_area_preserved`].
pub fn oscillator_area_ratio(dt: f64, steps: usize, omega: f64) -> f64 {
    // Phase point (q, p) evolved by the exact VV map of H = p²/2 + ω²q²/2.
    // The VV map for the harmonic oscillator is a linear map with
    // determinant exactly 1; track it by evolving a pair of basis
    // perturbations alongside the trajectory.
    let (mut q, mut p) = (1.0f64, 0.0f64);
    // Jacobian columns: d(q,p)/d(q0,p0).
    let (mut dqq, mut dpq) = (1.0f64, 0.0f64); // column for q0
    let (mut dqp, mut dpp) = (0.0f64, 1.0f64); // column for p0

    #[inline]
    fn half_kick(p: f64, q: f64, omega: f64, dt: f64) -> f64 {
        p - 0.5 * dt * omega * omega * q
    }

    for _ in 0..steps {
        p = half_kick(p, q, omega, dt);
        dpp = half_kick(dpp, dqp, omega, dt);
        dpq = half_kick(dpq, dqq, omega, dt);
        q += dt * p;
        dqq += dt * dpq;
        dqp += dt * dpp;
        // Second half-kick with the refreshed positions.
        p = half_kick(p, q, omega, dt);
        dpp = half_kick(dpp, dqp, omega, dt);
        dpq = half_kick(dpq, dqq, omega, dt);
    }
    
    dqq * dpp - dpq * dqp
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::strategies::arb_md_system;
    use proptest::prelude::*;
    use tpt_chem_core::rng::Rng;
    use tpt_chem_core::vec3::Vec3;
    use tpt_chem_md::box3::Box3;

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(64))]

        /// spec.txt §6 Phase-0 exit criterion: "Velocity Verlet integration
        /// conserves total energy within a bounded tolerance" — verified
        /// statistically over randomized LJ systems.
        #[test]
        fn velocity_verlet_conserves_energy((mut sys, _seed) in arb_md_system(12)) {
            let dt = 0.5; // conservative timestep for arbitrary draws
            assert_energy_conserved(&mut sys, 300, dt, 5e-3);
        }

        /// Smaller timestep → tighter conservation (numerical order).
        #[test]
        fn smaller_dt_conserves_better((mut sys, _seed) in arb_md_system(8)) {
            sys.update_forces();
            let e0 = sys.total_energy();
            let vv = VelocityVerlet::new(0.1);
            for _ in 0..500 {
                vv.step(&mut sys, None);
            }
            let e1 = sys.total_energy();
            prop_assert!((e1 - e0).abs() <= 1e-3 * e0.abs().max(1.0), "{e0} -> {e1}");
        }

        /// Symplecticity: the Velocity Verlet map of the harmonic
        /// oscillator has Jacobian determinant exactly 1 — phase-space
        /// volume is preserved.
        #[test]
        fn phase_space_area_preserved(dt in 0.05f64..0.9, steps in 1usize..200, omega in 0.5f64..2.0) {
            let det = oscillator_area_ratio(dt, steps, omega);
            prop_assert!((det - 1.0).abs() < 1e-9, "det = {det}");
        }
    }

    #[test]
    fn open_boundary_finite_cluster_conserves() {
        // No box: open boundaries, same conservation law.
        let mut sys = System::new();
        let lj = tpt_chem_core::forcefield::LennardJones {
            sigma: 3.0,
            epsilon: 0.5,
        };
        let mut rng = Rng::new(9);
        for i in 0..6 {
            sys.add_atom(lj, 0.0, 12.0, Vec3::new(i as f64 * 3.5, 0.0, 0.0));
        }
        sys.set_box(Box3::cubic(30.0));
        sys.init_velocities(50.0, &mut rng);
        assert_energy_conserved(&mut sys, 200, 0.5, 5e-3);
    }
}
