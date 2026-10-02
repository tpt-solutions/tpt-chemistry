//! Thermostats and barostats for MD ensembles.

use crate::system::System;
use tpt_chem_core::num;
use tpt_chem_core::rng::Rng;
use tpt_chem_core::units::BOLTZMANN_KJ_PER_MOL_K;
use tpt_chem_core::vec3::Vec3;

/// Berendsen weak-coupling thermostat: rescales velocities each step so
/// that the temperature approaches `target` with time constant `tau` (fs).
///
/// λ² = 1 + (Δt/τ)(T₀/T − 1), v ← √λ v.
pub fn berendsen_thermostat(sys: &mut System, target_t: f64, tau: f64, dt: f64) {
    let current = sys.temperature();
    if current <= 1e-12 || tau <= 0.0 {
        return;
    }
    let lambda2 = 1.0 + (dt / tau) * (target_t / current - 1.0);
    let lambda = num::sqrt(lambda2.max(0.0));
    for v in &mut sys.vel {
        *v = *v * lambda;
    }
}

/// Berendsen weak-coupling barostat (isotropic): rescales positions and the
/// box toward the target pressure with time constant `tau` (fs).
///
/// The simple form scales fractional coordinates by μ³ where
/// μ = [1 − κ(Δt/τ)(P₀ − P)]^{1/3}, with the compressibility κ folded into
/// the caller-chosen `compressibility` (bar⁻¹); pressures in bar.
pub fn berendsen_barostat(
    sys: &mut System,
    target_p: f64,
    current_p: f64,
    tau: f64,
    dt: f64,
    compressibility: f64,
) {
    if tau <= 0.0 {
        return;
    }
    let mu3 = 1.0 - compressibility * (dt / tau) * (target_p - current_p);
    let mu = num::powf(mu3.max(0.0), 1.0 / 3.0);
    if let Some(b) = sys.box_.as_mut() {
        b.length = b.length * mu;
    }
    for p in &mut sys.pos {
        *p = *p * mu;
    }
}

/// Nosé-Hoover thermostat: evolves an extra thermostat friction variable
/// `xi` alongside the system (single chain node).
///
/// Q (amu·Å²) controls the thermostat mass; `Q = N_dof·k_B·T·τ²` is a
/// typical choice. Returns the updated `xi`.
pub fn nose_hoover_step(
    sys: &mut System,
    xi: &mut f64,
    target_t: f64,
    q_thermostat: f64,
    dt: f64,
) -> f64 {
    let dof = (3 * sys.len()) as f64 - 3.0;
    if dof <= 0.0 {
        return *xi;
    }
    // ξ̇ = (2K − N_dof k_B T)/Q  with K the instantaneous kinetic energy.
    let twice_k = 2.0 * sys.kinetic_energy();
    let force_xi = twice_k - dof * BOLTZMANN_KJ_PER_MOL_K * target_t;
    *xi += dt * force_xi / q_thermostat;
    // Velocity rescale: v̇ = v − ξ v.
    for v in &mut sys.vel {
        *v = *v * (1.0 - *xi * dt);
    }
    *xi
}

/// Langevin thermostat: random friction forces + drag, integrated with the
/// scheme `v ← a·v + b·σ·R`, `x ← x + b·v`, where
/// `a = 1 − γΔt`, `b` pairs with the fluctuation-dissipation relation.
///
/// `gamma` in fs⁻¹, `target_t` in K, `rng` provides standard normals.
pub fn langevin_step(sys: &mut System, gamma: f64, target_t: f64, dt: f64, rng: &mut Rng) {
    // k_B in the MD unit system: Å²·fs⁻²·K⁻¹·amu⁻¹
    // (= 8.314e−3 kJ/mol/K × 1e−4 kJ/mol → amu·Å²/fs²).
    const KB_MD: f64 = 8.314_462_618e-7;
    let a = (-gamma * dt).exp();
    // Fluctuation-dissipation: σ² = k_B T/m (1 − a²).
    for i in 0..sys.len() {
        let sigma = (KB_MD * target_t / sys.mass[i] * (1.0 - a * a)).sqrt();
        let comps = [sys.vel[i].get(0), sys.vel[i].get(1), sys.vel[i].get(2)];
        for (c, comp) in comps.into_iter().enumerate() {
            let noise = rng.normal(0.0, 1.0) * sigma;
            sys.vel[i].set(c, comp * a + noise);
        }
    }
}

/// Leapfrog integrator (equivalent to velocity Verlet, different storage
/// convention: velocities stored at half-steps).
pub struct Leapfrog {
    /// Time step (fs).
    pub dt: f64,
}

impl Leapfrog {
    /// Create a Leapfrog integrator.
    pub fn new(dt: f64) -> Self {
        Leapfrog { dt }
    }

    /// One Leapfrog step using stored forces (must be current).
    /// `v(t+½) = v(t−½) + a(t)Δt; x(t+1) = x(t) + v(t+½)Δt`.
    /// Returns the potential energy after the step.
    pub fn step(&self, sys: &mut System) -> f64 {
        let dt = self.dt;
        if sys.force.iter().all(|f| *f == Vec3::ZERO) && sys.potential_energy == 0.0 {
            sys.update_forces();
        }
        for i in 0..sys.len() {
            let a = sys.force[i] * (crate::ACC_PER_FORCE_OVER_MASS / sys.mass[i]);
            sys.vel[i] += a * dt;
            sys.pos[i] += sys.vel[i] * dt;
        }
        if let Some(b) = sys.box_ {
            for p in &mut sys.pos {
                *p = b.wrap(*p);
            }
        }
        sys.update_forces();
        sys.potential_energy
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::box3::Box3;
    use crate::integrator::VelocityVerlet;
    use crate::system::System;
    use tpt_chem_core::forcefield::LennardJones;
    use tpt_chem_core::rng::Rng;

    fn argon_cluster(n: usize, seed: u64) -> System {
        let mut sys = System::new();
        let lj = LennardJones {
            sigma: 3.4,
            epsilon: 0.997,
        };
        let mut rng = Rng::new(seed);
        let per_axis = (n as f64).cbrt().ceil() as usize;
        let spacing = 4.5;
        let b = Box3::cubic(per_axis as f64 * spacing + 4.0);
        let mut placed = 0;
        for i in 0..per_axis {
            for j in 0..per_axis {
                for k in 0..per_axis {
                    if placed >= n {
                        break;
                    }
                    let p = Vec3::new(
                        2.0 + i as f64 * spacing + rng.range(-0.5, 0.5),
                        2.0 + j as f64 * spacing + rng.range(-0.5, 0.5),
                        2.0 + k as f64 * spacing + rng.range(-0.5, 0.5),
                    );
                    sys.add_atom(lj, 0.0, 39.95, p);
                    placed += 1;
                }
            }
        }
        sys.set_box(b);
        sys.wrap_in_box();
        sys.init_velocities(94.0, &mut rng);
        sys
    }

    #[test]
    fn berendsen_pulls_temperature_to_target() {
        let mut sys = argon_cluster(16, 3);
        sys.init_velocities(300.0, &mut Rng::new(5));
        let vv = VelocityVerlet::new(1.0);
        for _ in 0..2000 {
            vv.step(&mut sys, None);
            berendsen_thermostat(&mut sys, 50.0, 100.0, 1.0);
        }
        let t = sys.temperature();
        assert!((t - 50.0).abs() < 15.0, "T = {t}");
    }

    #[test]
    fn langevin_thermostats_to_target() {
        let mut sys = argon_cluster(16, 4);
        let mut rng = Rng::new(6);
        sys.init_velocities(20.0, &mut rng);
        let vv = VelocityVerlet::new(1.0);
        for _ in 0..5000 {
            vv.step(&mut sys, None);
            langevin_step(&mut sys, 0.01, 120.0, 1.0, &mut rng);
        }
        let t = sys.temperature();
        assert!((t - 120.0).abs() < 40.0, "T = {t}");
    }

    #[test]
    fn nose_hoover_stabilizes() {
        let mut sys = argon_cluster(16, 8);
        let mut xi = 0.0;
        let q = 1e3;
        let vv = VelocityVerlet::new(1.0);
        for _ in 0..20000 {
            vv.step(&mut sys, None);
            xi = nose_hoover_step(&mut sys, &mut xi, 60.0, q, 1.0);
        }
        let t = sys.temperature();
        assert!(t > 20.0 && t < 140.0, "T = {t}");
    }

    #[test]
    fn leapfrog_conserves_energy_like_verlet() {
        let mut sys = argon_cluster(12, 9);
        sys.update_forces();
        let e0 = sys.total_energy();
        let lf = Leapfrog::new(0.5);
        for _ in 0..400 {
            lf.step(&mut sys);
        }
        let e1 = sys.total_energy();
        assert!((e1 - e0).abs() / e0.abs() < 1e-3, "{e0} → {e1}");
    }

    #[test]
    fn leapfrog_matches_verlet_for_one_step() {
        // v₀ stored at t−½ vs t: with forces initialized at t the first
        // Leapfrog step equals velocity Verlet's first step.
        let mut a = argon_cluster(8, 11);
        let mut b = a.clone();
        let dt = 1.0;
        a.update_forces();
        let lf = Leapfrog::new(dt);
        let vv = VelocityVerlet::new(dt);
        lf.step(&mut a);
        vv.step(&mut b, None);
        let diff = a
            .pos
            .iter()
            .zip(b.pos.iter())
            .map(|(p, q)| p.dist(*q))
            .fold(0.0f64, f64::max);
        // First-step drift differs by one acceleration·Δt term
        // (Leapfrog drifts with v(t+½), Verlet with v(t+½) after the same
        // kick — the Δ(a·Δt²) offset is O(dt²) and negligible).
        assert!(diff < 2e-6, "max |Δx| = {diff}");
    }
}
