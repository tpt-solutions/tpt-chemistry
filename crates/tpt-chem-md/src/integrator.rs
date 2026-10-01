//! Symplectic integrators. Phase 1: Velocity Verlet.

use tpt_chem_core::vec3::Vec3;

use crate::neighbors::VerletList;
use crate::system::System;
use crate::ACC_PER_FORCE_OVER_MASS;

/// Velocity Verlet (sweep form):
///
/// 1. `v(t+Δt/2) = v(t) + F(t)/2m·Δt`
/// 2. `r(t+Δt) = r(t) + v(t+Δt/2)·Δt`
/// 3. `F(t+Δt) = F(r(t+Δt))`
/// 4. `v(t+Δt) = v(t+Δt/2) + F(t+Δt)/2m·Δt`
///
/// Symplectic and time-reversible: total energy oscillates tightly around
/// the conserved value instead of drifting (verified statistically in
/// `tpt-chem-verify`).
#[derive(Clone, Copy, Debug)]
pub struct VelocityVerlet {
    /// Time step (fs).
    pub dt: f64,
}

impl VelocityVerlet {
    /// Create an integrator with the given time step (fs).
    pub fn new(dt: f64) -> Self {
        VelocityVerlet { dt }
    }

    /// Advance the system by one step using the all-pairs force kernel.
    ///
    /// `cutoff` selects the truncation (`None` = no truncation); the
    /// potential energy reported uses the same pair set as the forces.
    ///
    /// Returns the potential energy after the step (kJ·mol⁻¹); the forces
    /// stored on the system are the fresh `F(t+Δt)`.
    pub fn step(&self, sys: &mut System, cutoff: Option<f64>) -> f64 {
        let dt = self.dt;
        if sys.force.iter().all(|f| *f == Vec3::ZERO) && sys.potential_energy == 0.0 {
            let (f, e) = crate::forces::forces_all_pairs(sys, sys.box_, cutoff);
            sys.force = f;
            sys.potential_energy = e;
        }
        for i in 0..sys.len() {
            let a = sys.force[i] * (ACC_PER_FORCE_OVER_MASS / sys.mass[i]);
            sys.vel[i] += a * (0.5 * dt);
            sys.pos[i] += sys.vel[i] * dt;
        }
        if let Some(b) = sys.box_ {
            for p in &mut sys.pos {
                *p = b.wrap(*p);
            }
        }
        let (f, e) = crate::forces::forces_all_pairs(sys, sys.box_, cutoff);
        sys.force = f;
        sys.potential_energy = e;
        for i in 0..sys.len() {
            let a = sys.force[i] * (ACC_PER_FORCE_OVER_MASS / sys.mass[i]);
            sys.vel[i] += a * (0.5 * dt);
        }
        sys.potential_energy
    }

    /// Advance the system by one step using a Verlet neighbor list
    /// (rebuilding it when the skin criterion trips).
    ///
    /// Returns the potential energy after the step (kJ·mol⁻¹).
    pub fn step_with_neighbors(&self, sys: &mut System, neighbors: &mut VerletList) -> f64 {
        let dt = self.dt;
        if sys.force.iter().all(|f| *f == Vec3::ZERO) && sys.potential_energy == 0.0 {
            sys.update_forces();
        }
        if neighbors.needs_update(&sys.pos, sys.box_) {
            neighbors.update(&sys.pos, sys.box_);
        }
        for i in 0..sys.len() {
            let a = sys.force[i] * (ACC_PER_FORCE_OVER_MASS / sys.mass[i]);
            sys.vel[i] += a * (0.5 * dt);
            sys.pos[i] += sys.vel[i] * dt;
        }
        if let Some(b) = sys.box_ {
            for p in &mut sys.pos {
                *p = b.wrap(*p);
            }
        }
        let (f, e) = crate::forces::forces_neighbor_list(sys, sys.box_, neighbors);
        sys.force = f;
        sys.potential_energy = e;
        for i in 0..sys.len() {
            let a = sys.force[i] * (ACC_PER_FORCE_OVER_MASS / sys.mass[i]);
            sys.vel[i] += a * (0.5 * dt);
        }
        sys.potential_energy
    }

    /// Run `steps` steps with the all-pairs kernel; returns the total
    /// energy (KE + final PE) in kJ·mol⁻¹.
    pub fn run(&self, sys: &mut System, steps: usize, cutoff: Option<f64>) -> f64 {
        for _ in 0..steps {
            self.step(sys, cutoff);
        }
        sys.total_energy()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::box3::Box3;
    use tpt_chem_core::forcefield::LennardJones;
    use tpt_chem_core::rng::Rng;

    fn argon_cluster(n: usize, seed: u64) -> System {
        // Jittered simple-cubic lattice with ~4.5 Å spacing: no atom
        // overlaps, so the LJ dynamics start in a physically sane state.
        let mut sys = System::new();
        let lj = LennardJones {
            sigma: 3.4,
            epsilon: 0.997,
        };
        let mut rng = Rng::new(seed);
        let per_axis = (n as f64).cbrt().ceil() as usize;
        let spacing = 4.5;
        let l = per_axis as f64 * spacing + 4.0;
        let b = Box3::cubic(l);
        let mut placed = 0;
        for i in 0..per_axis {
            for j in 0..per_axis {
                for k in 0..per_axis {
                    if placed >= n {
                        break;
                    }
                    let jit = |rng: &mut Rng| rng.range(-0.5, 0.5);
                    let p = Vec3::new(
                        2.0 + i as f64 * spacing + jit(&mut rng),
                        2.0 + j as f64 * spacing + jit(&mut rng),
                        2.0 + k as f64 * spacing + jit(&mut rng),
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
    fn energy_conservation_over_short_run() {
        let mut sys = argon_cluster(20, 11);
        sys.update_forces();
        let e0 = sys.total_energy();
        let vv = VelocityVerlet::new(1.0);
        for _ in 0..500 {
            vv.step(&mut sys, None);
        }
        let e1 = sys.total_energy();
        // Relative drift over 0.5 ps must stay tiny for a symplectic
        // integrator at 1 fs.
        assert!((e1 - e0).abs() / e0.abs() < 1e-3, "drift: {e0} -> {e1}");
    }

    #[test]
    fn neighbor_list_path_matches_all_pairs() {
        let mut sys = argon_cluster(15, 21);
        sys.update_forces();
        let mut vl = VerletList::build(&sys.pos, 8.0, 1.5, sys.box_);
        let vv = VelocityVerlet::new(1.0);
        let mut a = sys.clone();
        let mut b = sys.clone();

        // One step: the two force kernels are identical up to summation
        // order, so trajectories agree to floating-point noise.
        vv.step(&mut a, Some(8.0));
        vv.step_with_neighbors(&mut b, &mut vl);
        let d1 = a
            .pos
            .iter()
            .zip(b.pos.iter())
            .map(|(p, q)| p.dist(*q))
            .fold(0.0f64, f64::max);
        assert!(d1 < 1e-10, "single-step position difference {d1}");

        // Long run: Lyapunov amplification of the fp-ordering noise makes
        // exact trajectory comparison meaningless, so compare energies
        // statistically. Both kernels must use the same cutoff: `step` is
        // uncutoffed, so the reference energy is re-evaluated at 8 Å.
        for _ in 0..100 {
            vv.step(&mut a, Some(8.0));
            vv.step_with_neighbors(&mut b, &mut vl);
        }
        let (_, ea_pot) = crate::forces::forces_all_pairs(&a, a.box_, Some(8.0));
        let ea = ea_pot + a.kinetic_energy();
        let eb = b.potential_energy + b.kinetic_energy();
        assert!(
            (ea - eb).abs() / eb.abs() < 1e-3,
            "total energies diverged: {ea} vs {eb}"
        );
    }

    #[test]
    fn time_reversibility() {
        // Reverse velocities, run back, recover the initial state.
        let mut sys = argon_cluster(10, 5);
        sys.update_forces();
        let pos0 = sys.pos.clone();
        let vv = VelocityVerlet::new(1.0);
        for _ in 0..100 {
            vv.step(&mut sys, None);
        }
        for v in &mut sys.vel {
            *v = -*v;
        }
        for _ in 0..100 {
            vv.step(&mut sys, None);
        }
        let drift: f64 = pos0
            .iter()
            .zip(sys.pos.iter())
            .map(|(a, b)| a.dist(*b))
            .fold(0.0f64, f64::max);
        // Time reversibility is exact up to floating point.
        assert!(drift < 1e-8, "reversibility drift {drift}");
    }
}
