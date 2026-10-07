//! Holonomic distance constraints: SHAKE (positions) and RATTLE
//! (velocities) wrapped around velocity Verlet.
//!
//! Constrained dynamics allow rigid bonds (e.g. X-H, rigid water) and so
//! larger time steps. Degrees of freedom removed by the constraints are
//! *not* subtracted by [`System::temperature`]; scale by
//! `3N / (3N - n_constraints)` when a constrained temperature is needed.

use tpt_chem_core::vec3::Vec3;

use crate::forces::ForceModel;
use crate::system::System;
use crate::ACC_PER_FORCE_OVER_MASS;

/// A fixed distance `d` (Å) between atoms `i` and `j`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DistanceConstraint {
    /// First atom.
    pub i: usize,
    /// Second atom.
    pub j: usize,
    /// Target separation (Å).
    pub d: f64,
}

/// Failure to satisfy the constraints within the iteration cap.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ConstraintError;

impl core::fmt::Display for ConstraintError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "SHAKE/RATTLE did not converge")
    }
}

impl std::error::Error for ConstraintError {}

/// A set of distance constraints solved iteratively.
#[derive(Clone, Debug)]
pub struct Constraints {
    /// The constraints.
    pub list: Vec<DistanceConstraint>,
    /// Relative tolerance on `|r² - d²| / d²`.
    pub tolerance: f64,
    /// Maximum sweeps.
    pub max_iter: usize,
}

impl Constraints {
    /// A constraint set with tolerance `1e-10` and 500 sweeps.
    pub fn new(list: Vec<DistanceConstraint>) -> Self {
        Constraints {
            list,
            tolerance: 1e-10,
            max_iter: 500,
        }
    }

    /// Constraints taken from the current geometry for each `(i, j)` pair.
    pub fn from_geometry(sys: &System, pairs: &[(usize, usize)]) -> Self {
        let list = pairs
            .iter()
            .map(|&(i, j)| DistanceConstraint {
                i,
                j,
                d: sep(sys, sys.pos[i], sys.pos[j]).norm(),
            })
            .collect();
        Constraints::new(list)
    }

    /// Number of constraints.
    pub fn len(&self) -> usize {
        self.list.len()
    }

    /// True when there are no constraints.
    pub fn is_empty(&self) -> bool {
        self.list.is_empty()
    }

    /// SHAKE: correct the unconstrained positions in `sys.pos` so every
    /// constraint holds, moving atoms along the *reference* bond vectors
    /// from `old_pos`. Velocities are corrected by `Δr / dt`.
    ///
    /// # Errors
    /// [`ConstraintError`] if the iteration cap is reached.
    pub fn shake(
        &self,
        sys: &mut System,
        old_pos: &[Vec3],
        dt: f64,
    ) -> Result<(), ConstraintError> {
        for _ in 0..self.max_iter {
            let mut done = true;
            for c in &self.list {
                let s = sep(sys, sys.pos[c.i], sys.pos[c.j]);
                let diff = s.dot(s) - c.d * c.d;
                if diff.abs() > self.tolerance * c.d * c.d {
                    done = false;
                    let s0 = sep(sys, old_pos[c.i], old_pos[c.j]);
                    let (wi, wj) = (1.0 / sys.mass[c.i], 1.0 / sys.mass[c.j]);
                    let g = diff / (2.0 * (wi + wj) * s.dot(s0));
                    let di = s0 * (-g * wi);
                    let dj = s0 * (g * wj);
                    sys.pos[c.i] += di;
                    sys.pos[c.j] += dj;
                    sys.vel[c.i] += di * (1.0 / dt);
                    sys.vel[c.j] += dj * (1.0 / dt);
                }
            }
            if done {
                return Ok(());
            }
        }
        Err(ConstraintError)
    }

    /// RATTLE velocity stage: remove velocity components along each
    /// constrained bond.
    ///
    /// # Errors
    /// [`ConstraintError`] if the iteration cap is reached.
    pub fn rattle_velocities(&self, sys: &mut System) -> Result<(), ConstraintError> {
        for _ in 0..self.max_iter {
            let mut done = true;
            for c in &self.list {
                let r = sep(sys, sys.pos[c.i], sys.pos[c.j]);
                let v = sys.vel[c.i] - sys.vel[c.j];
                let (wi, wj) = (1.0 / sys.mass[c.i], 1.0 / sys.mass[c.j]);
                let k = r.dot(v) / ((wi + wj) * c.d * c.d);
                if (r.dot(v) / (c.d * c.d)).abs() > self.tolerance {
                    done = false;
                    sys.vel[c.i] -= r * (k * wi);
                    sys.vel[c.j] += r * (k * wj);
                }
            }
            if done {
                return Ok(());
            }
        }
        Err(ConstraintError)
    }

    /// Largest relative bond-length violation `|r - d| / d`.
    pub fn max_violation(&self, sys: &System) -> f64 {
        self.list
            .iter()
            .map(|c| (sep(sys, sys.pos[c.i], sys.pos[c.j]).norm() - c.d).abs() / c.d)
            .fold(0.0, f64::max)
    }

    /// One constrained velocity-Verlet step (SHAKE + RATTLE) with `dt` in
    /// fs. Returns the potential energy after the step.
    ///
    /// # Errors
    /// [`ConstraintError`] if SHAKE or RATTLE fails to converge.
    pub fn step(
        &self,
        sys: &mut System,
        model: &ForceModel,
        dt: f64,
    ) -> Result<f64, ConstraintError> {
        if sys.force.len() != sys.len() || sys.force.iter().all(|f| *f == Vec3::ZERO) {
            let (f, e) = model.evaluate(sys);
            sys.force = f;
            sys.potential_energy = e;
        }
        let old_pos = sys.pos.clone();
        for i in 0..sys.len() {
            let a = sys.force[i] * (ACC_PER_FORCE_OVER_MASS / sys.mass[i]);
            sys.vel[i] += a * (0.5 * dt);
            sys.pos[i] += sys.vel[i] * dt;
        }
        self.shake(sys, &old_pos, dt)?;
        if let Some(b) = sys.box_ {
            for p in &mut sys.pos {
                *p = b.wrap(*p);
            }
        }
        let (f, e) = model.evaluate(sys);
        sys.force = f;
        sys.potential_energy = e;
        for i in 0..sys.len() {
            let a = sys.force[i] * (ACC_PER_FORCE_OVER_MASS / sys.mass[i]);
            sys.vel[i] += a * (0.5 * dt);
        }
        self.rattle_velocities(sys)?;
        Ok(e)
    }
}

fn sep(sys: &System, a: Vec3, b: Vec3) -> Vec3 {
    let d = a - b;
    match sys.box_ {
        Some(bx) => bx.min_image(d),
        None => d,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::box3::Box3;
    use tpt_chem_core::forcefield::LennardJones;
    use tpt_chem_core::rng::Rng;

    /// Rigid LJ dimers on a loose lattice in a periodic box.
    fn dimers() -> (System, Constraints) {
        let mut sys = System::new();
        let lj = LennardJones {
            sigma: 3.0,
            epsilon: 0.8,
        };
        // Non-interacting light atoms: no intramolecular LJ clash (the crate
        // has no exclusion lists yet).
        let lj_h = LennardJones {
            sigma: 1.0,
            epsilon: 0.0,
        };
        let mut pairs = Vec::new();
        for ix in 0..3 {
            for iy in 0..3 {
                for iz in 0..2 {
                    let o = Vec3::new(ix as f64 * 6.0, iy as f64 * 6.0, iz as f64 * 6.0);
                    let a = sys.add_atom(lj, 0.0, 16.0, o + Vec3::new(0.5, 0.5, 0.5));
                    let b = sys.add_atom(lj_h, 0.0, 1.0, o + Vec3::new(1.5, 0.5, 0.9));
                    pairs.push((a, b));
                }
            }
        }
        sys.set_box(Box3::cubic(18.0));
        let c = Constraints::from_geometry(&sys, &pairs);
        (sys, c)
    }

    #[test]
    fn bonds_stay_rigid_and_energy_is_conserved() {
        let (mut sys, cons) = dimers();
        let mut rng = Rng::new(21);
        sys.init_velocities(150.0, &mut rng);
        sys.remove_com_motion();
        cons.rattle_velocities(&mut sys).unwrap();
        let model = ForceModel::AllPairs { cutoff: Some(8.0) };
        let (f, e) = model.evaluate(&sys);
        sys.force = f;
        sys.potential_energy = e;
        let e0 = sys.total_energy();
        let mut worst = 0.0f64;
        for _ in 0..2000 {
            cons.step(&mut sys, &model, 1.0).unwrap();
            worst = worst.max(cons.max_violation(&sys));
        }
        assert!(worst < 1e-8, "bond violation {worst}");
        let drift = (sys.total_energy() - e0).abs();
        assert!(drift < 0.05 * e0.abs().max(10.0), "drift {drift}");
    }

    #[test]
    fn rattle_leaves_no_velocity_along_bonds() {
        let (mut sys, cons) = dimers();
        let mut rng = Rng::new(5);
        sys.init_velocities(300.0, &mut rng);
        cons.rattle_velocities(&mut sys).unwrap();
        for c in &cons.list {
            let r = sys.pos[c.i] - sys.pos[c.j];
            let v = sys.vel[c.i] - sys.vel[c.j];
            assert!(r.dot(v).abs() < 1e-8);
        }
    }
}
