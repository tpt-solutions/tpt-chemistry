//! The MD system state: per-atom positions, velocities, forces, and
//! parameters, in the MD unit system (Å, fs, amu, kJ·mol⁻¹).

use tpt_chem_core::forcefield::LennardJones;
use tpt_chem_core::rng::Rng;
use tpt_chem_core::vec3::Vec3;

use crate::box3::Box3;
use crate::forces;
use crate::KJ_MOL_PER_MDU;

/// A classical MD system (structure-of-arrays layout).
#[derive(Clone, Debug, Default)]
pub struct System {
    /// Positions (Å).
    pub pos: Vec<Vec3>,
    /// Velocities (Å/fs).
    pub vel: Vec<Vec3>,
    /// Forces (kJ·mol⁻¹·Å⁻¹), refreshed by the integrator.
    pub force: Vec<Vec3>,
    /// Lennard-Jones parameters per atom.
    pub lj: Vec<LennardJones>,
    /// Partial charges (e).
    pub charge: Vec<f64>,
    /// Masses (amu).
    pub mass: Vec<f64>,
    /// Optional atom names.
    pub names: Vec<String>,
    /// Optional periodic box (Å).
    pub box_: Option<Box3>,
    /// Optional bonded topology; when present, every force evaluation adds
    /// the bonded energy/forces on top of the nonbonded pair terms.
    pub bonded: Option<crate::bonded::Bonded>,
    /// Optional nonbonded exclusions / 1-4 scaling applied by every pair
    /// kernel (see [`crate::exclusions`]).
    pub exclusions: Option<crate::exclusions::Exclusions>,
    /// Last potential energy reported by a force evaluation (kJ·mol⁻¹).
    pub potential_energy: f64,
}

impl System {
    /// An empty system.
    pub fn new() -> Self {
        Self::default()
    }

    /// Number of atoms.
    pub fn len(&self) -> usize {
        self.pos.len()
    }

    /// True if the system has no atoms.
    pub fn is_empty(&self) -> bool {
        self.pos.is_empty()
    }

    /// Add an atom with Lennard-Jones parameters, charge (e), mass (amu),
    /// and position (Å).
    pub fn add_atom(&mut self, lj: LennardJones, charge: f64, mass: f64, pos: Vec3) -> usize {
        let idx = self.pos.len();
        self.pos.push(pos);
        self.vel.push(Vec3::ZERO);
        self.force.push(Vec3::ZERO);
        self.lj.push(lj);
        self.charge.push(charge);
        self.mass.push(mass);
        self.names.push(String::new());
        idx
    }

    /// Set the periodic box.
    pub fn set_box(&mut self, box_: Box3) {
        self.box_ = Some(box_);
    }

    /// Wrap all positions into the periodic box (no-op without one).
    pub fn wrap_in_box(&mut self) {
        if let Some(b) = self.box_ {
            for p in &mut self.pos {
                *p = b.wrap(*p);
            }
        }
    }

    /// Draw Maxwell–Boltzmann velocities for `temperature` (K), remove the
    /// center-of-mass motion, and return the realized temperature (K).
    ///
    /// The velocity scale uses `√(k_B·T/m)` per component with the exact
    /// `k_B` in MD units (`8.314462618e-3 kJ·mol⁻¹·K⁻¹` → per-amu factor).
    pub fn init_velocities(&mut self, temperature: f64, rng: &mut Rng) -> f64 {
        // k_B per amu in Å²·fs⁻²·K⁻¹: (k_B[kJ/mol/K] / 1e4) gives amu·Å²/fs²/K.
        let kb_per_amu = tpt_chem_core::units::BOLTZMANN_KJ_PER_MOL_K / KJ_MOL_PER_MDU;
        for i in 0..self.len() {
            let sigma = (kb_per_amu * temperature / self.mass[i]).sqrt();
            self.vel[i] = Vec3::new(
                rng.normal(0.0, sigma),
                rng.normal(0.0, sigma),
                rng.normal(0.0, sigma),
            );
        }
        self.remove_com_motion();
        self.temperature()
    }

    /// Subtract the center-of-mass velocity from every atom.
    pub fn remove_com_motion(&mut self) {
        let total_mass: f64 = self.mass.iter().sum();
        if total_mass <= 0.0 || self.is_empty() {
            return;
        }
        let mut com = Vec3::ZERO;
        for i in 0..self.len() {
            com += self.vel[i] * self.mass[i];
        }
        com = com / total_mass;
        for v in &mut self.vel {
            *v -= com;
        }
    }

    /// Kinetic energy in kJ·mol⁻¹ (`Σ ½mv²` in MD units × 10⁴).
    pub fn kinetic_energy(&self) -> f64 {
        let mut ke = 0.0;
        for i in 0..self.len() {
            ke += 0.5 * self.mass[i] * self.vel[i].norm_sq();
        }
        ke * KJ_MOL_PER_MDU
    }

    /// Instantaneous temperature in K from `2·KE/(dof·k_B)` with
    /// `dof = 3N − 3` (center-of-mass constrained).
    pub fn temperature(&self) -> f64 {
        let dof = (3 * self.len()) as f64 - 3.0;
        if dof <= 0.0 {
            return 0.0;
        }
        2.0 * self.kinetic_energy() / (dof * tpt_chem_core::units::BOLTZMANN_KJ_PER_MOL_K)
    }

    /// Add the bonded topology's forces and energy (when present) to a
    /// nonbonded `(forces, energy)` result.
    pub fn with_bonded(&self, (mut f, mut e): (Vec<Vec3>, f64)) -> (Vec<Vec3>, f64) {
        if let Some(b) = &self.bonded {
            let (fb, eb) = crate::bonded::bonded_energy_forces(b, &self.pos);
            for (x, y) in f.iter_mut().zip(&fb) {
                *x += *y;
            }
            e += eb;
        }
        (f, e)
    }

    /// Recompute forces and potential energy (all-pairs reference path)
    /// and store the energy in [`System::potential_energy`].
    pub fn update_forces(&mut self) {
        let (f, e) = self.with_bonded(forces::forces_all_pairs(self, self.box_, None));
        self.force = f;
        self.potential_energy = e;
    }

    /// Total energy (kinetic + last potential) in kJ·mol⁻¹.
    pub fn total_energy(&self) -> f64 {
        self.kinetic_energy() + self.potential_energy
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::KJ_MOL_PER_MDU;

    #[test]
    fn md_unit_factor_is_exact_relation() {
        // amu·(1e5 m/s)²·N_A / 1000 = 1e4 kJ/mol (to 1e-4 relative).
        let expect = tpt_chem_core::units::AMU_KG * 1e10 * tpt_chem_core::units::AVOGADRO / 1000.0;
        assert!((KJ_MOL_PER_MDU - expect).abs() < 1e-8 * expect);
    }

    #[test]
    fn thermal_temperature_of_ar_gas() {
        let mut sys = System::new();
        let lj = LennardJones {
            sigma: 3.4,
            epsilon: 0.997,
        };
        let mut rng = Rng::new(3);
        for i in 0..50 {
            sys.add_atom(lj, 0.0, 39.95, Vec3::new(i as f64 * 6.0, 0.0, 0.0));
        }
        let t = sys.init_velocities(300.0, &mut rng);
        // Fluctuation at N=50 is √(2/dof) ≈ 19%; allow 40%.
        assert!((t - 300.0).abs() < 120.0, "T = {t}");
        // COM velocity removed.
        let mut p = Vec3::ZERO;
        for i in 0..sys.len() {
            p += sys.vel[i] * sys.mass[i];
        }
        assert!(p.norm() < 1e-10);
    }

    #[test]
    fn kinetic_energy_units() {
        // A single 1 amu atom moving at 0.1 Å/fs: KE = ½·1·0.01 = 5e-3 MDU
        // = 5e-3 × 1e4 = 50 kJ/mol... but 3N-3 = 0 dof, so T = 0.
        let mut sys = System::new();
        sys.add_atom(
            LennardJones {
                sigma: 3.0,
                epsilon: 0.1,
            },
            0.0,
            1.0,
            Vec3::ZERO,
        );
        sys.vel[0] = Vec3::new(0.1, 0.0, 0.0);
        assert!((sys.kinetic_energy() - 50.0).abs() < 1e-9);
        assert_eq!(sys.temperature(), 0.0);
    }
}
