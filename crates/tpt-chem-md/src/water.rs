//! Rigid three-site water models (SPC, TIP3P) and a lattice water-box
//! builder that wires up everything a rigid-water simulation needs:
//! charges and LJ sites, intramolecular [`Exclusions`], and SHAKE
//! [`Constraints`] (two O-H bonds plus the H-H distance).
//!
//! Parameters: SPC (Berendsen et al. 1981) and TIP3P (Jorgensen et al.
//! 1983), in Angstrom, kJ/mol and elementary charges.

use tpt_chem_core::forcefield::LennardJones;
use tpt_chem_core::vec3::Vec3;

use crate::box3::Box3;
use crate::constraints::{Constraints, DistanceConstraint};
use crate::exclusions::{Exclusions, PairScale};
use crate::system::System;

/// A rigid 3-site water model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WaterModel {
    /// Name.
    pub name: &'static str,
    /// Oxygen LJ sigma (Angstrom).
    pub sigma_o: f64,
    /// Oxygen LJ epsilon (kJ/mol).
    pub epsilon_o: f64,
    /// Oxygen charge (e); hydrogens carry `-q_o / 2`.
    pub q_o: f64,
    /// O-H bond length (Angstrom).
    pub r_oh: f64,
    /// H-O-H angle (degrees).
    pub theta_hoh: f64,
}

/// SPC water.
pub const SPC: WaterModel = WaterModel {
    name: "SPC",
    sigma_o: 3.1656,
    epsilon_o: 0.650,
    q_o: -0.82,
    r_oh: 1.0,
    theta_hoh: 109.47,
};

/// TIP3P water.
pub const TIP3P: WaterModel = WaterModel {
    name: "TIP3P",
    sigma_o: 3.15061,
    epsilon_o: 0.6364,
    q_o: -0.834,
    r_oh: 0.9572,
    theta_hoh: 104.52,
};

const MASS_O: f64 = 15.9994;
const MASS_H: f64 = 1.008;
/// Volume per molecule of liquid water at ~0.997 g/cm^3 (Angstrom^3).
const VOLUME_PER_MOLECULE: f64 = 29.9;

impl WaterModel {
    /// H-H distance implied by the geometry (Angstrom).
    pub fn r_hh(&self) -> f64 {
        2.0 * self.r_oh * (self.theta_hoh.to_radians() / 2.0).sin()
    }

    /// Site offsets from the oxygen: `[O, H1, H2]` in the molecular frame
    /// (C2 axis along +z, molecule in the xz plane).
    pub fn geometry(&self) -> [Vec3; 3] {
        let half = self.theta_hoh.to_radians() / 2.0;
        [
            Vec3::ZERO,
            Vec3::new(self.r_oh * half.sin(), 0.0, self.r_oh * half.cos()),
            Vec3::new(-self.r_oh * half.sin(), 0.0, self.r_oh * half.cos()),
        ]
    }
}

/// A ready-to-run rigid-water system.
#[derive(Clone, Debug)]
pub struct WaterBox {
    /// The particles (atoms ordered O, H, H per molecule). `exclusions` is
    /// already set.
    pub system: System,
    /// SHAKE constraints holding every molecule rigid.
    pub constraints: Constraints,
    /// Number of molecules.
    pub n_molecules: usize,
}

/// Build `n_per_axis^3` molecules on a simple-cubic lattice at liquid
/// density in a periodic cubic box. All molecules share one orientation;
/// follow with equilibration (thermostatted constrained dynamics) before
/// measuring anything.
///
/// # Panics
///
/// Panics if `n_per_axis` is zero.
pub fn water_box(model: &WaterModel, n_per_axis: usize) -> WaterBox {
    assert!(
        n_per_axis > 0,
        "water_box: need at least one molecule per axis"
    );
    let spacing = VOLUME_PER_MOLECULE.cbrt();
    let l = spacing * n_per_axis as f64;
    let geom = model.geometry();
    let lj_o = LennardJones {
        sigma: model.sigma_o,
        epsilon: model.epsilon_o,
    };
    let lj_h = LennardJones {
        sigma: 1.0,
        epsilon: 0.0,
    };
    let q_h = -model.q_o / 2.0;

    let mut sys = System::new();
    let mut ex = Exclusions::new();
    let mut cons = Vec::new();
    let full = PairScale {
        lj: 0.0,
        coulomb: 0.0,
    };
    let mut n_mol = 0;
    for ix in 0..n_per_axis {
        for iy in 0..n_per_axis {
            for iz in 0..n_per_axis {
                let o = Vec3::new(ix as f64, iy as f64, iz as f64) * spacing
                    + Vec3::new(0.5, 0.5, 0.5) * (spacing * 0.5);
                let a = sys.add_atom(lj_o, model.q_o, MASS_O, o + geom[0]);
                let h1 = sys.add_atom(lj_h, q_h, MASS_H, o + geom[1]);
                let h2 = sys.add_atom(lj_h, q_h, MASS_H, o + geom[2]);
                for (i, j) in [(a, h1), (a, h2), (h1, h2)] {
                    ex.set(i, j, full);
                }
                cons.push(DistanceConstraint {
                    i: a,
                    j: h1,
                    d: model.r_oh,
                });
                cons.push(DistanceConstraint {
                    i: a,
                    j: h2,
                    d: model.r_oh,
                });
                cons.push(DistanceConstraint {
                    i: h1,
                    j: h2,
                    d: model.r_hh(),
                });
                n_mol += 1;
            }
        }
    }
    sys.set_box(Box3::cubic(l));
    sys.wrap_in_box();
    sys.exclusions = Some(ex);
    WaterBox {
        system: sys,
        constraints: Constraints::new(cons),
        n_molecules: n_mol,
    }
}

/// Temperature (K) of a rigid-water system, using `3N - n_constraints - 3`
/// degrees of freedom.
pub fn constrained_temperature(wb: &WaterBox) -> f64 {
    let n = wb.system.len() as f64;
    let dof = 3.0 * n - wb.constraints.len() as f64 - 3.0;
    2.0 * wb.system.kinetic_energy() / (dof * tpt_chem_core::units::BOLTZMANN_KJ_PER_MOL_K)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forces::ForceModel;
    use tpt_chem_core::rng::Rng;

    #[test]
    fn models_have_consistent_geometry() {
        for m in [SPC, TIP3P] {
            let g = m.geometry();
            assert!(((g[1] - g[0]).norm() - m.r_oh).abs() < 1e-12);
            assert!(((g[1] - g[2]).norm() - m.r_hh()).abs() < 1e-12);
        }
        // Textbook TIP3P H-H distance: 1.5139 A.
        assert!((TIP3P.r_hh() - 1.5139).abs() < 1e-3);
    }

    #[test]
    fn box_is_neutral_dense_and_rigid() {
        let wb = water_box(&TIP3P, 3);
        assert_eq!(wb.n_molecules, 27);
        assert_eq!(wb.system.len(), 81);
        assert_eq!(wb.constraints.len(), 81);
        let q: f64 = wb.system.charge.iter().sum();
        assert!(q.abs() < 1e-12);
        let mass: f64 = wb.system.mass.iter().sum();
        let vol = wb.system.box_.unwrap().volume();
        // g/cm^3 = amu / A^3 * 1.66054.
        let rho = mass / vol * 1.660_539;
        assert!((rho - 0.997).abs() < 0.01, "{rho}");
        assert!(wb.constraints.max_violation(&wb.system) < 1e-12);
        assert_eq!(wb.system.exclusions.as_ref().unwrap().len(), 27 * 3);
    }

    #[test]
    fn rigid_water_dynamics_keeps_constraints() {
        let mut wb = water_box(&TIP3P, 3);
        let mut rng = Rng::new(11);
        wb.system.init_velocities(150.0, &mut rng);
        wb.constraints.rattle_velocities(&mut wb.system).unwrap();
        let model = ForceModel::AllPairs { cutoff: Some(4.0) };
        let mut worst = 0.0f64;
        for _ in 0..50 {
            let e = wb.constraints.step(&mut wb.system, &model, 0.5).unwrap();
            assert!(e.is_finite());
            worst = worst.max(wb.constraints.max_violation(&wb.system));
        }
        assert!(worst < 1e-8, "{worst}");
        let t = constrained_temperature(&wb);
        assert!(t.is_finite() && t > 0.0);
    }
}
