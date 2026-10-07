//! Topology builder: turn a [`Molecule`] (with or without bonds) into a
//! runnable [`System`] — bond perception from geometry, geometry-referenced
//! bonded terms, element-based LJ parameters, and exclusions with 1-4
//! scaling.
//!
//! This is a *generic* force field for relaxing and sampling around a given
//! structure, not a transferable parameter set like OPLS or AMBER: bond
//! lengths and angles take their equilibrium values from the input
//! geometry (so the input is by construction a minimum of the bonded
//! terms), while force constants come from a small table keyed by bond
//! order and the presence of hydrogen. Atomic charges are taken from the
//! molecule (`Atom::charge`, zero by default). No torsions are generated.

use tpt_chem_core::forcefield::{HarmonicAngle, HarmonicBond, LennardJones};
use tpt_chem_core::molecule::{BondOrder, Molecule};

use crate::bonded::Bonded;
use crate::exclusions::Exclusions;
use crate::system::System;

/// Tunable generic parameters.
#[derive(Clone, Copy, Debug)]
pub struct TopologySettings {
    /// Perceive bonds from geometry when the molecule has none; the
    /// covalent-radius tolerance passed to [`Molecule::infer_bonds`].
    pub bond_tolerance: f64,
    /// Bond force constant for X-H bonds (kJ/mol/A^2).
    pub k_bond_xh: f64,
    /// Bond force constant for single heavy-atom bonds (kJ/mol/A^2).
    pub k_bond_single: f64,
    /// Angle force constant (kJ/mol/rad^2).
    pub k_angle: f64,
    /// 1-4 Lennard-Jones scale.
    pub fudge_lj: f64,
    /// 1-4 Coulomb scale.
    pub fudge_qq: f64,
}

impl Default for TopologySettings {
    fn default() -> Self {
        TopologySettings {
            bond_tolerance: 1.25,
            k_bond_xh: 2900.0,
            k_bond_single: 2200.0,
            k_angle: 400.0,
            fudge_lj: 0.5,
            fudge_qq: 0.5,
        }
    }
}

/// Element-based Lennard-Jones parameters `(sigma A, epsilon kJ/mol)`,
/// OPLS-AA-like values for common elements with a generic fallback.
pub fn element_lj(z: u8) -> LennardJones {
    let (sigma, epsilon) = match z {
        1 => (2.50, 0.126),
        2 => (2.64, 0.084),
        6 => (3.50, 0.276),
        7 => (3.25, 0.711),
        8 => (3.12, 0.711),
        9 => (2.95, 0.255),
        15 => (3.74, 0.838),
        16 => (3.55, 1.046),
        17 => (3.50, 1.109),
        35 => (3.60, 1.50),
        _ => (3.50, 0.300),
    };
    LennardJones { sigma, epsilon }
}

/// Build a [`System`] from a molecule. Positions are the molecule's
/// (Angstrom); the system has no periodic box.
pub fn build_system(mol: &Molecule, settings: &TopologySettings) -> System {
    let mut mol = mol.clone();
    if mol.n_bonds() == 0 {
        mol.infer_bonds(settings.bond_tolerance);
    }
    let n = mol.len();
    let mut sys = System::new();
    for a in mol.atoms() {
        sys.add_atom(element_lj(a.z), a.charge, a.mass(), a.pos);
    }

    let mut bonded = Bonded::new();
    let mut neighbours: Vec<Vec<usize>> = vec![Vec::new(); n];
    for b in mol.bonds() {
        let (i, j) = (b.a.index(), b.b.index());
        neighbours[i].push(j);
        neighbours[j].push(i);
        let r0 = sys.pos[i].dist(sys.pos[j]);
        let has_h = mol.atom(b.a).z == 1 || mol.atom(b.b).z == 1;
        let order_scale = match b.order {
            BondOrder::Single => 1.0,
            BondOrder::Double => 1.8,
            BondOrder::Triple => 2.6,
            _ => 1.4,
        };
        let k = if has_h {
            settings.k_bond_xh
        } else {
            settings.k_bond_single * order_scale
        };
        bonded.bond(i, j, HarmonicBond { k, r0 });
    }
    for (center, nb) in neighbours.iter().enumerate() {
        for x in 0..nb.len() {
            for y in (x + 1)..nb.len() {
                let (a, c) = (nb[x], nb[y]);
                let u = sys.pos[a] - sys.pos[center];
                let v = sys.pos[c] - sys.pos[center];
                let cos = (u.dot(v) / (u.norm() * v.norm())).clamp(-1.0, 1.0);
                bonded.angle(
                    a,
                    center,
                    c,
                    HarmonicAngle {
                        k: settings.k_angle,
                        theta0: cos.acos(),
                    },
                );
            }
        }
    }
    sys.exclusions = Some(Exclusions::from_bonded(
        &bonded,
        n,
        settings.fudge_lj,
        settings.fudge_qq,
    ));
    sys.bonded = Some(bonded);
    sys
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forces::ForceModel;
    use crate::minimiser::minimise;
    use tpt_chem_core::rng::Rng;
    use tpt_chem_core::vec3::Vec3;

    /// Methane with C-H = 1.09 A in a tetrahedral geometry (no bonds given).
    fn methane() -> Molecule {
        let mut m = Molecule::new("CH4");
        m.add_atom::<6>(Vec3::ZERO);
        let d = 1.09 / 3f64.sqrt();
        for (sx, sy, sz) in [
            (1.0, 1.0, 1.0),
            (1.0, -1.0, -1.0),
            (-1.0, 1.0, -1.0),
            (-1.0, -1.0, 1.0),
        ] {
            m.add_atom::<1>(Vec3::new(sx * d, sy * d, sz * d));
        }
        m
    }

    #[test]
    fn perceives_bonds_angles_and_exclusions() {
        let sys = build_system(&methane(), &TopologySettings::default());
        let b = sys.bonded.as_ref().unwrap();
        assert_eq!(b.bonds.len(), 4);
        assert_eq!(b.angles.len(), 6);
        // 4 C-H (1-2) + 6 H-C-H (1-3), nothing at 1-4.
        assert_eq!(sys.exclusions.as_ref().unwrap().len(), 10);
        // Tetrahedral reference angle.
        assert!((b.angles[0].1.theta0.to_degrees() - 109.4712).abs() < 1e-3);
        assert!((b.bonds[0].1.r0 - 1.09).abs() < 1e-12);
    }

    #[test]
    fn input_geometry_is_a_stationary_point() {
        let sys = build_system(&methane(), &TopologySettings::default());
        let (f, e) = ForceModel::AllPairs { cutoff: None }.evaluate(&sys);
        assert!(e.abs() < 1e-9, "{e}");
        assert!(f.iter().all(|v| v.norm() < 1e-8));
    }

    #[test]
    fn perturbed_molecule_relaxes_back() {
        let mut sys = build_system(&methane(), &TopologySettings::default());
        let reference: Vec<f64> = (1..5).map(|i| sys.pos[0].dist(sys.pos[i])).collect();
        let mut rng = Rng::new(3);
        for p in &mut sys.pos {
            *p += Vec3::new(
                rng.normal(0.0, 0.08),
                rng.normal(0.0, 0.08),
                rng.normal(0.0, 0.08),
            );
        }
        let model = ForceModel::AllPairs { cutoff: None };
        let rep = minimise(&mut sys, &model, 5000, 1e-3, 0.01);
        assert!(rep.final_energy < 1e-3, "{rep:?}");
        let _ = (rep.steps, rep.final_max_force);
        for (i, r0) in reference.iter().enumerate() {
            assert!(
                (sys.pos[0].dist(sys.pos[i + 1]) - r0).abs() < 5e-3,
                "{i}: {} vs {r0} ({rep:?})",
                sys.pos[0].dist(sys.pos[i + 1])
            );
        }
    }
}

#[cfg(test)]
mod dynamics_tests {
    use super::*;
    use crate::integrator::VelocityVerlet;
    use tpt_chem_core::rng::Rng;
    use tpt_chem_core::vec3::Vec3;

    /// Regression: the integrator once ignored `System::bonded`, so a
    /// flexible molecule flew apart. Methane must stay bound and conserve
    /// energy under plain velocity Verlet.
    #[test]
    fn methane_stays_bound_and_conserves_energy() {
        let mut m = Molecule::new("CH4");
        m.add_atom::<6>(Vec3::ZERO);
        let d = 1.09 / 3f64.sqrt();
        for (sx, sy, sz) in [
            (1.0, 1.0, 1.0),
            (1.0, -1.0, -1.0),
            (-1.0, 1.0, -1.0),
            (-1.0, -1.0, 1.0),
        ] {
            m.add_atom::<1>(Vec3::new(sx * d, sy * d, sz * d));
        }
        let mut sys = build_system(&m, &TopologySettings::default());
        sys.init_velocities(300.0, &mut Rng::new(8));
        sys.remove_com_motion();
        sys.update_forces();
        let vv = VelocityVerlet::new(0.25);
        let e0 = sys.total_energy();
        let mut max_dev = 0.0f64;
        let mut drift = 0.0f64;
        for _ in 0..4000 {
            vv.step(&mut sys, None);
            drift = drift.max((sys.total_energy() - e0).abs());
            for i in 1..5 {
                max_dev = max_dev.max((sys.pos[0].dist(sys.pos[i]) - 1.09).abs());
            }
        }
        assert!(max_dev < 0.4, "molecule dissociated: {max_dev}");
        assert!(drift < 0.02 * e0.abs().max(10.0), "drift {drift}, e0 {e0}");
    }
}
