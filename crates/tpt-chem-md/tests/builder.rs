//! `System::builder` validation and initialization.

use tpt_chem_core::forcefield::{HarmonicBond, LennardJones};
use tpt_chem_core::vec3::Vec3;
use tpt_chem_md::bonded::Bonded;
use tpt_chem_md::box3::Box3;
use tpt_chem_md::exclusions::{Exclusions, PairScale};
use tpt_chem_md::prelude::*;
use tpt_chem_md::system::SystemError;

const LJ: LennardJones = LennardJones {
    sigma: 3.4,
    epsilon: 0.4,
};

fn three_atoms() -> tpt_chem_md::system::SystemBuilder {
    System::builder()
        .atom("O", Vec3::new(0.0, 0.0, 0.0), 16.0, -0.8, LJ)
        .atom("H", Vec3::new(0.96, 0.0, 0.0), 1.0, 0.4, LJ)
        .atom("H", Vec3::new(-0.24, 0.93, 0.0), 1.0, 0.4, LJ)
}

#[test]
fn atoms_keep_order_names_and_parameters() {
    let sys = three_atoms().build().unwrap();
    assert_eq!(sys.len(), 3);
    assert_eq!(sys.names, ["O", "H", "H"]);
    assert_eq!(sys.mass[0], 16.0);
    assert_eq!(sys.charge[0], -0.8);
    assert_eq!(sys.pos[1], Vec3::new(0.96, 0.0, 0.0));
}

#[test]
fn topology_index_out_of_range_is_rejected() {
    let mut bonded = Bonded::new();
    bonded.bond(0, 5, HarmonicBond { k: 500.0, r0: 0.96 });
    let err = three_atoms().bonded(bonded).build().unwrap_err();
    assert_eq!(
        err,
        SystemError::IndexOutOfRange {
            index: 5,
            n_atoms: 3
        }
    );
}

#[test]
fn exclusion_index_out_of_range_is_rejected() {
    let mut ex = Exclusions::new();
    ex.set(
        0,
        3,
        PairScale {
            lj: 0.0,
            coulomb: 0.0,
        },
    );
    assert!(three_atoms().exclusions(ex).build().is_err());
}

#[test]
fn valid_topology_is_accepted() {
    let mut bonded = Bonded::new();
    bonded.bond(0, 1, HarmonicBond { k: 500.0, r0: 0.96 });
    bonded.bond(0, 2, HarmonicBond { k: 500.0, r0: 0.96 });
    let sys = three_atoms()
        .bonded(bonded)
        .periodic(Box3::cubic(20.0))
        .build()
        .unwrap();
    assert_eq!(sys.bonded.as_ref().unwrap().bonds.len(), 2);
    assert!(sys.box_.is_some());
}

#[test]
fn thermal_velocities_are_seeded_and_drift_free() {
    let a = three_atoms().thermal_velocities(300.0, 42).build().unwrap();
    let b = three_atoms().thermal_velocities(300.0, 42).build().unwrap();
    assert_eq!(a.vel, b.vel, "same seed, same velocities");
    assert!(a.vel.iter().any(|v| v.norm_sq() > 0.0));
    // Total momentum is removed.
    let p = a
        .vel
        .iter()
        .zip(&a.mass)
        .fold(Vec3::ZERO, |acc, (v, m)| acc + *v * *m);
    assert!(p.norm_sq() < 1e-20, "momentum {p:?}");
}
