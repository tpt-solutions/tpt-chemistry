//! Serde round trip for MD state (requires the `serde` feature).

#![cfg(feature = "serde")]

use tpt_chem_core::forcefield::{HarmonicAngle, HarmonicBond, LennardJones};
use tpt_chem_core::vec3::Vec3;
use tpt_chem_md::bonded::Bonded;
use tpt_chem_md::box3::Box3;
use tpt_chem_md::exclusions::{Exclusions, PairScale};
use tpt_chem_md::system::System;

// Field-by-field setup keeps each parameter visible; `System` has no builder.
#[allow(clippy::field_reassign_with_default)]
fn sample() -> System {
    let mut sys = System::default();
    sys.pos = vec![
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
    ];
    sys.vel = vec![Vec3::default(); 3];
    sys.force = vec![Vec3::default(); 3];
    sys.lj = vec![
        LennardJones {
            sigma: 3.4,
            epsilon: 0.4
        };
        3
    ];
    sys.charge = vec![-0.8, 0.4, 0.4];
    sys.mass = vec![16.0, 1.0, 1.0];
    sys.names = vec!["O".into(), "H".into(), "H".into()];
    sys.box_ = Some(Box3::cubic(20.0));
    let mut bonded = Bonded::new();
    bonded.bond(0, 1, HarmonicBond { k: 500.0, r0: 0.96 });
    bonded.bond(0, 2, HarmonicBond { k: 500.0, r0: 0.96 });
    bonded.angles.push((
        [1, 0, 2],
        HarmonicAngle {
            k: 50.0,
            theta0: 1.82,
        },
    ));
    sys.bonded = Some(bonded);
    let mut ex = Exclusions::new();
    // Insert in a deliberately unsorted order.
    ex.set(
        2,
        1,
        PairScale {
            lj: 0.0,
            coulomb: 0.0,
        },
    );
    ex.set(
        0,
        2,
        PairScale {
            lj: 0.5,
            coulomb: 0.83,
        },
    );
    ex.set(
        0,
        1,
        PairScale {
            lj: 0.0,
            coulomb: 0.0,
        },
    );
    sys.exclusions = Some(ex);
    sys.potential_energy = -12.5;
    sys
}

#[test]
fn system_round_trips_with_topology_and_box() {
    let sys = sample();
    let json = serde_json::to_string(&sys).unwrap();
    let back: System = serde_json::from_str(&json).unwrap();
    assert_eq!(back.pos, sys.pos);
    assert_eq!(back.names, sys.names);
    assert_eq!(back.box_, sys.box_);
    assert_eq!(back.potential_energy, -12.5);
    assert_eq!(back.bonded.as_ref().unwrap().bonds.len(), 2);
    assert_eq!(back.bonded.as_ref().unwrap().angles.len(), 1);
    let ex = back.exclusions.unwrap();
    assert_eq!(ex.len(), 3);
    assert_eq!(
        ex.get(2, 1),
        Some(PairScale {
            lj: 0.0,
            coulomb: 0.0
        })
    );
    assert_eq!(
        ex.get(0, 2),
        Some(PairScale {
            lj: 0.5,
            coulomb: 0.83
        })
    );
}

#[test]
fn exclusion_serialization_is_deterministic() {
    // HashMap iteration order is random per map instance; the sorted repr
    // must make the output identical across instances.
    let a = serde_json::to_string(&sample()).unwrap();
    let b = serde_json::to_string(&sample()).unwrap();
    assert_eq!(a, b);
}
