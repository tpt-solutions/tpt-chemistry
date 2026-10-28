//! Serde round trip for quantum results (requires the `serde` feature).

#![cfg(feature = "serde")]

use tpt_chem_core::molecule::Molecule;
use tpt_chem_core::vec3::Vec3;
use tpt_chem_quantum::hf::{rhf, HfResult};

#[test]
fn hf_result_round_trips_through_json() {
    let mut h2 = Molecule::new("H2");
    h2.add_atom_raw(1, Vec3::new(0.0, 0.0, 0.0));
    h2.add_atom_raw(1, Vec3::new(0.0, 0.0, 0.7414));
    let res = rhf(&h2).unwrap();

    let json = serde_json::to_string(&res).unwrap();
    let back: HfResult = serde_json::from_str(&json).unwrap();
    assert_eq!(back.energy, res.energy);
    assert_eq!(back.n_basis, res.n_basis);
    assert_eq!(back.density, res.density);
    assert_eq!(back.orbital_energies, res.orbital_energies);
    // H2/STO-3G reference energy, Hartree.
    assert!((back.energy + 1.1168).abs() < 1e-3);
}
