//! Serde round trip for crystal types (requires the `serde` feature).

#![cfg(feature = "serde")]

use tpt_chem_crystal::bravais::Lattice;

#[test]
fn lattice_round_trips() {
    let lat = Lattice::cubic(4.0);
    let back: Lattice = serde_json::from_str(&serde_json::to_string(&lat).unwrap()).unwrap();
    assert_eq!(back, lat);
}
