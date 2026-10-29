//! cc-pVDZ beyond neon: argon (Z = 18) has 4s, 3p and 1d shells, which is
//! 4 + 9 + 5 = 18 spherical basis functions.
//!
//! This test is slow in debug builds (minutes), so it is ignored by default.
//! Run it with `cargo test --release -p tpt-chem-quantum --test cc_pvdz_argon -- --ignored`.

use tpt_chem_core::molecule::Molecule;
use tpt_chem_core::vec3::Vec3;
use tpt_chem_quantum::basis::BasisSet;
use tpt_chem_quantum::hf::rhf_with_basis;

#[test]
#[ignore = "slow in debug builds; run with --release -- --ignored"]
fn argon_atom_basis_and_energy() {
    let mut ar = Molecule::new("Ar");
    ar.add_atom_raw(18, Vec3::ZERO);
    let res = rhf_with_basis(&ar, BasisSet::CcPvdz).unwrap();
    assert_eq!(res.n_basis, 18);
    // Closed-shell Ar. The Hartree–Fock limit is about −526.8 Eh and this
    // contracted basis sits a little above it, so check a range only.
    assert!(res.energy < -526.0 && res.energy > -527.5, "{}", res.energy);
}
