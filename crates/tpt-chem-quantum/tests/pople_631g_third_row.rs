//! 6-31G for Na–Ar (third row), from Basis Set Exchange.
//!
//! No measured or literature total energy was available to check against, so
//! these tests cover what can be checked here: SCF convergence for an
//! argon atom, and a sanity window. The window comes from the requirement
//! that a variational energy lies above the Hartree–Fock limit (about
//! −526.82 Eₕ for Ar), together with the size of the basis-set error that is
//! typical for 6-31G on a third-row atom.

use tpt_chem_core::molecule::Molecule;
use tpt_chem_core::vec3::Vec3;
use tpt_chem_quantum::basis::BasisSet;
use tpt_chem_quantum::hf::rhf_with_basis;

#[test]
fn argon_atom_6_31g_converges_in_sanity_window() {
    let mut ar = Molecule::new("Ar");
    ar.add_atom::<18>(Vec3::ZERO);
    let res = rhf_with_basis(&ar, BasisSet::Pople631g).expect("SCF converges");
    assert_eq!(res.n_basis, 13);
    assert!(
        res.energy > -526.82 && res.energy < -526.70,
        "Ar 6-31G energy {:.8}",
        res.energy
    );
}
