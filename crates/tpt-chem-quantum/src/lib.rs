//! `tpt-chem-quantum` — electronic structure primitives
//! (spec.txt §4, `tpt-chem-quantum`).
//!
//! Gaussian-type-orbital basis sets (STO-3G for H–Ne), from-scratch
//! Obara–Saika overlap/kinetic and McMurchie–Davidson attraction/ERI
//! evaluation of one- and two-electron integrals, and a
//! restricted closed-shell Hartree–Fock self-consistent-field driver
//! solving the Roothaan–Hall equations with a from-scratch cyclic Jacobi
//! eigensolver.
//!
//! Everything runs in **atomic units** (Bohr, Hartree); molecular
//! geometries from `tpt-chem-core` (Ångström) are converted at the API
//! boundary.
//!
//! # Examples
//!
//! Restricted Hartree–Fock for H₂ at its experimental bond length:
//!
//! ```
//! use tpt_chem_core::molecule::Molecule;
//! use tpt_chem_core::units::BOHR_ANGSTROM;
//! use tpt_chem_core::vec3::Vec3;
//! use tpt_chem_quantum::hf::rhf_energy;
//!
//! let mut h2 = Molecule::new("H2");
//! let r = 1.4 * BOHR_ANGSTROM; // 1.4 Bohr in Angstrom
//! h2.add_atom::<1>(Vec3::new(0.0, 0.0, -r / 2.0));
//! h2.add_atom::<1>(Vec3::new(0.0, 0.0, r / 2.0));
//!
//! // HF/STO-3G literature value: −1.11675931 Eₕ (Schwartz & Schaad 1967).
//! let e = rhf_energy(&h2).unwrap();
//! assert!((e - (-1.11675931)).abs() < 5e-4);
//! ```

pub mod basis;
pub mod boys;
pub mod gaussian;
pub mod hf;
pub mod integrals;
pub mod linalg;
