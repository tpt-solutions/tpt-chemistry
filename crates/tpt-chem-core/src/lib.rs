//! `tpt-chem-core` — the foundation layer of `tpt-chemistry`.
//!
//! Molecular graphs, atom types, force-field parameterizations, physical
//! constants, and unit-safe coordinate types. Per spec.txt §4 this crate is
//! `no_std + alloc` compatible: disable the default `std` feature and it
//! builds for bare-metal targets (transcendental float math routes through
//! `libm`).
//!
//! # Units convention
//!
//! Every quantity crossing a public API is either explicitly unit-wrapped
//! (see [`units`]) or documented in one of two simulation unit systems:
//!
//! * **MD units** — length Å, time fs, mass amu, energy kJ·mol⁻¹.
//! * **Atomic units** — length Bohr, time ħ/Eₕ, mass mₑ, energy Eₕ
//!   (Hartree).
//!
//! # Examples
//!
//! Build a water molecule and measure its O–H distance:
//!
//! ```
//! use tpt_chem_core::molecule::{BondOrder, Molecule};
//!
//! let mut water = Molecule::new("water");
//! let o = water.add_atom::<8>([0.000, 0.000, 0.117].into());
//! let h1 = water.add_atom::<1>([0.000, 0.757, -0.469].into());
//! let h2 = water.add_atom::<1>([0.000, -0.757, -0.469].into());
//! water.add_bond(o, h1, BondOrder::Single).unwrap();
//! water.add_bond(o, h2, BondOrder::Single).unwrap();
//!
//! assert!((water.distance(o, h1) - 0.9577).abs() < 1e-3);
//! assert_eq!(water.formula(), "H2O");
//! ```

#![cfg_attr(not(feature = "std"), no_std)]
#![allow(clippy::needless_range_loop)]

extern crate alloc;

pub mod coords;
pub mod element;
pub mod forcefield;
pub mod molecule;
pub mod num;
pub mod potential;
pub mod rng;
pub mod special;
pub mod units;
pub mod vec3;

/// Convenient re-exports of the most-used core types.
pub mod prelude {
    pub use crate::coords::{Cartesian, Internal, Position};
    pub use crate::element::Element;
    pub use crate::forcefield::{
        Coulomb, HarmonicAngle, HarmonicBond, LennardJones, RyckaertBellemans,
    };
    pub use crate::molecule::{AtomId, BondOrder, Molecule};
    pub use crate::units::*;
    pub use crate::vec3::Vec3;
}
