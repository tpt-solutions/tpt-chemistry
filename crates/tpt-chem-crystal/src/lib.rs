//! `tpt-chem-crystal` — solid-state physics and crystallography
//! (spec.txt §4, `tpt-chem-crystal`).
//!
//! Bravais lattices, Miller indices, Wigner-Seitz cell construction,
//! space-group symmetry operations with Wyckoff orbit generation,
//! reciprocal-lattice k-point grids, and powder XRD simulation from
//! structure factors.
//!
//! Atomic positions are fractional in the direct lattice; the lattice
//! vectors are row vectors of the conventional cell.

pub mod bravais;
pub mod kpoints;
pub mod miller;
pub mod symmetry;
pub mod wigner;
pub mod xrd;

pub use bravais::{Bravais, Lattice};
pub use kpoints::{KGrid, KPoint};
pub use miller::Miller;
pub use xrd::{Reflection, XrdSimulator};
