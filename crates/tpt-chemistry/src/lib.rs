//! `tpt-chemistry` — feature-gated umbrella crate re-exporting the
//! `tpt-chem-*` computational-chemistry crates.
//!
//! Enable the crates you need via Cargo features:
//!
//! | feature | crate | contents |
//! |---|---|---|
//! | `core` | `tpt-chem-core` | molecular graphs, elements, force fields, unit-safe types (no_std) |
//! | `io` | `tpt-chem-io` | XYZ / PDB / MOL2 / CIF readers + writers, trajectory streaming |
//! | `md` | `tpt-chem-md` | integrators, thermostats, cell/neighbor lists, Ewald |
//! | `quantum` | `tpt-chem-quantum` | GTO basis sets, ERIs, SCF / Hartree-Fock |
//! | `kinetics` | `tpt-chem-kinetics` | mass-action networks, SSA, Arrhenius/Eyring |
//! | `crystal` | `tpt-chem-crystal` | Bravais lattices, space groups, XRD |
//! | `verify` | `tpt-chem-verify` | proptest strategies, conservation checks, Kani proofs |
//!
//! ```
//! # #[cfg(feature = "core")]
//! # {
//! use tpt_chemistry::prelude::*;
//!
//! let mut water = Molecule::new("water");
//! let o = water.add_atom::<8>([0.0, 0.0, 0.117].into());
//! let h1 = water.add_atom::<1>([0.0, 0.757, -0.469].into());
//! let h2 = water.add_atom::<1>([0.0, -0.757, -0.469].into());
//! water.add_bond(o, h1, BondOrder::Single).unwrap();
//! water.add_bond(o, h2, BondOrder::Single).unwrap();
//! assert_eq!(water.formula(), "H2O");
//! # }
//! ```

#[cfg(feature = "core")]
pub use tpt_chem_core;
#[cfg(feature = "crystal")]
pub use tpt_chem_crystal;
#[cfg(feature = "io")]
pub use tpt_chem_io;
#[cfg(feature = "kinetics")]
pub use tpt_chem_kinetics;
#[cfg(feature = "md")]
pub use tpt_chem_md;
#[cfg(feature = "quantum")]
pub use tpt_chem_quantum;
#[cfg(feature = "verify")]
pub use tpt_chem_verify;

/// Short aliases for the constituent crates (`tpt_chemistry::md`, …).
#[cfg(feature = "core")]
pub use tpt_chem_core as core;
#[cfg(feature = "crystal")]
pub use tpt_chem_crystal as crystal;
#[cfg(feature = "io")]
pub use tpt_chem_io as io;
#[cfg(feature = "kinetics")]
pub use tpt_chem_kinetics as kinetics;
#[cfg(feature = "md")]
pub use tpt_chem_md as md;
#[cfg(feature = "quantum")]
pub use tpt_chem_quantum as quantum;
#[cfg(feature = "verify")]
pub use tpt_chem_verify as verify;

/// Re-exports of the most-used types (behind the `core` feature).
#[cfg(feature = "core")]
pub mod prelude {
    pub use tpt_chem_core::prelude::*;
}

#[cfg(test)]
mod tests {
    #[test]
    fn core_feature_reexports() {
        #[cfg(feature = "core")]
        {
            use crate::prelude::*;
            let m = Molecule::new("t");
            assert_eq!(m.len(), 0);
        }
    }
}
