//! `tpt-chem-io` — interoperability with standard chemical file formats
//! (spec.txt §4, `tpt-chem-io`).
//!
//! Parsers and writers for XYZ, PDB, MOL2, and CIF, plus streaming
//! XYZ trajectory I/O for large MD outputs. All parsers hand-rolled — no
//! regex, no external dependencies.
//!
//! # Examples
//!
//! ```
//! use tpt_chem_io::{parse_xyz, write_xyz};
//!
//! let source = "3\nwater (comment line)\nO  0.000  0.000  0.117\nH  0.000  0.757 -0.469\nH  0.000 -0.757 -0.469\n";
//! let mol = parse_xyz(source).unwrap();
//! assert_eq!(mol.formula(), "H2O");
//!
//! let text = write_xyz(&mol, "round trip");
//! let again = parse_xyz(&text).unwrap();
//! assert_eq!(again.atom(0usize.into()).pos, mol.atom(0usize.into()).pos);
//! ```

pub mod cif;
pub mod error;
pub mod gro;
pub mod mol2;
pub mod pdb;
pub mod trajectory;
pub mod xyz;

pub use error::{IoError, Result};
pub use trajectory::{XyzTrajectoryReader, XyzTrajectoryWriter};
pub use xyz::{parse_xyz, parse_xyz_all, write_xyz};
