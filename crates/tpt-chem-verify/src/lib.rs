//! `tpt-chem-verify` — verification harnesses for `tpt-chemistry`
//! (spec.txt §4, `tpt-chem-verify`).
//!
//! Property-based testing strategies (`proptest`) for generating valid
//! molecular geometries and force-field parameters, plus reusable
//! conservation-law checks. The Kani bounded-model-checking harnesses live
//! behind `#[cfg(kani)]` and run on Linux via CI (Kani does not support
//! Windows hosts).
//!
//! # Examples
//!
//! Reuse the energy-conservation check on a custom system:
//!
//! ```
//! use tpt_chem_core::forcefield::LennardJones;
//! use tpt_chem_core::rng::Rng;
//! use tpt_chem_core::vec3::Vec3;
//! use tpt_chem_md::box3::Box3;
//! use tpt_chem_md::system::System;
//! use tpt_chem_verify::conservation::assert_energy_conserved;
//!
//! let mut sys = System::new();
//! let lj = LennardJones { sigma: 3.4, epsilon: 0.997 };
//! let mut rng = Rng::new(1);
//! for i in 0..8 {
//!     sys.add_atom(lj, 0.0, 39.95, Vec3::new(i as f64 * 4.5 + 2.0, 2.0, 2.0));
//! }
//! sys.set_box(Box3::cubic(40.0));
//! sys.init_velocities(94.0, &mut rng);
//! // 1 fs Velocity Verlet on a clean geometry conserves energy tightly.
//! assert_energy_conserved(&mut sys, 200, 1.0, 5e-3);
//! ```

pub mod conservation;
pub mod strategies;
