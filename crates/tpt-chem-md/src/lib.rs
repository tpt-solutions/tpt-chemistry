//! `tpt-chem-md` — the classical molecular dynamics engine.
//!
//! Symplectic Velocity Verlet integration, from-scratch cell lists and
//! Verlet neighbor lists for O(N) force evaluation, Lennard-Jones and
//! Coulomb interactions, Ewald summation, and Particle Mesh Ewald
//! long-range electrostatics on an in-house radix-2 FFT
//! (spec.txt §4, `tpt-chem-md`).
//!
//! # Units
//!
//! MD unit system throughout: positions Å, time fs, mass amu, energies
//! kJ·mol⁻¹, forces kJ·mol⁻¹·Å⁻¹.
//!
//! # Examples
//!
//! NVE dynamics of a small Lennard-Jones cluster:
//!
//! ```
//! use tpt_chem_core::forcefield::LennardJones;
//! use tpt_chem_core::rng::Rng;
//! use tpt_chem_md::{integrator::VelocityVerlet, system::System};
//!
//! let mut sys = System::new();
//! let lj = LennardJones { sigma: 3.4, epsilon: 0.997 };
//! for i in 0..13 {
//!     sys.add_atom(lj, 0.0, 39.95, [i as f64 * 4.0, 0.0, 0.0].into());
//! }
//! let mut rng = Rng::new(42);
//! sys.init_velocities(94.4, &mut rng); // ~94 K: Ar thermal energy
//! let vv = VelocityVerlet::new(1.0);   // 1 fs steps
//!
//! let _ = vv.run(&mut sys, 100, None);
//! let e1 = vv.run(&mut sys, 100, None);
//! // Symplectic integration: total energy drift stays bounded.
//! let e_pot = sys.potential_energy;
//! let e_tot = e_pot + sys.kinetic_energy();
//! assert!((e_tot - e1).abs() / e1.abs() < 1e-2);
//! ```

pub mod bonded;
pub mod box3;
pub mod bspline;
pub mod cell;
pub mod constraints;
pub mod ensemble;
pub mod ewald;
pub mod exclusions;
pub mod fft;
pub mod forces;
pub mod integrator;
pub mod minimiser;
pub mod neighbors;
pub mod observables;
pub mod pme;
pub mod system;
pub mod thermostat;
pub mod topology;
pub mod water;

pub use thermostat::{
    berendsen_barostat, berendsen_thermostat, langevin_step, nose_hoover_step, Leapfrog,
};

/// Conversion factor: kinetic energy per particle in
/// amu·Å²·fs⁻² to kJ·mol⁻¹ (`= amu·(10⁵ m·s⁻¹)²·N_A/1000`).
pub const KJ_MOL_PER_MDU: f64 = 1.0e4;

/// Conversion factor: acceleration in Å·fs⁻² from
/// force/mass in kJ·mol⁻¹·Å⁻¹·amu⁻¹ (`= 10³·10¹⁰/(N_A·amu·10²⁰)`).
pub const ACC_PER_FORCE_OVER_MASS: f64 = 1.0e-4;
