//! `tpt-chem-kinetics` — reaction kinetics and population dynamics
//! (spec.txt §4, `tpt-chem-kinetics`).
//!
//! Deterministic mass-action kinetics (stiff ODE systems, solved with a
//! from-scratch Rosenbrock-type exponential integrator and classic
//! Runge-Kutta), the stochastic simulation algorithm (Gillespie SSA) for
//! low-copy-number regimes, Arrhenius and Eyring temperature-dependent
//! rate constants, and equilibrium / steady-state analysis.
//!
//! # Examples
//!
//! First-order decay A → B integrated deterministically:
//!
//! ```
//! use tpt_chem_kinetics::{network::ReactionNetwork, solver};
//!
//! // A → B with rate constant k = 1.0 (time units are whatever the
//! // rate constants are expressed in).
//! let net = ReactionNetwork::builder()
//!     .species(&["A", "B"])
//!     .reaction(&[("A", 1.0)], &[], 1.0)
//!     .build();
//! let y0 = [10.0, 0.0];
//!
//! let traj = solver::integrate(&net, &y0, 0.0, 5.0, 500);
//! let a_final = traj.last().unwrap()[0];
//! assert!((a_final - 10.0 * (-5.0_f64).exp()).abs() < 1e-6);
//! ```

pub mod arrhenius;
pub mod network;
pub mod solver;
