//! `tpt-chem-aqueous` — aqueous-solution thermodynamics: activity
//! coefficients, solubility products and chelator-controlled calcium
//! phosphate precipitation.
//!
//! * [`activity`] — [`ActivityModel`](activity::ActivityModel) implementations
//!   for multi-component electrolytes: limiting-law / extended Debye–Hückel,
//!   Davies, Pitzer virial equations and a single-solvent electrolyte NRTL.
//! * [`equilibrium`] — van't Hoff temperature dependence of equilibrium
//!   constants.
//! * [`solubility`] — solid phases, ion-activity products and saturation
//!   indices.
//! * [`phosphate`] — the H₃PO₄ / calcium-phosphate constant set and phases.
//! * [`chelator`] — multi-site Ca²⁺-binding ligands (phosphopeptides, …).
//! * [`precipitation`] — pH × temperature scans and onset-pH search that
//!   predict when a calcium phosphate phase becomes supersaturated
//!   ("chalkiness").
//!
//! # Units convention
//!
//! Temperatures in kelvin, concentrations as molalities (mol per kg of
//! water), enthalpies in J·mol⁻¹, pH on the activity scale
//! (`a(H⁺) = 10^(−pH)`), equilibrium constants as `log₁₀ K` at 25 °C.
//!
//! # Scope and caveats
//!
//! Everything here is *thermodynamic*: a saturation index above zero means
//! the solution is supersaturated with respect to a phase, not that the
//! phase has nucleated. Tabulated constants shipped as defaults are
//! literature-typical illustrative values and **must be replaced with
//! measured data** for quantitative work; every value is overridable.
//!
//! # Examples
//!
//! Mean activity coefficient of 1 m NaCl from the Pitzer model:
//!
//! ```
//! use tpt_chem_aqueous::activity::{ActivityModel, pitzer::{Pitzer, PitzerPair}};
//! use tpt_chem_aqueous::species::Species;
//!
//! let mut pitzer = Pitzer::new();
//! pitzer.add_pair(PitzerPair::new("Na+", "Cl-", 0.0765, 0.2664, 0.0, 0.00127));
//! let sol = [Species::new("Na+", 1, 1.0), Species::new("Cl-", -1, 1.0)];
//! let ln_g = pitzer.ln_gamma(&sol, 298.15).unwrap();
//! let gamma_pm = ((ln_g[0] + ln_g[1]) / 2.0).exp();
//! assert!((gamma_pm - 0.657).abs() < 0.005);
//! ```

pub mod activity;
pub mod chelator;
pub mod constants;
pub mod equilibrium;
pub mod phosphate;
pub mod precipitation;
pub mod solubility;
pub mod solvent;
mod solver;
pub mod species;

use std::fmt;

/// Errors produced by this crate.
#[derive(Debug, Clone, PartialEq)]
pub enum AqueousError {
    /// An input was outside its valid range.
    InvalidInput(&'static str),
    /// An iterative solver did not converge.
    NoConvergence(&'static str),
    /// A root-finding bracket did not contain a sign change.
    NotBracketed,
    /// A computation produced NaN or infinity.
    NonFinite(&'static str),
    /// The requested operation is not supported by this model.
    Unsupported(&'static str),
}

impl fmt::Display for AqueousError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidInput(m) => write!(f, "invalid input: {m}"),
            Self::NoConvergence(m) => write!(f, "no convergence: {m}"),
            Self::NotBracketed => write!(f, "root not bracketed"),
            Self::NonFinite(m) => write!(f, "non-finite value in {m}"),
            Self::Unsupported(m) => write!(f, "unsupported: {m}"),
        }
    }
}

impl std::error::Error for AqueousError {}

/// Result alias for this crate.
pub type Result<T> = std::result::Result<T, AqueousError>;
