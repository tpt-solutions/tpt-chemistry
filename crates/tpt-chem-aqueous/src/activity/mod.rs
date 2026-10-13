//! Activity-coefficient models for multi-component aqueous electrolytes.
//!
//! All models implement [`ActivityModel`] and return **molality-scale,
//! unsymmetric** (infinite-dilution reference) activity coefficients, so they
//! can be used interchangeably by the speciation and solubility code.
//! Neutral solutes are assigned `γ = 1` by every model in this crate.

pub mod debye_huckel;
pub mod enrtl;
pub mod pitzer;

use crate::species::Species;
use crate::{AqueousError, Result};

/// A model for the activity coefficients of dissolved species.
pub trait ActivityModel {
    /// Natural logarithm of the molality-scale activity coefficient of every
    /// species in `species` (same order) at temperature `temperature` (K).
    fn ln_gamma(&self, species: &[Species], temperature: f64) -> Result<Vec<f64>>;

    /// Osmotic coefficient `φ` of the solution (water activity follows from
    /// `ln a_w = −φ M_w Σ mᵢ`). Models without a consistent definition return
    /// [`AqueousError::Unsupported`].
    fn osmotic_coefficient(&self, _species: &[Species], _temperature: f64) -> Result<f64> {
        Err(AqueousError::Unsupported("osmotic coefficient"))
    }
}

pub(crate) fn check_state(species: &[Species], temperature: f64) -> Result<()> {
    if !(temperature.is_finite() && temperature > 200.0 && temperature < 700.0) {
        return Err(AqueousError::InvalidInput(
            "temperature must be in 200–700 K",
        ));
    }
    if species
        .iter()
        .any(|s| !s.molality.is_finite() || s.molality < 0.0)
    {
        return Err(AqueousError::InvalidInput(
            "molalities must be finite and non-negative",
        ));
    }
    Ok(())
}
