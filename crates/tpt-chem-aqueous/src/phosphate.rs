//! Constant set for the H₃PO₄ / Ca²⁺ / calcium-phosphate system.
//!
//! # Provenance
//!
//! The defaults are *literature-typical* values assembled for illustration
//! (pKa's of phosphoric acid, Ca–phosphate ion-pair constants and solubility
//! products of brushite, amorphous calcium phosphate, octacalcium phosphate
//! and hydroxyapatite). Reaction enthalpies of the solid phases and ion pairs
//! are rough estimates. They are **not** a vetted database: replace them with
//! measured values (all fields are public) for quantitative predictions.

use crate::equilibrium::Equilibrium;
use crate::solubility::{Component, SolidPhase};

/// Equilibrium constants for the Ca–phosphate–water system.
#[derive(Debug, Clone, PartialEq)]
pub struct PhosphateConstants {
    /// `H₃PO₄ ⇌ H⁺ + H₂PO₄⁻` (dissociation).
    pub ka1: Equilibrium,
    /// `H₂PO₄⁻ ⇌ H⁺ + HPO₄²⁻`.
    pub ka2: Equilibrium,
    /// `HPO₄²⁻ ⇌ H⁺ + PO₄³⁻`.
    pub ka3: Equilibrium,
    /// `H₂O ⇌ H⁺ + OH⁻`.
    pub kw: Equilibrium,
    /// `Ca²⁺ + H₂PO₄⁻ ⇌ CaH₂PO₄⁺` (association).
    pub ca_h2po4: Equilibrium,
    /// `Ca²⁺ + HPO₄²⁻ ⇌ CaHPO₄(aq)`.
    pub ca_hpo4: Equilibrium,
    /// `Ca²⁺ + PO₄³⁻ ⇌ CaPO₄⁻`.
    pub ca_po4: Equilibrium,
}

impl Default for PhosphateConstants {
    fn default() -> Self {
        Self {
            ka1: Equilibrium::new(-2.15, -8.0e3),
            ka2: Equilibrium::new(-7.20, 3.6e3),
            ka3: Equilibrium::new(-12.35, 16.0e3),
            kw: Equilibrium::new(-13.995, 55.8e3),
            ca_h2po4: Equilibrium::new(1.4, 6.0e3),
            ca_hpo4: Equilibrium::new(2.74, 12.0e3),
            ca_po4: Equilibrium::new(6.5, 15.0e3),
        }
    }
}

/// Brushite, `CaHPO₄·2H₂O ⇌ Ca²⁺ + HPO₄²⁻`.
pub fn brushite() -> SolidPhase {
    SolidPhase::new(
        "brushite",
        &[(Component::Ca, 1.0), (Component::Hpo4, 1.0)],
        Equilibrium::new(-6.59, -6.0e3),
    )
}

/// Amorphous calcium phosphate, `Ca₃(PO₄)₂·nH₂O ⇌ 3 Ca²⁺ + 2 PO₄³⁻`.
pub fn amorphous_calcium_phosphate() -> SolidPhase {
    SolidPhase::new(
        "amorphous calcium phosphate",
        &[(Component::Ca, 3.0), (Component::Po4, 2.0)],
        Equilibrium::new(-25.5, -10.0e3),
    )
}

/// Octacalcium phosphate, `Ca₈H₂(PO₄)₆·5H₂O ⇌ 8 Ca²⁺ + 2 H⁺ + 6 PO₄³⁻`.
pub fn octacalcium_phosphate() -> SolidPhase {
    SolidPhase::new(
        "octacalcium phosphate",
        &[
            (Component::Ca, 8.0),
            (Component::H, 2.0),
            (Component::Po4, 6.0),
        ],
        Equilibrium::new(-96.6, -20.0e3),
    )
}

/// Hydroxyapatite written per `Ca₅(PO₄)₃OH` (half the usual `Ca₁₀(PO₄)₆(OH)₂`
/// formula unit, so `log Ksp ≈ −58.4`): `Ca₅(PO₄)₃OH ⇌ 5 Ca²⁺ + 3 PO₄³⁻ + OH⁻`.
pub fn hydroxyapatite() -> SolidPhase {
    SolidPhase::new(
        "hydroxyapatite",
        &[
            (Component::Ca, 5.0),
            (Component::Po4, 3.0),
            (Component::Oh, 1.0),
        ],
        Equilibrium::new(-58.4, -20.0e3),
    )
}

/// Default candidate phases for precipitation onset: the kinetically
/// accessible precursors brushite, amorphous calcium phosphate and
/// octacalcium phosphate (Ostwald step rule). Hydroxyapatite is the
/// thermodynamic sink and is supersaturated across almost the whole
/// near-neutral range, so including it makes the "onset" meaningless as a
/// chalkiness predictor; use [`all_phases`] to include it deliberately.
pub fn default_phases() -> Vec<SolidPhase> {
    vec![
        brushite(),
        amorphous_calcium_phosphate(),
        octacalcium_phosphate(),
    ]
}

/// [`default_phases`] plus hydroxyapatite.
pub fn all_phases() -> Vec<SolidPhase> {
    let mut v = default_phases();
    v.push(hydroxyapatite());
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pka_ordering_is_sensible() {
        let c = PhosphateConstants::default();
        assert!(c.ka1.log_k_298 > c.ka2.log_k_298 && c.ka2.log_k_298 > c.ka3.log_k_298);
    }

    #[test]
    fn default_phases_are_named_and_distinct() {
        let phases = default_phases();
        assert_eq!(phases.len(), 3);
        assert_eq!(all_phases().len(), 4);
        assert!(phases.iter().all(|p| p.ksp.log_k_298 < 0.0));
    }
}
