//! Solid phases, ion-activity products and saturation indices.

use crate::equilibrium::Equilibrium;

/// A dissolved component that can appear in a solubility product.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Component {
    /// Ca²⁺
    Ca,
    /// H⁺
    H,
    /// OH⁻
    Oh,
    /// H₂PO₄⁻
    H2po4,
    /// HPO₄²⁻
    Hpo4,
    /// PO₄³⁻
    Po4,
}

/// Activities (not logs) of the components, in the same order as
/// [`Component`] variants are declared.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ComponentLogActivities {
    /// `log₁₀ a(Ca²⁺)`
    pub ca: f64,
    /// `log₁₀ a(H⁺)`
    pub h: f64,
    /// `log₁₀ a(OH⁻)`
    pub oh: f64,
    /// `log₁₀ a(H₂PO₄⁻)`
    pub h2po4: f64,
    /// `log₁₀ a(HPO₄²⁻)`
    pub hpo4: f64,
    /// `log₁₀ a(PO₄³⁻)`
    pub po4: f64,
}

impl ComponentLogActivities {
    /// Looks up the value for a component.
    pub fn get(&self, c: Component) -> f64 {
        match c {
            Component::Ca => self.ca,
            Component::H => self.h,
            Component::Oh => self.oh,
            Component::H2po4 => self.h2po4,
            Component::Hpo4 => self.hpo4,
            Component::Po4 => self.po4,
        }
    }
}

/// A sparingly soluble solid written as a dissolution reaction
/// `solid ⇌ Σ νᵢ componentᵢ`.
#[derive(Debug, Clone, PartialEq)]
pub struct SolidPhase {
    /// Phase name.
    pub name: String,
    /// Stoichiometric coefficients of the dissolved components.
    pub stoichiometry: Vec<(Component, f64)>,
    /// Solubility product: `log₁₀ Ksp` at 25 °C and dissolution enthalpy.
    pub ksp: Equilibrium,
}

impl SolidPhase {
    /// Creates a phase.
    pub fn new(name: &str, stoichiometry: &[(Component, f64)], ksp: Equilibrium) -> Self {
        Self {
            name: name.to_string(),
            stoichiometry: stoichiometry.to_vec(),
            ksp,
        }
    }

    /// `log₁₀` of the ion-activity product.
    pub fn log_iap(&self, a: &ComponentLogActivities) -> f64 {
        self.stoichiometry
            .iter()
            .map(|&(c, nu)| nu * a.get(c))
            .sum()
    }

    /// Saturation index `SI = log₁₀(IAP/Ksp(T))`; positive means
    /// supersaturated.
    pub fn saturation_index(&self, a: &ComponentLogActivities, t: f64) -> f64 {
        self.log_iap(a) - self.ksp.log_k(t)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn saturation_index_is_zero_at_equilibrium() {
        // CaCO3-like 1:1 salt using Ca and HPO4: Ksp = 1e-7.
        let phase = SolidPhase::new(
            "test",
            &[(Component::Ca, 1.0), (Component::Hpo4, 1.0)],
            Equilibrium::new(-7.0, 0.0),
        );
        let a = ComponentLogActivities {
            ca: -3.5,
            h: -7.0,
            oh: -7.0,
            h2po4: -9.0,
            hpo4: -3.5,
            po4: -15.0,
        };
        assert!(phase.saturation_index(&a, 298.15).abs() < 1e-12);
        let more = ComponentLogActivities { ca: -3.0, ..a };
        assert!((phase.saturation_index(&more, 298.15) - 0.5).abs() < 1e-12);
    }

    #[test]
    fn retrograde_solid_becomes_less_soluble_when_hot() {
        let phase = SolidPhase::new("r", &[(Component::Ca, 1.0)], Equilibrium::new(-6.0, -20e3));
        assert!(phase.ksp.log_k(333.15) < phase.ksp.log_k(298.15));
    }
}
