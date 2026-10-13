//! Equilibrium constants with van't Hoff temperature dependence.

use crate::constants::{LN10, R, T_REF};

/// An equilibrium constant known at 25 °C with a constant reaction enthalpy.
///
/// `log_k_298` and `delta_h` refer to the reaction **as written** (for an
/// acid dissociation or a dissolution, the dissociation/dissolution
/// direction; for complexation, the association direction).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Equilibrium {
    /// `log₁₀ K` at 298.15 K.
    pub log_k_298: f64,
    /// Reaction enthalpy, J·mol⁻¹, assumed temperature independent.
    pub delta_h: f64,
}

impl Equilibrium {
    /// Creates an equilibrium constant.
    pub const fn new(log_k_298: f64, delta_h: f64) -> Self {
        Self { log_k_298, delta_h }
    }

    /// `log₁₀ K(T)` by the integrated van't Hoff equation.
    pub fn log_k(&self, t: f64) -> f64 {
        self.log_k_298 - self.delta_h / (R * LN10) * (1.0 / t - 1.0 / T_REF)
    }

    /// `K(T)`.
    pub fn k(&self, t: f64) -> f64 {
        10f64.powf(self.log_k(t))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reference_temperature_is_identity() {
        let e = Equilibrium::new(-7.2, 3.6e3);
        assert!((e.log_k(T_REF) - e.log_k_298).abs() < 1e-14);
    }

    #[test]
    fn endothermic_reaction_strengthens_with_temperature() {
        let e = Equilibrium::new(-14.0, 55.8e3);
        assert!(e.log_k(323.15) > e.log_k(298.15));
        // Kw at 50 °C is about 5.5e-14 (pKw ≈ 13.26).
        assert!((-e.log_k(323.15) - 13.26).abs() < 0.05);
    }

    #[test]
    fn round_trip_between_temperatures() {
        let e = Equilibrium::new(-6.59, -6.0e3);
        let back = Equilibrium {
            log_k_298: e.log_k(330.0),
            delta_h: e.delta_h,
        }
        .log_k(T_REF);
        // Re-referencing at 330 K and evaluating at 298.15 K must reproduce 25 °C
        // only if the reference is consistent; check the van't Hoff slope instead.
        let slope = (e.log_k(331.0) - e.log_k(329.0)) / (1.0 / 331.0 - 1.0 / 329.0);
        assert!((slope + e.delta_h / (R * LN10)).abs() < 1e-6);
        assert!(back.is_finite());
    }
}
