//! Debye–Hückel family: limiting law, extended (Debye–Hückel with ion-size
//! term) and Davies. Valid roughly up to `I ≈ 0.01 m`, `0.1 m` and `0.5 m`
//! respectively; they depend only on ionic charge and ionic strength.

use super::{check_state, ActivityModel};
use crate::constants::LN10;
use crate::solvent::{debye_huckel_a_gamma, debye_huckel_b};
use crate::species::{ionic_strength, Species};
use crate::Result;

/// Debye–Hückel limiting law, `ln γ = −A_γ z² √I`.
#[derive(Debug, Clone, Copy, Default)]
pub struct LimitingLaw;

impl ActivityModel for LimitingLaw {
    fn ln_gamma(&self, species: &[Species], t: f64) -> Result<Vec<f64>> {
        check_state(species, t)?;
        let sqrt_i = ionic_strength(species).sqrt();
        let a = debye_huckel_a_gamma(t);
        Ok(species
            .iter()
            .map(|s| -a * f64::from(s.charge * s.charge) * sqrt_i)
            .collect())
    }
}

/// Extended Debye–Hückel, `ln γ = −A_γ z² √I / (1 + B å √I)`.
#[derive(Debug, Clone, Copy)]
pub struct ExtendedDebyeHuckel {
    /// Ion-size parameter å, Å (common default 4.5).
    pub ion_size: f64,
}

impl Default for ExtendedDebyeHuckel {
    fn default() -> Self {
        Self { ion_size: 4.5 }
    }
}

impl ActivityModel for ExtendedDebyeHuckel {
    fn ln_gamma(&self, species: &[Species], t: f64) -> Result<Vec<f64>> {
        check_state(species, t)?;
        let sqrt_i = ionic_strength(species).sqrt();
        let a = debye_huckel_a_gamma(t);
        let denom = 1.0 + debye_huckel_b(t) * self.ion_size * sqrt_i;
        Ok(species
            .iter()
            .map(|s| -a * f64::from(s.charge * s.charge) * sqrt_i / denom)
            .collect())
    }
}

/// Davies equation,
/// `log₁₀ γ = −A z² (√I/(1+√I) − c I)` with `A = A_γ/ln 10` and `c ≈ 0.3`.
#[derive(Debug, Clone, Copy)]
pub struct Davies {
    /// Linear-term coefficient (0.3 is standard, 0.2 is also used).
    pub c: f64,
}

impl Default for Davies {
    fn default() -> Self {
        Self { c: 0.3 }
    }
}

impl ActivityModel for Davies {
    fn ln_gamma(&self, species: &[Species], t: f64) -> Result<Vec<f64>> {
        check_state(species, t)?;
        let i = ionic_strength(species);
        let sqrt_i = i.sqrt();
        let a_ln = debye_huckel_a_gamma(t);
        let f = sqrt_i / (1.0 + sqrt_i) - self.c * i;
        Ok(species
            .iter()
            .map(|s| -a_ln * f64::from(s.charge * s.charge) * f)
            .collect())
    }
}

/// `log₁₀` helper: converts a natural-log coefficient to base 10.
pub fn log10_from_ln(ln_gamma: f64) -> f64 {
    ln_gamma / LN10
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nacl(m: f64) -> Vec<Species> {
        vec![Species::new("Na+", 1, m), Species::new("Cl-", -1, m)]
    }

    #[test]
    fn gamma_tends_to_one_at_infinite_dilution() {
        for model in [
            &LimitingLaw as &dyn ActivityModel,
            &ExtendedDebyeHuckel::default(),
            &Davies::default(),
        ] {
            let g = model.ln_gamma(&nacl(1e-10), 298.15).unwrap();
            assert!(g[0].abs() < 1e-4);
        }
    }

    #[test]
    fn limiting_law_slope() {
        // log10 γ± = −0.509 √I for 1:1 electrolytes at 25 °C.
        let g = LimitingLaw.ln_gamma(&nacl(0.001), 298.15).unwrap();
        let log_g = log10_from_ln(g[0]);
        assert!((log_g + 0.509 * 0.001_f64.sqrt()).abs() < 2e-4);
    }

    #[test]
    fn davies_for_tenth_molal_1_1_electrolyte() {
        // Davies: log10 γ = −0.509 (0.3162/1.3162 − 0.03) = −0.1079 → γ ≈ 0.780.
        let g = Davies::default().ln_gamma(&nacl(0.1), 298.15).unwrap();
        assert!((g[0].exp() - 0.780).abs() < 0.004);
    }

    #[test]
    fn divalent_ions_scale_with_charge_squared() {
        let sol = [Species::new("Ca2+", 2, 0.01), Species::new("Cl-", -1, 0.02)];
        let g = Davies::default().ln_gamma(&sol, 298.15).unwrap();
        assert!((g[0] / g[1] - 4.0).abs() < 1e-12);
    }

    #[test]
    fn rejects_negative_molality() {
        let sol = [Species::new("Na+", 1, -1.0)];
        assert!(Davies::default().ln_gamma(&sol, 298.15).is_err());
    }
}
