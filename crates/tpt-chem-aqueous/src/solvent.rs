//! Temperature-dependent properties of liquid water at saturation pressure
//! and the Debye–Hückel parameters derived from them.

use crate::constants::{EPS0, E_CHARGE, K_B, N_A};

/// Density of liquid water, kg·m⁻³, for `t` in kelvin (Kell 1975
/// correlation, valid 0–100 °C at 1 atm).
pub fn water_density(t: f64) -> f64 {
    let c = t - 273.15;
    let num = 999.839_52 + 16.945_176 * c - 7.987_040_1e-3 * c.powi(2) - 46.170_461e-6 * c.powi(3)
        + 105.563_02e-9 * c.powi(4)
        - 280.542_53e-12 * c.powi(5);
    num / (1.0 + 16.879_850e-3 * c)
}

/// Static relative permittivity of water (Malmberg & Maryott 1956, valid
/// 0–100 °C), `t` in kelvin.
pub fn relative_permittivity(t: f64) -> f64 {
    let c = t - 273.15;
    87.740 - 0.400_08 * c + 9.398e-4 * c.powi(2) - 1.410e-6 * c.powi(3)
}

/// Bjerrum length `e²/(4π ε₀ ε_r k_B T)` of water, m.
pub fn bjerrum_length(t: f64) -> f64 {
    E_CHARGE * E_CHARGE / (4.0 * std::f64::consts::PI * EPS0 * relative_permittivity(t) * K_B * t)
}

/// Debye–Hückel osmotic-coefficient slope `A_φ` (natural-log convention,
/// (kg/mol)^½). About 0.391 at 25 °C.
pub fn debye_huckel_a_phi(t: f64) -> f64 {
    let lb = bjerrum_length(t);
    (2.0 * std::f64::consts::PI * N_A * water_density(t)).sqrt() * lb.powf(1.5) / 3.0
}

/// Debye–Hückel activity-coefficient slope `A_γ = 3 A_φ` (natural-log
/// convention), so that `ln γ_i = −A_γ z_i² √I` in the limiting law.
pub fn debye_huckel_a_gamma(t: f64) -> f64 {
    3.0 * debye_huckel_a_phi(t)
}

/// Debye–Hückel `B`, Å⁻¹·(kg/mol)^½ (about 0.329 at 25 °C), such that the
/// inverse screening length is `κ = B √I` in Å⁻¹.
pub fn debye_huckel_b(t: f64) -> f64 {
    let kappa2_per_i = 8.0 * std::f64::consts::PI * bjerrum_length(t) * N_A * water_density(t);
    kappa2_per_i.sqrt() * 1e-10
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn water_at_25c() {
        assert!((water_density(298.15) - 997.05).abs() < 0.05);
        assert!((relative_permittivity(298.15) - 78.4).abs() < 0.1);
    }

    #[test]
    fn debye_huckel_constants_at_25c() {
        assert!((debye_huckel_a_phi(298.15) - 0.3915).abs() < 0.002);
        assert!((debye_huckel_a_gamma(298.15) / std::f64::consts::LN_10 - 0.5091).abs() < 0.003);
        assert!((debye_huckel_b(298.15) - 0.3291).abs() < 0.002);
    }

    #[test]
    fn a_phi_increases_with_temperature() {
        assert!(debye_huckel_a_phi(343.15) > debye_huckel_a_phi(298.15));
    }
}
