//! Dissolved species and bulk solution descriptors.

/// A dissolved species with a net charge and a molality.
#[derive(Debug, Clone, PartialEq)]
pub struct Species {
    /// Name used to look up model parameters (e.g. `"Na+"`, `"Cl-"`).
    pub name: String,
    /// Net charge in elementary charges (0 for neutral species).
    pub charge: i32,
    /// Molality, mol per kg of water.
    pub molality: f64,
}

impl Species {
    /// Creates a species.
    pub fn new(name: &str, charge: i32, molality: f64) -> Self {
        Self {
            name: name.to_string(),
            charge,
            molality,
        }
    }
}

/// Molal ionic strength `I = ½ Σ mᵢ zᵢ²`, mol/kg.
pub fn ionic_strength(species: &[Species]) -> f64 {
    0.5 * species
        .iter()
        .map(|s| s.molality * f64::from(s.charge * s.charge))
        .sum::<f64>()
}

/// Net charge concentration `Σ mᵢ zᵢ`, equivalents per kg (zero for an
/// electroneutral solution).
pub fn charge_imbalance(species: &[Species]) -> f64 {
    species
        .iter()
        .map(|s| s.molality * f64::from(s.charge))
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ionic_strength_of_cacl2() {
        // 0.1 m CaCl2: I = ½(0.1·4 + 0.2·1) = 0.3
        let sol = [Species::new("Ca2+", 2, 0.1), Species::new("Cl-", -1, 0.2)];
        assert!((ionic_strength(&sol) - 0.3).abs() < 1e-15);
        assert!(charge_imbalance(&sol).abs() < 1e-15);
    }
}
