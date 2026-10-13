//! Multi-site Ca²⁺-binding ligands (phosphopeptides and similar chelators).
//!
//! Each [`BindingSite`] type is treated as an independent, non-cooperative
//! site that can be protonated (`HS`), deprotonated (`S`) or bound to
//! calcium (`CaS`):
//!
//! ```text
//! HS ⇌ H⁺ + S      (K_a)        S + Ca²⁺ ⇌ CaS     (K_Ca)
//! ```
//!
//! so the fraction of sites carrying calcium is
//! `K_Ca a(Ca) / (1 + a(H)/K_a + K_Ca a(Ca))`, which captures the competition
//! between protons and calcium (pH dependence) and, through van't Hoff
//! scaling of both constants, the temperature dependence. Site species are
//! assigned activity coefficient 1, and only one Ca²⁺ binds per site.
//! Calcium bridging between sites and electrostatic cooperativity are not
//! modelled.

use crate::equilibrium::Equilibrium;

/// One class of identical, independent binding sites on a ligand.
#[derive(Debug, Clone, PartialEq)]
pub struct BindingSite {
    /// Label (e.g. `"phosphoserine"`).
    pub name: String,
    /// Number of such sites per ligand molecule.
    pub count: f64,
    /// Charge of the protonated form `HS`.
    pub protonated_charge: i32,
    /// Proton dissociation `HS ⇌ H⁺ + S` (use `log K = −pKa`).
    pub proton: Equilibrium,
    /// Calcium association `S + Ca²⁺ ⇌ CaS`.
    pub calcium: Equilibrium,
}

/// A chelating ligand made of one or more site classes.
#[derive(Debug, Clone, PartialEq)]
pub struct Ligand {
    /// Ligand name.
    pub name: String,
    /// Site classes.
    pub sites: Vec<BindingSite>,
}

/// Equilibrium populations of one site class, molality of sites.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SitePopulation {
    /// Protonated sites `HS`.
    pub protonated: f64,
    /// Deprotonated, calcium-free sites `S`.
    pub free: f64,
    /// Calcium-bound sites `CaS`.
    pub bound: f64,
}

impl BindingSite {
    /// Populations of this site class for `sites` molal sites at the given
    /// calcium activity, proton activity and temperature.
    pub fn population(&self, sites: f64, a_ca: f64, a_h: f64, t: f64) -> SitePopulation {
        let ka = self.proton.k(t);
        let kca = self.calcium.k(t);
        let free = sites / (1.0 + a_h / ka + kca * a_ca);
        SitePopulation {
            protonated: free * a_h / ka,
            free,
            bound: free * kca * a_ca,
        }
    }
}

impl Ligand {
    /// Creates a ligand from its sites.
    pub fn new(name: &str, sites: Vec<BindingSite>) -> Self {
        Self {
            name: name.to_string(),
            sites,
        }
    }

    /// Molality of Ca²⁺ bound to `ligand_molality` molal ligand.
    pub fn bound_calcium(&self, ligand_molality: f64, a_ca: f64, a_h: f64, t: f64) -> f64 {
        self.sites
            .iter()
            .map(|s| s.population(ligand_molality * s.count, a_ca, a_h, t).bound)
            .sum()
    }

    /// A casein-phosphopeptide-like ligand: five phosphoserines
    /// (second pKa ≈ 5.8, log K_Ca ≈ 2.0) and six carboxylates (pKa ≈ 4.3,
    /// log K_Ca ≈ 0.9). **Illustrative defaults only** — substitute measured
    /// constants for the peptide of interest.
    pub fn casein_phosphopeptide() -> Self {
        Self::new(
            "casein phosphopeptide",
            vec![
                BindingSite {
                    name: "phosphoserine".to_string(),
                    count: 5.0,
                    protonated_charge: -1,
                    proton: Equilibrium::new(-5.8, 2.0e3),
                    calcium: Equilibrium::new(2.0, 5.0e3),
                },
                BindingSite {
                    name: "carboxylate".to_string(),
                    count: 6.0,
                    protonated_charge: 0,
                    proton: Equilibrium::new(-4.3, 0.0),
                    calcium: Equilibrium::new(0.9, 3.0e3),
                },
            ],
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site() -> BindingSite {
        BindingSite {
            name: "s".into(),
            count: 1.0,
            protonated_charge: 0,
            proton: Equilibrium::new(-5.0, 0.0),
            calcium: Equilibrium::new(3.0, 0.0),
        }
    }

    #[test]
    fn populations_sum_to_total_sites() {
        let p = site().population(2e-3, 1e-4, 1e-6, 298.15);
        assert!((p.protonated + p.free + p.bound - 2e-3).abs() < 1e-15);
    }

    #[test]
    fn half_bound_when_k_times_activity_is_one_at_high_ph() {
        // pH 12 ≫ pKa: sites are deprotonated, so bound fraction = K a/(1+K a).
        let p = site().population(1.0, 1e-3, 1e-12, 298.15);
        assert!((p.bound / 1.0 - 0.5).abs() < 1e-6);
    }

    #[test]
    fn calcium_binding_increases_with_ph() {
        let l = Ligand::casein_phosphopeptide();
        let low = l.bound_calcium(1e-3, 1e-3, 10f64.powf(-5.5), 298.15);
        let high = l.bound_calcium(1e-3, 1e-3, 10f64.powf(-7.5), 298.15);
        assert!(high > low);
    }

    #[test]
    fn no_ligand_binds_no_calcium() {
        assert_eq!(
            Ligand::casein_phosphopeptide().bound_calcium(0.0, 1e-3, 1e-7, 298.15),
            0.0
        );
    }
}
