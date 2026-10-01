//! Reaction networks: species, reactions, and mass-action propensities.
//!
//! A [`ReactionNetwork`] holds `S` species and `R` elementary reactions.
//! Each reaction carries a stoichiometry vector (net change per species)
//! and a rate constant; the mass-action propensity of reaction `r` in
//! state `x` is `k_r × ∏_i x_i^{ν⁻_{r,i}}` where `ν⁻` counts reactant
//! stoichiometries.
//!
//! Rate constants may be temperature-independent constants or evaluated
//! at a temperature through [`crate::arrhenius::Arrhenius`] /
//! [`crate::arrhenius::Eyring`] laws (see [`ReactionInput`]).

use std::string::String;
use std::vec::Vec;

use crate::arrhenius::RateLaw;

/// How a reaction's rate constant is obtained.
#[derive(Clone, Debug)]
pub enum ReactionInput {
    /// Temperature-independent constant.
    Constant(f64),
    /// Temperature-dependent law evaluated at the system temperature.
    TemperatureDependent(RateLaw),
}

/// One elementary reaction.
#[derive(Clone, Debug)]
pub struct Reaction {
    /// Reactant stoichiometry: `(species index, multiplicity)` pairs.
    pub reactants: Vec<(usize, f64)>,
    /// Product stoichiometry.
    pub products: Vec<(usize, f64)>,
    /// Net stoichiometry change `ν_r` (products − reactants), length `S`.
    pub stoichiometry: Vec<f64>,
    /// Rate input (constant or temperature law).
    pub rate: ReactionInput,
}

impl Reaction {
    /// Mass-action propensity `k(x) × ∏ x^ν⁻` in state `x`.
    ///
    /// Returns 0 when any reactant concentration is non-positive.
    pub fn propensity(&self, x: &[f64]) -> f64 {
        let mut rate = match &self.rate {
            ReactionInput::Constant(k) => *k,
            ReactionInput::TemperatureDependent(law) => law.rate_at_300k(),
        };
        for &(idx, mult) in &self.reactants {
            let v = x[idx];
            if v <= 0.0 {
                return 0.0;
            }
            rate *= v.powf(mult);
        }
        rate
    }
}

/// Errors from network construction.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum NetworkError {
    /// A reaction references a species index that does not exist.
    BadSpeciesIndex,
    /// The builder was not finalized (`build` without `species`).
    Incomplete,
}

impl core::fmt::Display for NetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            NetworkError::BadSpeciesIndex => write!(f, "reaction references unknown species"),
            NetworkError::Incomplete => write!(f, "network builder incomplete"),
        }
    }
}

impl std::error::Error for NetworkError {}

/// A mass-action reaction network.
#[derive(Clone, Debug, Default)]
pub struct ReactionNetwork {
    /// Species names, in state-vector order.
    pub species: Vec<String>,
    /// Elementary reactions.
    pub reactions: Vec<Reaction>,
}

impl ReactionNetwork {
    /// An empty network.
    pub fn new() -> Self {
        Self::default()
    }

    /// Start a builder with the given species names (defines the
    /// state-vector order).
    pub fn builder() -> ReactionNetworkBuilder {
        ReactionNetworkBuilder::default()
    }

    /// Number of species.
    pub fn n_species(&self) -> usize {
        self.species.len()
    }

    /// Number of reactions.
    pub fn n_reactions(&self) -> usize {
        self.reactions.len()
    }

    /// All stoichiometry vectors as an `R × S` matrix (row-major).
    pub fn stoichiometry_matrix(&self) -> Vec<Vec<f64>> {
        self.reactions
            .iter()
            .map(|r| r.stoichiometry.clone())
            .collect()
    }

    /// Propensities of all reactions in state `x`.
    pub fn propensities(&self, x: &[f64]) -> Vec<f64> {
        self.reactions.iter().map(|r| r.propensity(x)).collect()
    }

    /// The mass-action derivative `dx/dt = Σ_r ν_r k_r(x)`.
    pub fn derivative(&self, x: &[f64]) -> Vec<f64> {
        let mut dx = vec![0.0; self.n_species()];
        for r in &self.reactions {
            let a = r.propensity(x);
            for (i, &nu) in r.stoichiometry.iter().enumerate() {
                dx[i] += nu * a;
            }
        }
        dx
    }
}

/// Builder for [`ReactionNetwork`].
#[derive(Clone, Debug, Default)]
pub struct ReactionNetworkBuilder {
    species: Vec<String>,
    reactions: Vec<Reaction>,
    index_of: Vec<(String, usize)>,
}

impl ReactionNetworkBuilder {
    /// Declare species by name (order defines the state vector).
    #[must_use]
    pub fn species(mut self, names: &[&str]) -> Self {
        for (i, n) in names.iter().enumerate() {
            self.species.push((*n).into());
            self.index_of.push(((*n).into(), i));
        }
        self
    }

    fn lookup(&self, name: &str) -> Option<usize> {
        self.index_of
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, i)| *i)
    }

    /// Add a reaction: reactant and product `(name, multiplicity)` pairs
    /// with a temperature-independent rate constant.
    #[must_use]
    pub fn reaction(self, reactants: &[(&str, f64)], products: &[(&str, f64)], rate: f64) -> Self {
        self.reaction_with(reactants, products, ReactionInput::Constant(rate))
    }

    /// Add a reaction with an explicit [`ReactionInput`] (e.g. an
    /// Arrhenius law).
    #[must_use]
    pub fn reaction_with(
        mut self,
        reactants: &[(&str, f64)],
        products: &[(&str, f64)],
        rate: ReactionInput,
    ) -> Self {
        let n = self.species.len();
        let mut stoich = vec![0.0; n];
        let mut r_list = Vec::new();
        for (name, mult) in reactants {
            if let Some(i) = self.lookup(name) {
                stoich[i] -= mult;
                r_list.push((i, *mult));
            }
        }
        let mut p_list = Vec::new();
        for (name, mult) in products {
            if let Some(i) = self.lookup(name) {
                stoich[i] += mult;
                p_list.push((i, *mult));
            }
        }
        self.reactions.push(Reaction {
            reactants: r_list,
            products: p_list,
            stoichiometry: stoich,
            rate,
        });
        self
    }

    /// Finalize the network.
    #[must_use]
    pub fn build(self) -> ReactionNetwork {
        ReactionNetwork {
            species: self.species,
            reactions: self.reactions,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_order_propensity_and_derivative() {
        let net = ReactionNetwork::builder()
            .species(&["A", "B"])
            .reaction(&[("A", 1.0)], &[("B", 1.0)], 2.0)
            .build();
        let prop = net.propensities(&[3.0, 1.0]);
        assert_eq!(prop, vec![6.0]);
        let dx = net.derivative(&[3.0, 1.0]);
        assert_eq!(dx[0], -6.0);
        assert_eq!(dx[1], 6.0);
    }

    #[test]
    fn bimolecular_propensity() {
        // 2H → H2 (rate k): propensity = k [H]^2.
        let net = ReactionNetwork::builder()
            .species(&["H", "H2"])
            .reaction(&[("H", 2.0)], &[("H2", 1.0)], 0.5)
            .build();
        assert_eq!(net.propensities(&[4.0, 1.0])[0], 0.5 * 16.0);
        let dx = net.derivative(&[4.0, 0.0]);
        // Stoichiometry for H is -2, so the rate acts twice.
        assert_eq!(dx[0], -16.0);
        assert_eq!(dx[1], 8.0);
    }

    #[test]
    fn reversible_equilibrium_constant() {
        // A ⇌ B with kf = 3, kr = 1: K = 3, equilibrium [B]/[A] = 3.
        let net = ReactionNetwork::builder()
            .species(&["A", "B"])
            .reaction(&[("A", 1.0)], &[("B", 1.0)], 3.0)
            .reaction(&[("B", 1.0)], &[("A", 1.0)], 1.0)
            .build();
        // Equilibrium: 3[A] = [B], total = 1 → [A] = 0.25.
        let mut x = [1.0, 0.0];
        for _ in 0..200 {
            let dx = net.derivative(&x);
            for (xi, dxi) in x.iter_mut().zip(&dx) {
                *xi += 0.05 * dxi;
            }
        }
        assert!((x[0] - 0.25).abs() < 0.01, "A = {}", x[0]);
        assert!((x[1] - 0.75).abs() < 0.01, "B = {}", x[1]);
    }

    #[test]
    fn zero_reactant_gives_zero_propensity() {
        let net = ReactionNetwork::builder()
            .species(&["A", "B"])
            .reaction(&[("A", 1.0)], &[("B", 1.0)], 1.0)
            .build();
        assert_eq!(net.propensities(&[0.0, 5.0])[0], 0.0);
    }
}
