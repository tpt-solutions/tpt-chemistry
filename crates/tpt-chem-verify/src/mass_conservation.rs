//! Mass conservation for closed reaction networks (spec.txt §6, Phase 1).
//!
//! A closed network's atom inventory is invariant: the total atomic mass
//! `Σ m_i x_i` must be constant along deterministic trajectories and in
//! stochastic-sample averages.

use std::vec::Vec;

use tpt_chem_kinetics::network::ReactionNetwork;

/// Total weighted mass `Σ m_i x_i` of a state vector given per-species
/// molar masses (any positive masses work: conservation is a property of
/// the stoichiometry, not the values).
pub fn total_mass(masses: &[f64], x: &[f64]) -> f64 {
    masses.iter().zip(x).map(|(m, xi)| m * xi).sum()
}

/// Build a mass-balanced test network set with per-species molar masses:
/// returns `(network, masses)` pairs whose reactions are atom-balanced.
pub fn sample_networks() -> Vec<(ReactionNetwork, Vec<f64>)> {
    vec![
        // Isomerization A ⇌ B (identical molar mass: same atoms).
        (
            ReactionNetwork::builder()
                .species(&["A", "B"])
                .reaction(&[("A", 1.0)], &[("B", 1.0)], 1.3)
                .reaction(&[("B", 1.0)], &[("A", 1.0)], 0.7)
                .build(),
            vec![12.011, 12.011],
        ),
        // Dimerization 2O ⇌ O2 (2 × 16.0 = 32.0).
        (
            ReactionNetwork::builder()
                .species(&["O", "O2"])
                .reaction(&[("O", 2.0)], &[("O2", 1.0)], 0.4)
                .reaction(&[("O2", 1.0)], &[("O", 2.0)], 0.1)
                .build(),
            vec![16.0, 32.0],
        ),
        // Hydrogenation H2 + O → H2O (2 × 1.008 + 16.0 = 18.016).
        (
            ReactionNetwork::builder()
                .species(&["H2", "O", "H2O"])
                .reaction(&[("H2", 1.0), ("O", 1.0)], &[("H2O", 1.0)], 0.25)
                .reaction(&[("H2O", 1.0)], &[("H2", 1.0), ("O", 1.0)], 0.05)
                .build(),
            vec![2.016, 16.0, 18.016],
        ),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    use tpt_chem_core::rng::Rng;
    use tpt_chem_kinetics::solver;

    #[test]
    fn rk4_conserves_total_mass() {
        for (net, masses) in sample_networks() {
            let y0: Vec<f64> = net.species.iter().map(|_| 1.0).collect();
            let m0 = total_mass(&masses, &y0);
            let traj = solver::integrate(&net, &y0, 0.0, 2.0, 400);
            for state in &traj {
                let m = total_mass(&masses, state);
                assert!(
                    (m - m0).abs() < 1e-9,
                    "{}: mass {m} vs {m0}",
                    net.species[0]
                );
            }
        }
    }

    #[test]
    fn stiff_conserves_total_mass() {
        for (net, masses) in sample_networks() {
            let y0: Vec<f64> = net.species.iter().map(|_| 1.0).collect();
            let m0 = total_mass(&masses, &y0);
            let traj = solver::integrate_stiff(&net, &y0, 0.0, 1.0, 100, 16);
            for state in &traj {
                let m = total_mass(&masses, state);
                assert!((m - m0).abs() < 1e-6, "mass {m} vs {m0}");
            }
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(32))]

        /// spec.txt §6 Phase 1 exit criterion: mass is strictly conserved
        /// in closed reaction networks, over randomized initial states.
        #[test]
        fn rk4_mass_conservation_random_states(
            a in 0.1f64..5.0,
            b in 0.1f64..5.0,
            c in 0.1f64..5.0,
        ) {
            for (net, masses) in sample_networks() {
                let n = net.n_species();
                let y0: Vec<f64> = match n {
                    2 => vec![a, b],
                    _ => vec![a, b, c],
                };
                let m0 = total_mass(&masses, &y0);
                let traj = solver::integrate(&net, &y0, 0.0, 1.5, 300);
                let last = traj.last().unwrap();
                let m = total_mass(&masses, last);
                prop_assert!((m - m0).abs() < 1e-8, "mass {m} vs {m0}");
            }
        }
    }

    #[test]
    fn ssa_conserves_total_mass_per_event() {
        // Hydrogenation: H2 + O → H2O keeps total mass fixed per event.
        let (net, masses) = sample_networks()[2].clone();
        let mut rng = Rng::new(4242);
        let y0 = [2.0, 2.0, 0.0];
        let m0 = total_mass(&masses, &y0);
        for _ in 0..50 {
            let final_state = solver::ssa(&net, &y0, 3.0, &mut rng);
            let m = total_mass(&masses, &final_state);
            assert!((m - m0).abs() < 1e-9, "mass {m} vs {m0}");
        }
    }
}
