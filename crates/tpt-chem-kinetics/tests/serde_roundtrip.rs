//! Serde round trip for kinetics networks and fit results (requires the
//! `serde` feature).

#![cfg(feature = "serde")]

use tpt_chem_kinetics::fit::{fit_rates, FitResult, FitSettings, Observation};
use tpt_chem_kinetics::network::ReactionNetwork;
use tpt_chem_kinetics::solver::integrate;

#[test]
fn network_round_trips_and_simulates_identically() {
    let net = ReactionNetwork::builder()
        .species(&["A", "B"])
        .reaction(&[("A", 1.0)], &[("B", 1.0)], 0.7)
        .build();
    let json = serde_json::to_string(&net).unwrap();
    let back: ReactionNetwork = serde_json::from_str(&json).unwrap();
    let y0 = [2.0, 0.0];
    assert_eq!(
        integrate(&net, &y0, 0.0, 3.0, 300),
        integrate(&back, &y0, 0.0, 3.0, 300)
    );
}

#[test]
fn fit_result_round_trips() {
    let net = ReactionNetwork::builder()
        .species(&["A", "B"])
        .reaction(&[("A", 1.0)], &[("B", 1.0)], 0.5)
        .build();
    let obs: Vec<Observation> = [0.5, 1.0, 2.0, 3.0]
        .iter()
        .map(|&t| Observation {
            time: t,
            species: 0,
            value: (-0.8_f64 * t).exp(),
        })
        .collect();
    let fit = fit_rates(&net, &[1.0, 0.0], &obs, &[0.5], &FitSettings::default()).unwrap();
    let back: FitResult = serde_json::from_str(&serde_json::to_string(&fit).unwrap()).unwrap();
    assert_eq!(back, fit);
    assert!((back.rates[0] - 0.8).abs() < 1e-6);
}
