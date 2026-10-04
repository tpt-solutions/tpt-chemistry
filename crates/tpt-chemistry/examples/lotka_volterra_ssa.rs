//! Predator–prey dynamics (logistic Lotka–Volterra) via Gillespie SSA
//! against the deterministic mass-action ODE.
//!
//! Run: `cargo run -p tpt-chemistry --features kinetics --example lotka_volterra_ssa`

use tpt_chemistry::core::rng::Rng;
use tpt_chemistry::kinetics::network::ReactionNetwork;
use tpt_chemistry::kinetics::solver::{integrate, ssa_events};

fn main() {
    // Logistic Lotka–Volterra (bounded, SSA-friendly):
    //   A → 2A         (prey births ∝ A, α)
    //   2A → ∅         (competition, κ)
    //   A + B → 2B     (predation, γ)
    //   B → ∅          (predator death, δ)
    // Fixed point: prey A* = δ/γ = 50, predator B* = (α − 2κA*)/γ = 75
    // (the factor 2: each competition event consumes two prey) — a damped
    // spiral, so orbits stay bounded. The pure LV system (κ = 0) is
    // neutrally stable: SSA noise diffuses to arbitrarily large orbits and
    // event counts explode; the competition term is what makes long
    // stochastic runs tractable.
    let (alpha, kappa, gamma, delta) = (1.0, 0.0025, 0.01, 0.5);
    let net = ReactionNetwork::builder()
        .species(&["prey", "predator"])
        .reaction(&[("prey", 1.0)], &[("prey", 2.0)], alpha)
        .reaction(&[("prey", 2.0)], &[], kappa)
        .reaction(
            &[("prey", 1.0), ("predator", 1.0)],
            &[("predator", 2.0)],
            gamma,
        )
        .reaction(&[("predator", 1.0)], &[], delta)
        .build();

    let y0 = [50.0, 75.0];
    let t_end = 50.0;

    // Deterministic reference.
    let traj = integrate(&net, &y0, 0.0, t_end, 500);
    let det = &traj[traj.len() - 1];
    println!("logistic predator-prey to t = {t_end} (fixed point 50, 75):");
    println!(
        "deterministic: prey = {:.2}, predator = {:.2}",
        det[0], det[1]
    );

    // Stochastic: mean over an ensemble of independent runs.
    const RUNS: usize = 500;
    let mut sum = [0.0f64; 2];
    let mut event_counts = Vec::new();
    for s in 0..RUNS {
        let mut rng = Rng::new(1000 + s as u64);
        let events = ssa_events(&net, &y0, t_end, &mut rng);
        event_counts.push(events.len() as f64);
        let final_state = &events[events.len() - 1].1;
        sum[0] += final_state[0];
        sum[1] += final_state[1];
    }
    println!(
        "SSA mean over {RUNS} runs: prey = {:.2}, predator = {:.2}          (mean event count {:.0})",
        sum[0] / RUNS as f64,
        sum[1] / RUNS as f64,
        event_counts.iter().sum::<f64>() / event_counts.len() as f64
    );
    println!(
        "
(the SSA ensemble scatters around the damped deterministic          spiral with 1/sqrt(N) sampling noise)"
    );
}
