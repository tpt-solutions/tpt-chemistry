//! Michaelis–Menten kinetics: deterministic mass-action ODE vs Gillespie
//! SSA at moderate copy numbers.
//!
//! Run: `cargo run -p tpt-chemistry --features kinetics --example michaelis_menten`

use tpt_chemistry::core::rng::Rng;
use tpt_chemistry::kinetics::network::ReactionNetwork;
use tpt_chemistry::kinetics::solver::{integrate, ssa};

fn main() {
    // E + S ⇌ ES → E + P (standard Michaelis–Menten mechanism).
    let net = ReactionNetwork::builder()
        .species(&["E", "S", "ES", "P"])
        .reaction(&[("E", 1.0), ("S", 1.0)], &[("ES", 1.0)], 1.0)
        .reaction(&[("ES", 1.0)], &[("E", 1.0), ("S", 1.0)], 0.5)
        .reaction(&[("ES", 1.0)], &[("E", 1.0), ("P", 1.0)], 0.1)
        .build();

    let y0 = [30.0, 120.0, 0.0, 0.0];
    let t_end = 100.0;

    // Deterministic.
    let traj = integrate(&net, &y0, 0.0, t_end, 2_000);
    let det = &traj[traj.len() - 1];
    println!("Michaelis-Menten to t = {t_end} (E0 = 30, S0 = 120):");
    println!("deterministic: P = {:.2}  (E free = {:.2})", det[3], det[0]);

    // Stochastic, ensemble mean.
    const RUNS: usize = 1_000;
    let mut p_sum = 0.0;
    for s in 0..RUNS {
        let mut rng = Rng::new(7 + s as u64);
        let final_state = ssa(&net, &y0, t_end, &mut rng);
        p_sum += final_state[3];
    }
    println!("SSA mean ({RUNS} runs): P = {:.2}", p_sum / RUNS as f64);
    println!(
        "\n(at these copy numbers both regimes agree to sampling noise; \
         shrink E0/S0 by 10× to watch the trajectories separate)"
    );
}
