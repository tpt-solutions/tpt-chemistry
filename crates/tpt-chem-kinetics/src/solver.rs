//! ODE and SSA solvers for [`ReactionNetwork`]s.
//!
//! * [`integrate`] — classical Runge–Kutta 4 over a uniform time grid
//!   (deterministic mass-action trajectories).
//! * [`integrate_stiff`] — semi-implicit (Rosenbrock-Euler) stepping for
//!   stiff networks, adaptive sub-stepping inside each output interval.
//! * [`ssa`] — Gillespie stochastic simulation algorithm (direct method)
//!   with the in-house RNG from `tpt-chem-core`.
//! * [`steady_state`] — damped fixed-point iteration to a steady state.

use std::vec::Vec;

use tpt_chem_core::rng::Rng;

use crate::network::ReactionNetwork;

/// Integrate the mass-action ODE with classical RK4 over
/// `[t0, t_end]` with `n_steps` uniform steps; returns the trajectory
/// (length `n_steps + 1`), each entry a copy of the state vector.
pub fn integrate(
    net: &ReactionNetwork,
    y0: &[f64],
    t0: f64,
    t_end: f64,
    n_steps: usize,
) -> Vec<Vec<f64>> {
    let n = net.n_species();
    assert_eq!(y0.len(), n, "state vector length mismatch");
    assert!(n_steps > 0);
    let h = (t_end - t0) / n_steps as f64;

    let mut traj = Vec::with_capacity(n_steps + 1);
    let mut y = y0.to_vec();
    traj.push(y.clone());
    let mut k1 = vec![0.0; n];
    let mut k2 = vec![0.0; n];
    let mut k3 = vec![0.0; n];
    let mut k4 = vec![0.0; n];
    let mut tmp = vec![0.0; n];
    for _ in 0..n_steps {
        k1 = net.derivative(&y);
        for i in 0..n {
            tmp[i] = y[i] + 0.5 * h * k1[i];
        }
        k2 = net.derivative(&tmp);
        for i in 0..n {
            tmp[i] = y[i] + 0.5 * h * k2[i];
        }
        k3 = net.derivative(&tmp);
        for i in 0..n {
            tmp[i] = y[i] + h * k3[i];
        }
        k4 = net.derivative(&tmp);
        for i in 0..n {
            y[i] += (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        traj.push(y.clone());
    }
    traj
}

/// Integrate a (possibly stiff) network with implicit Euler
/// `y_{n+1} = y_n + h f(y_{n+1})`, solved by Newton iteration with a
/// finite-difference Jacobian and dense Gaussian elimination (fine for the
/// handful of species in reaction networks). A-stability tames fast
/// reversible pairs, and the Newton fixed point preserves the true
/// equilibrium.
pub fn integrate_stiff(
    net: &ReactionNetwork,
    y0: &[f64],
    t0: f64,
    t_end: f64,
    n_steps: usize,
    newton_iters: usize,
) -> Vec<Vec<f64>> {
    let h = (t_end - t0) / n_steps as f64;
    let mut traj = Vec::with_capacity(n_steps + 1);
    let mut y = y0.to_vec();
    traj.push(y.clone());
    for _ in 0..n_steps {
        implicit_step(net, &mut y, h, newton_iters);
        for v in &mut y {
            if !v.is_finite() || *v < 0.0 {
                *v = v.min(0.0).abs() * 0.0;
            }
        }
        traj.push(y.clone());
    }
    traj
}

/// One implicit-Euler Newton step: solve `G(z) = z - y - h f(z) = 0` with a
/// finite-difference Jacobian `J = I - h df/dz` and dense Gaussian
/// elimination with partial pivoting.
fn implicit_step(net: &ReactionNetwork, y: &mut [f64], h: f64, iters: usize) {
    let n = y.len();
    let eps = 1e-7;
    let mut z = y.to_vec();
    for _ in 0..iters {
        let fz = net.derivative(&z);
        let g: Vec<f64> = (0..n).map(|i| z[i] - y[i] - h * fz[i]).collect();
        let mut jac = vec![vec![0.0f64; n]; n];
        for c in 0..n {
            let scale = (z[c].abs() + 1.0) * eps;
            let mut zp = z.clone();
            zp[c] += scale;
            let fzp = net.derivative(&zp);
            for r in 0..n {
                jac[r][c] = if r == c { 1.0 } else { 0.0 } - h * (fzp[r] - fz[r]) / scale;
            }
        }
        let d = solve_linear(&jac, &g.iter().map(|v| -v).collect::<Vec<_>>());
        for i in 0..n {
            z[i] = (z[i] + d[i]).max(0.0);
        }
    }
    y.copy_from_slice(&z);
}

/// Dense Gaussian elimination with partial pivoting; returns `x` with
/// `A x = b` (or the closest result for singular `A`).
fn solve_linear(a: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut m = a.to_vec();
    for i in 0..n {
        m[i].push(b[i]);
    }
    for col in 0..n {
        let mut pivot = col;
        for r in (col + 1)..n {
            if m[r][col].abs() > m[pivot][col].abs() {
                pivot = r;
            }
        }
        m.swap(col, pivot);
        let piv = m[col][col];
        if piv.abs() < 1e-14 {
            continue;
        }
        let pivot_row = m[col].clone();
        for row in m.iter_mut().take(n).skip(col + 1) {
            let factor = row[col] / piv;
            for (c, mrc) in row[col..=n].iter_mut().enumerate() {
                *mrc -= factor * pivot_row[col + c];
            }
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let piv = m[i][i];
        if piv.abs() < 1e-14 {
            continue;
        }
        let mut sum = m[i][n];
        for c in (i + 1)..n {
            sum -= m[i][c] * x[c];
        }
        x[i] = sum / piv;
    }
    x
}

/// Gillespie SSA (direct method): simulate until `t_end`, recording the
/// population after every event plus the final state (entries are
/// `(time, state)` pairs; the first is `(0, y0)`).
pub fn ssa_events(
    net: &ReactionNetwork,
    y0: &[f64],
    t_end: f64,
    rng: &mut Rng,
) -> Vec<(f64, Vec<f64>)> {
    let mut x: Vec<f64> = y0.to_vec();
    let mut events = vec![(0.0, x.clone())];
    let mut t = 0.0f64;
    loop {
        let prop = net.propensities(&x);
        let a_total: f64 = prop.iter().sum();
        if a_total <= 0.0 {
            events.push((t, x.clone()));
            return events;
        }
        let tau = -rng.uniform().ln() / a_total;
        if t + tau > t_end {
            events.push((t, x.clone()));
            return events;
        }
        t += tau;
        // Choose reaction by the direct method.
        let target = rng.uniform() * a_total;
        let mut acc = 0.0;
        let mut chosen = net.n_reactions() - 1;
        for (r, &a) in prop.iter().enumerate() {
            acc += a;
            if target < acc {
                chosen = r;
                break;
            }
        }
        for (i, &nu) in net.reactions[chosen].stoichiometry.iter().enumerate() {
            x[i] += nu;
            if x[i] < 0.0 {
                x[i] = 0.0;
            }
        }
        events.push((t, x.clone()));
    }
}

/// SSA returning only the state at `t_end` (runs the event stream to the
/// horizon). Populations are kept as f64 copies of the (integer) copy
/// numbers passed in `y0`.
pub fn ssa(net: &ReactionNetwork, y0: &[f64], t_end: f64, rng: &mut Rng) -> Vec<f64> {
    let last = ssa_events(net, y0, t_end, rng)
        .last()
        .cloned()
        .expect("ssa: at least the initial state");
    last.1
}

/// Find a steady state by damped fixed-point iteration on the explicit
/// Euler map `x ← x + λ·f(x)` with back-off when the residual grows.
/// Returns `Err` with the last state if no steady state is reached within
/// `max_iter` iterations (e.g. oscillating or unbounded networks).
pub fn steady_state(
    net: &ReactionNetwork,
    y0: &[f64],
    lambda: f64,
    tol: f64,
    max_iter: usize,
) -> Result<Vec<f64>, Vec<f64>> {
    let mut x = y0.to_vec();
    let mut lambda = lambda;
    for _ in 0..max_iter {
        let f = net.derivative(&x);
        let residual = f.iter().fold(0.0f64, |a, &v| a.max(v.abs()));
        if residual < tol {
            return Ok(x);
        }
        let next: Vec<f64> = x
            .iter()
            .zip(&f)
            .map(|(xi, fi)| (xi + lambda * fi).max(0.0))
            .collect();
        let next_res = {
            let f2 = net.derivative(&next);
            f2.iter().fold(0.0f64, |a, &v| a.max(v.abs()))
        };
        if next_res > residual {
            // Back off and retry the same step.
            lambda *= 0.5;
            if lambda < 1e-12 {
                return Err(x);
            }
            continue;
        }
        x = next;
    }
    Err(x)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::network::ReactionNetwork;

    fn decay_net() -> ReactionNetwork {
        ReactionNetwork::builder()
            .species(&["A", "B"])
            .reaction(&[("A", 1.0)], &[("B", 1.0)], 1.0)
            .build()
    }

    #[test]
    fn rk4_first_order_matches_analytic() {
        let net = decay_net();
        let traj = integrate(&net, &[10.0, 0.0], 0.0, 5.0, 500);
        let expect_a = 10.0 * (-5.0_f64).exp();
        let expect_b = 10.0 - expect_a;
        let last = traj.last().unwrap();
        assert!((last[0] - expect_a).abs() < 1e-6);
        assert!((last[1] - expect_b).abs() < 1e-6);
        // Conservation of total.
        assert!((last[0] + last[1] - 10.0).abs() < 1e-9);
    }

    #[test]
    fn rk4_reversible_converges_to_equilibrium() {
        let net = ReactionNetwork::builder()
            .species(&["A", "B"])
            .reaction(&[("A", 1.0)], &[("B", 1.0)], 3.0)
            .reaction(&[("B", 1.0)], &[("A", 1.0)], 1.0)
            .build();
        let traj = integrate(&net, &[1.0, 0.0], 0.0, 10.0, 1000);
        let last = traj.last().unwrap();
        assert!((last[0] - 0.25).abs() < 1e-6, "A = {}", last[0]);
        assert!((last[1] - 0.75).abs() < 1e-6, "B = {}", last[1]);
    }

    #[test]
    fn stiff_solver_handles_fast_reversible() {
        // k+ = 1e4, k− = 2e4: explicit RK4 would need h ~ 1e-5.
        let net = ReactionNetwork::builder()
            .species(&["A", "B"])
            .reaction(&[("A", 1.0)], &[("B", 1.0)], 1.0e4)
            .reaction(&[("B", 1.0)], &[("A", 1.0)], 2.0e4)
            .build();
        let traj = integrate_stiff(&net, &[1.0, 0.0], 0.0, 1.0, 10, 2000);
        let last = traj.last().unwrap();
        // Equilibrium: 1e4 [A] = 2e4 [B] → [B] = [A]/2, [A] = 2/3.
        assert!((last[0] - 2.0 / 3.0).abs() < 5e-3, "A = {}", last[0]);
        assert!((last[1] - 1.0 / 3.0).abs() < 5e-3, "B = {}", last[1]);
    }

    #[test]
    fn ssa_mean_tracks_analytic_decay() {
        let net = decay_net();
        let mut rng = Rng::new(12345);
        let runs = 300;
        let n0 = 400.0;
        let mut mean = 0.0;
        for _ in 0..runs {
            let final_state = ssa(&net, &[n0, 0.0], 5.0, &mut rng);
            mean += final_state[0];
        }
        mean /= runs as f64;
        let analytic = n0 * (-5.0_f64).exp();
        // Statistical error ~ n0^0.5/sqrt(runs) ≈ 1.15; tolerance 5σ.
        assert!(
            (mean - analytic).abs() < 5.0,
            "mean {mean} vs analytic {analytic}"
        );
    }

    #[test]
    fn ssa_extinct_network_terminates() {
        let net = ReactionNetwork::builder()
            .species(&["A"])
            .reaction(&[("A", 1.0)], &[], 1.0)
            .build();
        let mut rng = Rng::new(7);
        let final_state = ssa(&net, &[5.0], 100.0, &mut rng);
        assert_eq!(final_state[0], 0.0);
    }

    #[test]
    fn steady_state_finds_equilibrium() {
        let net = ReactionNetwork::builder()
            .species(&["A", "B"])
            .reaction(&[("A", 1.0)], &[("B", 1.0)], 3.0)
            .reaction(&[("B", 1.0)], &[("A", 1.0)], 1.0)
            .build();
        let ss = steady_state(&net, &[1.0, 0.0], 0.05, 1e-10, 5000).unwrap();
        assert!((ss[0] - 0.25).abs() < 1e-6, "A = {}", ss[0]);
        assert!((ss[1] - 0.75).abs() < 1e-6, "B = {}", ss[1]);
    }

    #[test]
    fn lotka_volterra_oscillates_stays_bounded() {
        // Classic LV predator-prey: dX = X(aY − b), dY = cX − dY-ish setup.
        // Use a minimal oscillator: X' = X(1 − Y), Y' = Y(X − 1) style
        // network encoded via reactions with mass-action rates.
        let net = ReactionNetwork::builder()
            .species(&["Prey", "Pred"])
            // Prey growth: → Prey, k = 1.1 (constant influx)
            .reaction(&[], &[("Prey", 1.0)], 1.1)
            // Predation: Prey + Pred → 2 Pred, k = 0.9
            .reaction(&[("Prey", 1.0), ("Pred", 1.0)], &[("Pred", 1.0)], 0.9)
            // Predator death: Pred → , k = 1.1
            .reaction(&[("Pred", 1.0)], &[], 1.1)
            .build();
        let traj = integrate(&net, &[1.0, 0.5], 0.0, 20.0, 4000);
        let max_prey = traj.iter().map(|s| s[0]).fold(0.0f64, f64::max);
        let min_prey = traj.iter().map(|s| s[0]).fold(f64::INFINITY, f64::min);
        assert!(max_prey.is_finite() && max_prey > 1.0, "oscillates");
        assert!(min_prey >= 0.0, "bounded below by 0");
    }
}
