//! Rate-constant fitting and local sensitivity analysis for mass-action
//! networks.
//!
//! [`fit_rates`] estimates every reaction's rate constant from observed
//! species concentrations by unweighted least squares. It uses
//! Levenberg–Marquardt on `θ = ln k`, so rates stay positive. The model is
//! integrated with classical RK4, and the Jacobian comes from central finite
//! differences. [`sensitivity`] reports the normalized local sensitivity
//! `k_j ∂y_s(t)/∂k_j`, using the same integrator.
//!
//! Fitting replaces every reaction's rate input with a constant, including
//! temperature-dependent laws. Fit at one temperature, then reapply the
//! temperature law if needed.
//!
//! # Examples
//!
//! Recover the rate constant of A → B from noise-free data:
//!
//! ```
//! use tpt_chem_kinetics::fit::{fit_rates, FitSettings, Observation};
//! use tpt_chem_kinetics::network::ReactionNetwork;
//!
//! let net = ReactionNetwork::builder()
//!     .species(&["A", "B"])
//!     .reaction(&[("A", 1.0)], &[("B", 1.0)], 1.0)
//!     .build();
//! let y0 = [1.0, 0.0];
//! let obs: Vec<Observation> = (1..=5)
//!     .map(|i| {
//!         let t = f64::from(i) * 0.5;
//!         Observation { time: t, species: 0, value: (-0.8 * t).exp() }
//!     })
//!     .collect();
//! let fit = fit_rates(&net, &y0, &obs, &[0.3], &FitSettings::default()).unwrap();
//! assert!((fit.rates[0] - 0.8).abs() < 1e-6);
//! ```

use crate::network::{ReactionInput, ReactionNetwork};

/// One measured concentration.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Observation {
    /// Time of the measurement (≥ 0, same units as the rate constants).
    pub time: f64,
    /// Index of the observed species in the network.
    pub species: usize,
    /// Measured concentration.
    pub value: f64,
}

/// Tuning parameters for [`fit_rates`] and [`sensitivity`].
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FitSettings {
    /// Maximum Levenberg–Marquardt iterations.
    pub max_iterations: usize,
    /// Stop when the largest log-rate step falls below this.
    pub tolerance: f64,
    /// Largest RK4 step used to integrate between observation times.
    pub h_max: f64,
    /// Central-difference step in `ln k` for Jacobians and sensitivities.
    pub fd_step: f64,
}

impl Default for FitSettings {
    fn default() -> Self {
        Self {
            max_iterations: 100,
            tolerance: 1e-10,
            h_max: 0.01,
            fd_step: 1e-6,
        }
    }
}

/// Result of [`fit_rates`].
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FitResult {
    /// Fitted rate constants, one per reaction, in network order.
    pub rates: Vec<f64>,
    /// Sum of squared residuals at the fitted rates.
    pub sse: f64,
    /// Iterations performed.
    pub iterations: usize,
    /// `true` if the log-rate step fell below the tolerance.
    pub converged: bool,
    /// Standard error of `ln k` for each reaction, from the linearized
    /// covariance `σ² (JᵀJ)⁻¹` with `σ² = SSE / (m − n)`. A relative error
    /// on `k` is approximately this value. `None` when there are no spare
    /// observations (`m ≤ n`) or the normal matrix is singular.
    pub log_std_errors: Option<Vec<f64>>,
}

/// Errors from fitting and sensitivity analysis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FitError {
    /// An input was invalid (empty data, bad index, non-positive rate, …).
    InvalidInput(&'static str),
    /// The normal equations were singular, so the parameters are not all
    /// identifiable from the data.
    Singular,
}

impl core::fmt::Display for FitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            FitError::InvalidInput(m) => write!(f, "invalid input: {m}"),
            FitError::Singular => write!(f, "normal equations are singular"),
        }
    }
}

impl std::error::Error for FitError {}

/// A copy of `net` in which reaction `j` uses constant rate `rates[j]`.
///
/// # Panics
/// Panics if `rates.len() != net.reactions.len()`.
pub fn with_rates(net: &ReactionNetwork, rates: &[f64]) -> ReactionNetwork {
    assert_eq!(rates.len(), net.reactions.len(), "one rate per reaction");
    let mut out = net.clone();
    for (reaction, &k) in out.reactions.iter_mut().zip(rates) {
        reaction.rate = ReactionInput::Constant(k);
    }
    out
}

/// State at each of `times` (sorted ascending, all ≥ 0), starting from
/// `y0` at `t = 0`. Each interval is split into RK4 steps no longer than
/// `h_max`.
pub fn trajectory_at(
    net: &ReactionNetwork,
    y0: &[f64],
    times: &[f64],
    h_max: f64,
) -> Vec<Vec<f64>> {
    let mut out = Vec::with_capacity(times.len());
    let mut y = y0.to_vec();
    let mut t = 0.0;
    for &target in times {
        let span = target - t;
        if span > 0.0 {
            let steps = (span / h_max).ceil().max(1.0) as usize;
            let h = span / steps as f64;
            for _ in 0..steps {
                rk4_step(net, &mut y, h);
            }
            t = target;
        }
        out.push(y.clone());
    }
    out
}

fn rk4_step(net: &ReactionNetwork, y: &mut [f64], h: f64) {
    let n = y.len();
    let k1 = net.derivative(y);
    let tmp: Vec<f64> = (0..n).map(|i| y[i] + 0.5 * h * k1[i]).collect();
    let k2 = net.derivative(&tmp);
    let tmp: Vec<f64> = (0..n).map(|i| y[i] + 0.5 * h * k2[i]).collect();
    let k3 = net.derivative(&tmp);
    let tmp: Vec<f64> = (0..n).map(|i| y[i] + h * k3[i]).collect();
    let k4 = net.derivative(&tmp);
    for i in 0..n {
        y[i] += (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
}

/// Validated observation data: sorted unique times and, per observation,
/// its index into those times.
struct Data {
    times: Vec<f64>,
    time_index: Vec<usize>,
}

fn prepare(net: &ReactionNetwork, y0: &[f64], obs: &[Observation]) -> Result<Data, FitError> {
    if y0.len() != net.n_species() {
        return Err(FitError::InvalidInput(
            "y0 length must equal the species count",
        ));
    }
    if obs.is_empty() {
        return Err(FitError::InvalidInput("no observations"));
    }
    for o in obs {
        if o.species >= net.n_species() {
            return Err(FitError::InvalidInput(
                "observation species index out of range",
            ));
        }
        if !(o.time.is_finite() && o.time >= 0.0 && o.value.is_finite()) {
            return Err(FitError::InvalidInput(
                "observation times and values must be finite, times ≥ 0",
            ));
        }
    }
    let mut times: Vec<f64> = obs.iter().map(|o| o.time).collect();
    times.sort_by(|a, b| a.total_cmp(b));
    times.dedup();
    let time_index = obs
        .iter()
        .map(|o| {
            times
                .binary_search_by(|t| t.total_cmp(&o.time))
                .expect("time present")
        })
        .collect();
    Ok(Data { times, time_index })
}

fn residuals(
    base: &ReactionNetwork,
    y0: &[f64],
    data: &Data,
    obs: &[Observation],
    theta: &[f64],
    h_max: f64,
) -> Vec<f64> {
    let rates: Vec<f64> = theta.iter().map(|t| t.exp()).collect();
    let net = with_rates(base, &rates);
    let traj = trajectory_at(&net, y0, &data.times, h_max);
    obs.iter()
        .zip(&data.time_index)
        .map(|(o, &ti)| traj[ti][o.species] - o.value)
        .collect()
}

fn sum_squares(r: &[f64]) -> f64 {
    r.iter().map(|x| x * x).sum()
}

/// Central-difference Jacobian of `f` at `theta`, `m × n` row-major.
fn jacobian<F: Fn(&[f64]) -> Vec<f64>>(f: &F, theta: &[f64], step: f64) -> Vec<Vec<f64>> {
    let n = theta.len();
    let mut cols: Vec<Vec<f64>> = Vec::with_capacity(n);
    for j in 0..n {
        let mut plus = theta.to_vec();
        let mut minus = theta.to_vec();
        plus[j] += step;
        minus[j] -= step;
        let rp = f(&plus);
        let rm = f(&minus);
        cols.push(
            rp.iter()
                .zip(&rm)
                .map(|(a, b)| (a - b) / (2.0 * step))
                .collect(),
        );
    }
    let m = cols.first().map_or(0, Vec::len);
    (0..m)
        .map(|i| cols.iter().map(|c| c[i]).collect())
        .collect()
}

/// Solves `a x = b` by Gaussian elimination with partial pivoting.
fn solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Result<Vec<f64>, FitError> {
    let n = b.len();
    for col in 0..n {
        let pivot = (col..n)
            .max_by(|&i, &j| a[i][col].abs().total_cmp(&a[j][col].abs()))
            .expect("non-empty range");
        if a[pivot][col].abs() < 1e-300 {
            return Err(FitError::Singular);
        }
        a.swap(col, pivot);
        b.swap(col, pivot);
        let pivot_row = a[col].clone();
        for row in col + 1..n {
            let f = a[row][col] / pivot_row[col];
            for (x, p) in a[row][col..].iter_mut().zip(&pivot_row[col..]) {
                *x -= f * p;
            }
            b[row] -= f * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let s: f64 = (row + 1..n).map(|k| a[row][k] * x[k]).sum();
        x[row] = (b[row] - s) / a[row][row];
    }
    Ok(x)
}

/// Least-squares estimate of every reaction's rate constant.
///
/// `initial` gives one positive starting rate per reaction. Levenberg–
/// Marquardt minimizes the unweighted sum of squared residuals between the
/// integrated model and `obs`.
///
/// # Errors
/// [`FitError::InvalidInput`] for mismatched or non-finite inputs, and
/// [`FitError::Singular`] when the normal equations cannot be solved.
pub fn fit_rates(
    net: &ReactionNetwork,
    y0: &[f64],
    obs: &[Observation],
    initial: &[f64],
    settings: &FitSettings,
) -> Result<FitResult, FitError> {
    if initial.len() != net.reactions.len() {
        return Err(FitError::InvalidInput("one initial rate per reaction"));
    }
    if initial.iter().any(|k| !(k.is_finite() && *k > 0.0)) {
        return Err(FitError::InvalidInput("initial rates must be positive"));
    }
    let data = prepare(net, y0, obs)?;
    let h = settings.h_max;
    let r_of = |theta: &[f64]| residuals(net, y0, &data, obs, theta, h);

    let mut theta: Vec<f64> = initial.iter().map(|k| k.ln()).collect();
    let mut r = r_of(&theta);
    let mut sse = sum_squares(&r);
    let mut lambda = 1e-3;
    let mut converged = false;
    let mut iterations = 0;

    for iter in 1..=settings.max_iterations {
        iterations = iter;
        let jac = jacobian(&r_of, &theta, settings.fd_step);
        let n = theta.len();
        let mut jtj = vec![vec![0.0; n]; n];
        let mut g = vec![0.0; n];
        for (row, res) in jac.iter().zip(&r) {
            for a in 0..n {
                g[a] += row[a] * res;
                for b in 0..n {
                    jtj[a][b] += row[a] * row[b];
                }
            }
        }

        // Marquardt damping: scale the diagonal, so a small λ is a
        // Gauss–Newton step and a large λ is a scaled gradient step.
        let mut accepted = None;
        while lambda < 1e12 {
            let mut a = jtj.clone();
            for (i, row) in a.iter_mut().enumerate() {
                row[i] += lambda * jtj[i][i].max(1e-12);
            }
            let step = solve(a, g.iter().map(|v| -v).collect())?;
            let trial: Vec<f64> = theta.iter().zip(&step).map(|(t, s)| t + s).collect();
            let r_trial = r_of(&trial);
            let sse_trial = sum_squares(&r_trial);
            if sse_trial < sse {
                accepted = Some((trial, r_trial, sse_trial, step));
                lambda = (lambda / 10.0).max(1e-12);
                break;
            }
            lambda *= 10.0;
        }

        match accepted {
            Some((trial, r_trial, sse_trial, step)) => {
                theta = trial;
                r = r_trial;
                sse = sse_trial;
                let largest = step.iter().fold(0.0_f64, |m, s| m.max(s.abs()));
                if largest < settings.tolerance {
                    converged = true;
                    break;
                }
            }
            // No damping level reduced the error: at a minimum (or a flat
            // valley), so stop.
            None => {
                converged = true;
                break;
            }
        }
    }

    let log_std_errors = log_standard_errors(&r_of, &theta, sse, obs.len(), settings.fd_step);
    Ok(FitResult {
        rates: theta.iter().map(|t| t.exp()).collect(),
        sse,
        iterations,
        converged,
        log_std_errors,
    })
}

/// Standard errors of the log-rates at `theta` from the linearized
/// covariance. See [`FitResult::log_std_errors`].
fn log_standard_errors<F: Fn(&[f64]) -> Vec<f64>>(
    f: &F,
    theta: &[f64],
    sse: f64,
    m_obs: usize,
    step: f64,
) -> Option<Vec<f64>> {
    let n = theta.len();
    if m_obs <= n {
        return None;
    }
    let sigma2 = sse / (m_obs - n) as f64;
    let jac = jacobian(f, theta, step);
    let mut jtj = vec![vec![0.0; n]; n];
    for row in &jac {
        for a in 0..n {
            for b in 0..n {
                jtj[a][b] += row[a] * row[b];
            }
        }
    }
    // Diagonal of (JᵀJ)⁻¹: solve against each unit vector in turn.
    let mut errors = Vec::with_capacity(n);
    for j in 0..n {
        let mut unit = vec![0.0; n];
        unit[j] = 1.0;
        let col = solve(jtj.clone(), unit).ok()?;
        errors.push((sigma2 * col[j]).max(0.0).sqrt());
    }
    Some(errors)
}

/// Normalized local sensitivities of one species' concentration.
///
/// Returns `s[m][j] = k_j ∂y_species(times[m]) / ∂k_j` for the reactions of
/// `net` at its current rates, where `k_j` is the reaction's rate
/// constant at `T = 300 K`. Each entry says how much a relative change in
/// `k_j` moves the concentration. `times` must be sorted ascending.
///
/// # Errors
/// [`FitError::InvalidInput`] for an unsorted or negative time grid, an
/// out-of-range species index, or a non-positive rate constant.
pub fn sensitivity(
    net: &ReactionNetwork,
    y0: &[f64],
    times: &[f64],
    species: usize,
    settings: &FitSettings,
) -> Result<Vec<Vec<f64>>, FitError> {
    if species >= net.n_species() {
        return Err(FitError::InvalidInput("species index out of range"));
    }
    if times.windows(2).any(|w| w[0] > w[1]) || times.iter().any(|t| *t < 0.0) {
        return Err(FitError::InvalidInput(
            "times must be sorted and non-negative",
        ));
    }
    let rates: Vec<f64> = net
        .reactions
        .iter()
        .map(|r| match &r.rate {
            ReactionInput::Constant(k) => *k,
            ReactionInput::TemperatureDependent(law) => law.rate_at_300k(),
        })
        .collect();
    if rates.iter().any(|k| !(k.is_finite() && *k > 0.0)) {
        return Err(FitError::InvalidInput("rate constants must be positive"));
    }
    let theta: Vec<f64> = rates.iter().map(|k| k.ln()).collect();
    let y_of = |th: &[f64]| -> Vec<f64> {
        let k: Vec<f64> = th.iter().map(|t| t.exp()).collect();
        let traj = trajectory_at(&with_rates(net, &k), y0, times, settings.h_max);
        traj.iter().map(|state| state[species]).collect()
    };
    // Each column is ∂y/∂ln k_j, which is already the normalized sensitivity.
    let jac = jacobian(&y_of, &theta, settings.fd_step);
    Ok(jac)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn decay(k: f64) -> ReactionNetwork {
        ReactionNetwork::builder()
            .species(&["A", "B"])
            .reaction(&[("A", 1.0)], &[("B", 1.0)], k)
            .build()
    }

    #[test]
    fn sensitivity_matches_analytic_decay() {
        // A(t) = e^{-kt}, so k ∂A/∂k = -k t e^{-kt}.
        let k = 0.8;
        let net = decay(k);
        let times = [0.5, 1.0, 2.0, 3.0];
        let s = sensitivity(&net, &[1.0, 0.0], &times, 0, &FitSettings::default()).unwrap();
        for (row, &t) in s.iter().zip(&times) {
            let exact = -k * t * (-k * t).exp();
            assert!(
                (row[0] - exact).abs() < 1e-5,
                "t={t}: {} vs {exact}",
                row[0]
            );
        }
    }

    #[test]
    fn recovers_two_consecutive_rates() {
        // A → B → C with k1 = 0.8, k2 = 0.3.
        let truth = ReactionNetwork::builder()
            .species(&["A", "B", "C"])
            .reaction(&[("A", 1.0)], &[("B", 1.0)], 0.8)
            .reaction(&[("B", 1.0)], &[("C", 1.0)], 0.3)
            .build();
        let y0 = [1.0, 0.0, 0.0];
        let times: Vec<f64> = (1..=12).map(|i| f64::from(i) * 0.5).collect();
        let model = trajectory_at(&truth, &y0, &times, 0.01);
        let mut obs = Vec::new();
        for (state, &t) in model.iter().zip(&times) {
            obs.push(Observation {
                time: t,
                species: 1,
                value: state[1],
            });
            obs.push(Observation {
                time: t,
                species: 2,
                value: state[2],
            });
        }
        let net = ReactionNetwork::builder()
            .species(&["A", "B", "C"])
            .reaction(&[("A", 1.0)], &[("B", 1.0)], 0.4)
            .reaction(&[("B", 1.0)], &[("C", 1.0)], 0.5)
            .build();
        let fit = fit_rates(&net, &y0, &obs, &[0.4, 0.5], &FitSettings::default()).unwrap();
        assert!(fit.converged);
        assert!((fit.rates[0] - 0.8).abs() < 1e-5, "k1 = {}", fit.rates[0]);
        assert!((fit.rates[1] - 0.3).abs() < 1e-5, "k2 = {}", fit.rates[1]);
        assert!(fit.sse < 1e-10, "sse = {}", fit.sse);
    }

    #[test]
    fn rejects_bad_input() {
        let net = decay(1.0);
        let y0 = [1.0, 0.0];
        let s = FitSettings::default();
        let ok = Observation {
            time: 1.0,
            species: 1,
            value: 0.1,
        };
        assert!(matches!(
            fit_rates(&net, &y0, &[], &[1.0], &s),
            Err(FitError::InvalidInput(_))
        ));
        assert!(matches!(
            fit_rates(&net, &y0, std::slice::from_ref(&ok), &[1.0, 2.0], &s),
            Err(FitError::InvalidInput(_))
        ));
        assert!(matches!(
            fit_rates(&net, &y0, std::slice::from_ref(&ok), &[-1.0], &s),
            Err(FitError::InvalidInput(_))
        ));
        let bad_species = Observation { species: 5, ..ok };
        assert!(matches!(
            fit_rates(&net, &y0, &[bad_species], &[1.0], &s),
            Err(FitError::InvalidInput(_))
        ));
        assert!(sensitivity(&net, &y0, &[2.0, 1.0], 0, &s).is_err());
    }

    #[test]
    fn standard_errors_match_replicate_spread() {
        use tpt_chem_core::rng::Rng;
        let truth = ReactionNetwork::builder()
            .species(&["A", "B", "C"])
            .reaction(&[("A", 1.0)], &[("B", 1.0)], 0.8)
            .reaction(&[("B", 1.0)], &[("C", 1.0)], 0.3)
            .build();
        let y0 = [1.0, 0.0, 0.0];
        let times: Vec<f64> = (1..=12).map(|i| f64::from(i) * 0.5).collect();
        let model = trajectory_at(&truth, &y0, &times, 0.01);
        let net = ReactionNetwork::builder()
            .species(&["A", "B", "C"])
            .reaction(&[("A", 1.0)], &[("B", 1.0)], 0.4)
            .reaction(&[("B", 1.0)], &[("C", 1.0)], 0.5)
            .build();
        let settings = FitSettings::default();
        let (mut fitted, mut reported) = (Vec::new(), Vec::new());
        for seed in 0..60 {
            let mut rng = Rng::new(seed + 1);
            let mut obs = Vec::new();
            for (state, &t) in model.iter().zip(&times) {
                for species in [1, 2] {
                    obs.push(Observation {
                        time: t,
                        species,
                        value: state[species] + rng.normal(0.0, 0.01),
                    });
                }
            }
            let fit = fit_rates(&net, &y0, &obs, &[0.4, 0.5], &settings).unwrap();
            assert!(fit.converged, "replicate {seed} did not converge");
            fitted.push(fit.rates[0].ln());
            reported.push(fit.log_std_errors.expect("24 obs, 2 params")[0]);
        }
        let n = fitted.len() as f64;
        let mean = fitted.iter().sum::<f64>() / n;
        let sd = (fitted.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt();
        let mean_se = reported.iter().sum::<f64>() / n;
        let ratio = mean_se / sd;
        assert!(
            (0.7..1.4).contains(&ratio),
            "reported SE {mean_se} vs spread {sd}"
        );
        // The estimate is unbiased: its mean sits within 3 standard errors
        // of the mean of the replicates.
        assert!(
            (mean - 0.8_f64.ln()).abs() < 3.0 * sd / n.sqrt(),
            "mean ln k1 = {mean}"
        );
    }

    #[test]
    fn no_error_estimate_without_spare_observations() {
        let net = ReactionNetwork::builder()
            .species(&["A", "B", "C"])
            .reaction(&[("A", 1.0)], &[("B", 1.0)], 0.4)
            .reaction(&[("B", 1.0)], &[("C", 1.0)], 0.5)
            .build();
        let y0 = [1.0, 0.0, 0.0];
        let obs = [Observation {
            time: 1.0,
            species: 1,
            value: 0.5,
        }];
        let fit = fit_rates(&net, &y0, &obs, &[0.4, 0.5], &FitSettings::default()).unwrap();
        assert!(fit.log_std_errors.is_none());
    }

    #[test]
    fn with_rates_replaces_constants() {
        let net = decay(1.0);
        let out = with_rates(&net, &[2.5]);
        assert!(matches!(out.reactions[0].rate, ReactionInput::Constant(k) if k == 2.5));
        assert!((out.propensities(&[1.0, 0.0])[0] - 2.5).abs() < 1e-12);
    }
}
