//! The Boys function `F_n(T)`, the workhorse of Gaussian integral
//! evaluation: `F_n(T) = ∫₀¹ t^{2n} e^{−Tt²} dt`.
//!
//! Three regimes, each numerically stable (see Obara & Saika 1986 and
//! standard references):
//!
//! * `T ≈ 0`: Taylor limit `F_n(0) = 1/(2n+1)`.
//! * small/medium `T ≤ 12`: the convergent power series
//!   `F_n(T) = Σ_k (−T)^k / (k!(2n+2k+1))` (alternating; cancellation
//!   stays harmless below T ≈ 12).
//! * large `T > 12`: upward recurrence
//!   `F_n = [(2n−1)F_{n−1} − e^{−T}]/(2T)` from
//!   `F_0 = ½√(π/T)·erf(√T)` — the division by `2T` shrinks errors, so
//!   the upward direction is stable.

use tpt_chem_core::num;
use tpt_chem_core::special;

/// `F_0(T) = ½√(π/T)·erf(√T)`, using the erf from `tpt-chem-core`.
pub fn boys_0(t: f64) -> f64 {
    if t < 1e-12 {
        return 1.0;
    }
    0.5 * num::sqrt(core::f64::consts::PI / t) * special::erf(num::sqrt(t))
}

/// `F_n(T)` for `n` in `0..=max_n`, returned as a vector of length
/// `max_n + 1`.
pub fn boys_array(t: f64, max_n: usize) -> Vec<f64> {
    let mut out = vec![0.0; max_n + 1];
    if t < 1e-10 {
        // Taylor limit.
        for (n, v) in out.iter_mut().enumerate() {
            *v = 1.0 / (2 * n + 1) as f64;
        }
        return out;
    }
    if t <= 12.0 {
        // Power series per n; alternating terms with harmless cancellation
        // below T ≈ 12 (largest term ~10⁴, keeping ≥ 12 digits).
        for (n, v) in out.iter_mut().enumerate() {
            let nf = n as f64;
            // k = 0 term, then term = (−T)^k / k!.
            let mut sum = 1.0 / (2.0 * nf + 1.0);
            let mut term = 1.0f64;
            for k in 1..200usize {
                term *= -t / (k as f64);
                let contrib = term / (2.0 * (n + k) as f64 + 1.0);
                sum += contrib;
                if contrib.abs() <= 1e-17 * sum.abs() {
                    break;
                }
            }
            *v = sum;
        }
        return out;
    }
    // Large T: F_0 from erf, then upward recurrence
    // F_n = [(2n−1)F_{n−1} − e^{−T}]/(2T) — division by 2T shrinks
    // errors, so the upward direction is numerically stable here.
    out[0] = boys_0(t);
    let exp_t = num::exp(-t);
    for n in 1..=max_n {
        out[n] = ((2.0 * n as f64 - 1.0) * out[n - 1] - exp_t) / (2.0 * t);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn boys_0_reference_values() {
        assert!(close(boys_0(0.0), 1.0, 1e-12));
        // F_0(1) = 0.7468241328124271 (½√π·erf(1)).
        assert!(close(boys_0(1.0), 0.746_824_132_812_427_1, 1e-13));
        // F_0(4) = ½√(π/4)·erf(2), checked against the erf implementation.
        let want4 = 0.5 * (core::f64::consts::PI / 4.0).sqrt() * special::erf(2.0);
        assert!(close(boys_0(4.0), want4, 1e-14));
    }

    #[test]
    fn boys_recurrence_identity() {
        // Integration by parts gives 2T·F_n = (2n−1)·F_{n−1} − e^{−T}
        // for every n ≥ 1.
        for &t in &[0.5, 5.0, 11.9, 12.1, 100.0, 1e4] {
            let f = boys_array(t, 6);
            for n in 1..=6 {
                let lhs = (2.0 * n as f64 - 1.0) * f[n - 1] - num::exp(-t);
                let rhs = 2.0 * t * f[n];
                assert!(
                    (lhs - rhs).abs() < 1e-9 * (1.0 + t),
                    "T = {t}, n = {n}: {lhs} vs {rhs}"
                );
            }
        }
    }

    #[test]
    fn small_t_taylor_limit() {
        let f = boys_array(1e-14, 3);
        assert!(close(f[0], 1.0, 1e-12));
        assert!(close(f[1], 1.0 / 3.0, 1e-12));
        assert!(close(f[2], 0.2, 1e-12));
    }

    #[test]
    fn series_and_recurrence_agree_across_regimes() {
        // Both regimes must reproduce the erf-based F_0 on their own side
        // of the T = 12 boundary.
        for &t in &[11.99, 12.01] {
            let f = boys_array(t, 1);
            assert!(
                (f[0] - boys_0(t)).abs() < 1e-10,
                "T = {t}: {} vs {}",
                f[0],
                boys_0(t)
            );
        }
    }

    #[test]
    fn f1_reference_value() {
        // F_1(0.5) from the analytic integral = 0.24909373217951547.
        assert!(close(
            boys_array(0.5, 1)[1],
            0.249_093_732_179_515_47,
            1e-13
        ));
        // F_1(1) = 0.18947232570982622.
        assert!(close(
            boys_array(1.0, 1)[1],
            0.189_472_345_820_492_33,
            1e-13
        ));
    }
}
