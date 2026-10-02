//! Cubic (order-4) cardinal B-spline values and derivatives — the charge-
//! assignment (interpolation) functions for [`crate::pme`]
//! (Essmann et al., J. Chem. Phys. 103, 8577 (1995)).
//!
//! `M4(w)` is supported on `[0, 4)`, symmetric about `w = 2`, and satisfies
//! the partition of unity `Σ_k M4(u − k) = 1` over integer `k`.

/// Order of the B-spline (number of grid points a charge spreads over).
pub const ORDER: usize = 4;

/// Cubic B-spline value `M4(w)` for `w` in `[0, 4)`.
///
/// # Panics
/// Panics (debug) if `w` is outside `[0, 4)`.
pub fn m4(w: f64) -> f64 {
    debug_assert!((0.0..4.0).contains(&w), "M4 argument {w} out of range");
    let w2 = w * w;
    let w3 = w2 * w;
    if w < 1.0 {
        w3 / 6.0
    } else if w < 2.0 {
        (-3.0 * w3 + 12.0 * w2 - 12.0 * w + 4.0) / 6.0
    } else if w < 3.0 {
        (3.0 * w3 - 24.0 * w2 + 60.0 * w - 44.0) / 6.0
    } else {
        let t = 4.0 - w;
        t * t * t / 6.0
    }
}

/// Derivative `M4'(w)` for `w` in `[0, 4)`.
///
/// # Panics
/// Panics (debug) if `w` is outside `[0, 4)`.
pub fn m4_derivative(w: f64) -> f64 {
    debug_assert!((0.0..4.0).contains(&w), "M4 argument {w} out of range");
    let w2 = w * w;
    if w < 1.0 {
        w2 / 2.0
    } else if w < 2.0 {
        (-9.0 * w2 + 24.0 * w - 12.0) / 6.0
    } else if w < 3.0 {
        (9.0 * w2 - 48.0 * w + 60.0) / 6.0
    } else {
        let t = 4.0 - w;
        -t * t / 2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn partition_of_unity() {
        // Σ_k M4(u − k) = 1 for any u.
        let mut u: f64 = 0.0;
        while u < 20.0 {
            let k0 = u.floor() as i64 - ORDER as i64;
            let sum: f64 = (k0..k0 + 2 * ORDER as i64)
                .map(|k| {
                    let w = u - k as f64;
                    if (0.0..4.0).contains(&w) {
                        m4(w)
                    } else {
                        0.0
                    }
                })
                .sum();
            assert!((sum - 1.0).abs() < 1e-12, "u = {u}: sum = {sum}");
            u += 0.037;
        }
    }

    #[test]
    fn derivative_matches_finite_difference() {
        let h = 1e-6;
        let mut w = 0.01;
        while w < 3.99 {
            let fd = (m4(w + h) - m4(w - h)) / (2.0 * h);
            let an = m4_derivative(w);
            assert!(
                (fd - an).abs() < 1e-7,
                "w = {w}: fd = {fd}, analytic = {an}"
            );
            w += 0.013;
        }
    }

    #[test]
    fn symmetry_about_center() {
        let mut w = 0.005;
        while w < 2.0 {
            assert!(
                (m4(w) - m4(4.0 - w)).abs() < 1e-14,
                "w = {w}: {m4:?} not symmetric",
                m4 = m4(w)
            );
            w += 0.043;
        }
    }

    #[test]
    fn integrates_to_one() {
        // ∫ M4 dw = 1 (numerically).
        let n = 400_000;
        let h = 4.0 / n as f64;
        let s: f64 = (0..n).map(|i| m4((i as f64 + 0.5) * h)).sum::<f64>() * h;
        assert!((s - 1.0).abs() < 1e-10, "integral = {s}");
    }
}
