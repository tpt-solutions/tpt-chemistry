//! Special functions: the error function family, needed by Ewald
//! electrostatics (`erfc` screening) and statistics.
//!
//! Implemented from scratch (spec.txt §2): power series for small arguments,
//! the classical continued fraction (Abramowitz & Stegun 7.1.14) for large.
//! Accuracy is ~1e-14 relative, far below the tolerances used by the Ewald
//! convergence tests in `tpt-chem-md`.

use crate::num;

/// Gauss' error function `erf(x)`.
pub fn erf(x: f64) -> f64 {
    if x < 0.0 {
        -erf(-x)
    } else if x < 1.5 {
        erf_series(x)
    } else {
        1.0 - erfc_cf(x)
    }
}

/// Complementary error function `erfc(x) = 1 - erf(x)`.
pub fn erfc(x: f64) -> f64 {
    if x < 0.0 {
        2.0 - erfc(-x)
    } else if x < 1.5 {
        1.0 - erf_series(x)
    } else {
        erfc_cf(x)
    }
}

/// `erf` by the Maclaurin series, accurate for |x| < ~1.5.
fn erf_series(x: f64) -> f64 {
    let two_over_sqrt_pi = 2.0 / crate::num::sqrt(core::f64::consts::PI);
    let x2 = x * x;
    let mut term = x; // x^(2n+1)/n! before division by (2n+1)
    let mut sum = x;
    let mut n: u32 = 1;
    while n < 200 {
        term *= -x2 / f64::from(n);
        let contrib = term / f64::from(2 * n + 1);
        sum += contrib;
        if contrib.abs() <= 1e-17 * sum.abs() {
            break;
        }
        n += 1;
    }
    two_over_sqrt_pi * sum
}

/// `erfc` by the continued fraction
/// `√π·e^{x²}·erfc(x) = 1/(x + (1/2)/(x + (2/2)/(x + (3/2)/(x + …))))`,
/// evaluated backwards. Accurate for x ≥ ~1.5 where the series loses
/// precision to cancellation.
fn erfc_cf(x: f64) -> f64 {
    const TERMS: usize = 300;
    let mut t = x;
    for k in (1..=TERMS).rev() {
        t = x + (f64::from(k as u32) * 0.5) / t;
    }
    num::exp(-x * x) / (crate::num::sqrt(core::f64::consts::PI) * t)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn erf_reference_values() {
        assert!(close(erf(0.0), 0.0, 1e-15));
        assert!(close(erf(0.5), 0.520_499_877_813_046_5, 1e-13));
        assert!(close(erf(1.0), 0.842_700_792_949_714_9, 1e-13));
        assert!(close(erf(-1.0), -0.842_700_792_949_714_9, 1e-13));
        assert!(close(erf(3.0), 0.999_977_909_503_001_4, 1e-13));
    }

    #[test]
    fn erfc_reference_values() {
        assert!(close(erfc(0.0), 1.0, 1e-15));
        assert!(close(erfc(0.5), 0.479_500_122_186_953_5, 1e-13));
        assert!(close(erfc(1.0), 0.157_299_207_050_285_13, 1e-13));
        assert!(close(erfc(2.0), 4.677_734_981_047_266e-3, 1e-16));
        assert!(close(erfc(3.0), 2.209_049_699_858_544e-5, 1e-19));
        assert!(close(erfc(-1.0), 1.842_700_792_949_714_8, 1e-13));
        // Deep tail: relative accuracy.
        assert!(close(erfc(6.0), 2.151_973_671_249_891_3e-17, 1e-31));
    }

    #[test]
    fn erfc_erf_consistent() {
        for i in -20..=20 {
            let x = f64::from(i) * 0.25;
            assert!((erf(x) + erfc(x) - 1.0).abs() < 1e-13, "x = {x}");
        }
    }
}
