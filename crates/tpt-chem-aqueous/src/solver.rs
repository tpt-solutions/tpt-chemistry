//! Small scalar root finder (Brent's method) used by the speciation and
//! onset-pH searches.

use crate::{AqueousError, Result};

/// Finds a root of `f` in the bracket `[a, b]` (Brent's method).
///
/// `f(a)` and `f(b)` must have opposite signs (or one must be exactly zero).
pub(crate) fn brent<F: FnMut(f64) -> f64>(
    mut f: F,
    mut a: f64,
    mut b: f64,
    tol: f64,
    max_iter: usize,
) -> Result<f64> {
    let mut fa = f(a);
    let mut fb = f(b);
    if fa == 0.0 {
        return Ok(a);
    }
    if fb == 0.0 {
        return Ok(b);
    }
    if !(fa.is_finite() && fb.is_finite()) {
        return Err(AqueousError::NonFinite("root-finder endpoint"));
    }
    if fa * fb > 0.0 {
        return Err(AqueousError::NotBracketed);
    }
    let (mut c, mut fc) = (a, fa);
    let mut d = b - a;
    let mut e = d;
    for _ in 0..max_iter {
        if fb * fc > 0.0 {
            c = a;
            fc = fa;
            d = b - a;
            e = d;
        }
        if fc.abs() < fb.abs() {
            a = b;
            b = c;
            c = a;
            fa = fb;
            fb = fc;
            fc = fa;
        }
        let tol1 = 2.0 * f64::EPSILON * b.abs() + 0.5 * tol;
        let xm = 0.5 * (c - b);
        if xm.abs() <= tol1 || fb == 0.0 {
            return Ok(b);
        }
        if e.abs() >= tol1 && fa.abs() > fb.abs() {
            let s = fb / fa;
            let (mut p, mut q);
            if a == c {
                p = 2.0 * xm * s;
                q = 1.0 - s;
            } else {
                let qq = fa / fc;
                let r = fb / fc;
                p = s * (2.0 * xm * qq * (qq - r) - (b - a) * (r - 1.0));
                q = (qq - 1.0) * (r - 1.0) * (s - 1.0);
            }
            if p > 0.0 {
                q = -q;
            }
            p = p.abs();
            if 2.0 * p < (3.0 * xm * q - (tol1 * q).abs()).min((e * q).abs()) {
                e = d;
                d = p / q;
            } else {
                d = xm;
                e = d;
            }
        } else {
            d = xm;
            e = d;
        }
        a = b;
        fa = fb;
        b += if d.abs() > tol1 { d } else { tol1.copysign(xm) };
        fb = f(b);
        if !fb.is_finite() {
            return Err(AqueousError::NonFinite("root-finder iterate"));
        }
    }
    Err(AqueousError::NoConvergence("brent"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_sqrt2() {
        let r = brent(|x| x * x - 2.0, 0.0, 2.0, 1e-14, 100).unwrap();
        assert!((r - 2.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn rejects_unbracketed() {
        assert_eq!(
            brent(|x| x * x + 1.0, -1.0, 1.0, 1e-12, 50),
            Err(AqueousError::NotBracketed)
        );
    }
}
