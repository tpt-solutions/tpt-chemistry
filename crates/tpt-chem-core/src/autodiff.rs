//! Forward-mode automatic differentiation with dual numbers.
//!
//! Write an energy against the [`Real`] trait, then evaluate it with [`Dual`]
//! numbers to get exact derivatives, with no finite differences. Evaluating
//! with `f64` gives the plain energy, so the same code serves both. Forward
//! mode costs one evaluation per input, which suits molecular inputs (a few
//! coordinates per atom). It is not a substitute for reverse mode over large
//! parameter sets.
//!
//! ```
//! use tpt_chem_core::autodiff::{gradient, Real};
//!
//! // f(x, y) = x² y + sqrt(y), generic over the number type.
//! fn f<T: Real>(v: &[T]) -> T {
//!     v[0].powi(2) * v[1] + v[1].sqrt()
//! }
//! let (value, grad) = gradient(f, &[2.0, 4.0]);
//! assert!((value - (4.0 * 4.0 + 2.0)).abs() < 1e-12);
//! assert!((grad[0] - 2.0 * 2.0 * 4.0).abs() < 1e-12); // ∂f/∂x = 2xy
//! assert!((grad[1] - (4.0 + 0.25)).abs() < 1e-12); // ∂f/∂y = x² + 1/(2√y)
//! ```

use alloc::vec::Vec;
use core::ops::{Add, Div, Mul, Neg, Sub};

use crate::num;

/// A scalar that energies can be written against: plain `f64` or [`Dual`].
pub trait Real:
    Copy
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Neg<Output = Self>
{
    /// A constant with no derivative.
    fn constant(x: f64) -> Self;
    /// The underlying value.
    fn value(self) -> f64;
    /// Square root.
    fn sqrt(self) -> Self;
    /// Exponential.
    fn exp(self) -> Self;
    /// Arc cosine (for angles; the argument must lie strictly inside (−1, 1)).
    fn acos(self) -> Self;
    /// Integer power.
    fn powi(self, n: i32) -> Self;
}

impl Real for f64 {
    fn constant(x: f64) -> Self {
        x
    }
    fn value(self) -> f64 {
        self
    }
    fn sqrt(self) -> Self {
        num::sqrt(self)
    }
    fn exp(self) -> Self {
        num::exp(self)
    }
    fn acos(self) -> Self {
        num::acos(self)
    }
    fn powi(self, n: i32) -> Self {
        num::powi(self, n)
    }
}

/// A number carrying its value `re` and its derivative `eps` with respect to
/// the chosen input.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Dual {
    /// Value.
    pub re: f64,
    /// Derivative.
    pub eps: f64,
}

impl Dual {
    /// A value with an explicit derivative.
    pub const fn new(re: f64, eps: f64) -> Self {
        Dual { re, eps }
    }

    /// An input variable: derivative 1 with respect to itself.
    pub const fn variable(re: f64) -> Self {
        Dual { re, eps: 1.0 }
    }
}

impl Add for Dual {
    type Output = Dual;
    fn add(self, o: Dual) -> Dual {
        Dual::new(self.re + o.re, self.eps + o.eps)
    }
}

impl Sub for Dual {
    type Output = Dual;
    fn sub(self, o: Dual) -> Dual {
        Dual::new(self.re - o.re, self.eps - o.eps)
    }
}

impl Mul for Dual {
    type Output = Dual;
    fn mul(self, o: Dual) -> Dual {
        Dual::new(self.re * o.re, self.re * o.eps + self.eps * o.re)
    }
}

impl Div for Dual {
    type Output = Dual;
    fn div(self, o: Dual) -> Dual {
        Dual::new(
            self.re / o.re,
            (self.eps * o.re - self.re * o.eps) / (o.re * o.re),
        )
    }
}

impl Neg for Dual {
    type Output = Dual;
    fn neg(self) -> Dual {
        Dual::new(-self.re, -self.eps)
    }
}

impl Real for Dual {
    fn constant(x: f64) -> Self {
        Dual::new(x, 0.0)
    }
    fn value(self) -> f64 {
        self.re
    }
    fn sqrt(self) -> Self {
        let r = num::sqrt(self.re);
        Dual::new(r, self.eps / (2.0 * r))
    }
    fn exp(self) -> Self {
        let e = num::exp(self.re);
        Dual::new(e, self.eps * e)
    }
    fn acos(self) -> Self {
        let d = -1.0 / num::sqrt(1.0 - self.re * self.re);
        Dual::new(num::acos(self.re), self.eps * d)
    }
    fn powi(self, n: i32) -> Self {
        if n == 0 {
            return Dual::constant(1.0);
        }
        let p = num::powi(self.re, n);
        let dp = f64::from(n) * num::powi(self.re, n - 1);
        Dual::new(p, self.eps * dp)
    }
}

/// Value and gradient of `f` at `x`: one forward pass per input.
///
/// `f` must be generic over [`Real`], so it can take [`Dual`] inputs.
pub fn gradient<F>(f: F, x: &[f64]) -> (f64, Vec<f64>)
where
    F: Fn(&[Dual]) -> Dual,
{
    let value = if x.is_empty() {
        f(&[]).re
    } else {
        let vars: Vec<Dual> = x.iter().map(|&v| Dual::constant(v)).collect();
        f(&vars).re
    };
    let mut grad = Vec::with_capacity(x.len());
    for k in 0..x.len() {
        let vars: Vec<Dual> = x
            .iter()
            .enumerate()
            .map(|(i, &v)| {
                if i == k {
                    Dual::variable(v)
                } else {
                    Dual::constant(v)
                }
            })
            .collect();
        grad.push(f(&vars).eps);
    }
    (value, grad)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arithmetic_matches_calculus() {
        // x·(x+1)/x equals x + 1, so at x = 3 the value is 4 and the derivative is 1.
        let x = Dual::variable(3.0);
        let y = x * (x + Dual::constant(1.0)) / x;
        assert!((y.re - 4.0).abs() < 1e-12);
        assert!((y.eps - 1.0).abs() < 1e-12);
    }

    #[test]
    fn transcendental_derivatives() {
        let x = Dual::variable(0.7);
        let e = (x * x).exp();
        // d/dx e^{x²} = 2x e^{x²}
        assert!((e.eps - 2.0 * 0.7 * (0.49f64).exp()).abs() < 1e-12);
        let r = x.sqrt();
        assert!((r.eps - 0.5 / 0.7f64.sqrt()).abs() < 1e-12);
        let p = x.powi(-2);
        assert!((p.eps - (-2.0) * 0.7f64.powi(-3)).abs() < 1e-12);
    }

    #[test]
    fn gradient_of_a_generic_function() {
        fn f<T: Real>(v: &[T]) -> T {
            (v[0] * v[1]).exp() + v[1].powi(3)
        }
        let (value, grad) = gradient(f, &[0.5, 1.5]);
        let (x, y) = (0.5f64, 1.5f64);
        assert!((value - ((x * y).exp() + y.powi(3))).abs() < 1e-12);
        assert!((grad[0] - y * (x * y).exp()).abs() < 1e-12);
        assert!((grad[1] - (x * (x * y).exp() + 3.0 * y * y)).abs() < 1e-12);
    }

    #[test]
    fn plain_f64_evaluation_is_unchanged() {
        fn f<T: Real>(v: &[T]) -> T {
            v[0].sqrt() + v[1]
        }
        assert!((f(&[4.0f64, 1.0]) - 3.0).abs() < 1e-12);
    }
}
