//! Reverse-mode automatic differentiation (tape-based).
//!
//! The same generic code as [`crate::autodiff`] works here: write the energy
//! against [`Real`](crate::autodiff::Real), and [`reverse_gradient`] records
//! the operations on a tape during one forward pass, then propagates adjoints
//! backwards. The gradient with respect to *all* inputs costs one forward and
//! one backward sweep, independent of the number of inputs. That is what makes
//! it suit large parameter sets, where forward mode would need one pass per
//! parameter.
//!
//! The tape is thread-local, so this module needs the `std` feature. Forward
//! mode in [`crate::autodiff`] stays `no_std`.
//!
//! ```
//! use tpt_chem_core::autodiff::Real;
//! use tpt_chem_core::autodiff_reverse::reverse_gradient;
//!
//! // f(x, y) = x² y + sqrt(y), the same function as the forward-mode example.
//! fn f<T: Real>(v: &[T]) -> T {
//!     v[0].powi(2) * v[1] + v[1].sqrt()
//! }
//! let (value, grad) = reverse_gradient(f, &[2.0, 4.0]);
//! assert!((value - 18.0).abs() < 1e-12);
//! assert!((grad[0] - 16.0).abs() < 1e-12); // ∂f/∂x = 2xy
//! assert!((grad[1] - 4.25).abs() < 1e-12); // ∂f/∂y = x² + 1/(2√y)
//! ```
//!
//! Calls must not nest: a function passed to [`reverse_gradient`] must not
//! call it again, because both use the same tape.

use std::cell::RefCell;
use std::vec::Vec;

use crate::autodiff::Real;

/// One recorded operation: its value, up to two parents on the tape and the
/// local partial derivative with respect to each parent.
#[derive(Clone, Copy, Debug)]
struct Node {
    value: f64,
    parents: [usize; 2],
    partials: [f64; 2],
    arity: u8,
}

thread_local! {
    static TAPE: RefCell<Vec<Node>> = const { RefCell::new(Vec::new()) };
}

/// A value recorded on the tape. It is an index, so it is `Copy` like `f64`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Var(usize);

fn push(node: Node) -> Var {
    TAPE.with(|t| {
        let mut t = t.borrow_mut();
        t.push(node);
        Var(t.len() - 1)
    })
}

fn unary(a: Var, value: f64, partial: f64) -> Var {
    push(Node {
        value,
        parents: [a.0, 0],
        partials: [partial, 0.0],
        arity: 1,
    })
}

fn binary(a: Var, b: Var, value: f64, partials: [f64; 2]) -> Var {
    push(Node {
        value,
        parents: [a.0, b.0],
        partials,
        arity: 2,
    })
}

fn node_value(v: Var) -> f64 {
    TAPE.with(|t| t.borrow()[v.0].value)
}

impl core::ops::Add for Var {
    type Output = Var;
    fn add(self, o: Var) -> Var {
        binary(self, o, node_value(self) + node_value(o), [1.0, 1.0])
    }
}

impl core::ops::Sub for Var {
    type Output = Var;
    fn sub(self, o: Var) -> Var {
        binary(self, o, node_value(self) - node_value(o), [1.0, -1.0])
    }
}

impl core::ops::Mul for Var {
    type Output = Var;
    fn mul(self, o: Var) -> Var {
        let (a, b) = (node_value(self), node_value(o));
        binary(self, o, a * b, [b, a])
    }
}

impl core::ops::Div for Var {
    type Output = Var;
    fn div(self, o: Var) -> Var {
        let (a, b) = (node_value(self), node_value(o));
        binary(self, o, a / b, [1.0 / b, -a / (b * b)])
    }
}

impl core::ops::Neg for Var {
    type Output = Var;
    fn neg(self) -> Var {
        unary(self, -node_value(self), -1.0)
    }
}

impl Real for Var {
    fn constant(x: f64) -> Self {
        push(Node {
            value: x,
            parents: [0, 0],
            partials: [0.0, 0.0],
            arity: 0,
        })
    }
    fn value(self) -> f64 {
        node_value(self)
    }
    fn sqrt(self) -> Self {
        let r = node_value(self).sqrt();
        unary(self, r, 0.5 / r)
    }
    fn exp(self) -> Self {
        let e = node_value(self).exp();
        unary(self, e, e)
    }
    fn acos(self) -> Self {
        let x = node_value(self);
        unary(self, x.acos(), -1.0 / (1.0 - x * x).sqrt())
    }
    fn powi(self, n: i32) -> Self {
        if n == 0 {
            return Var::constant(1.0);
        }
        let x = node_value(self);
        let p = x.powi(n);
        unary(self, p, f64::from(n) * x.powi(n - 1))
    }
}

/// Value and gradient of `f` at `x` by reverse-mode differentiation.
///
/// `f` must be generic over [`Real`]; it is called once, with [`Var`] inputs.
pub fn reverse_gradient<F>(f: F, x: &[f64]) -> (f64, Vec<f64>)
where
    F: Fn(&[Var]) -> Var,
{
    TAPE.with(|t| t.borrow_mut().clear());
    let inputs: Vec<Var> = x.iter().map(|&v| Var::constant(v)).collect();
    let out = f(&inputs);
    let value = node_value(out);

    let grad = TAPE.with(|t| {
        let tape = t.borrow();
        let mut adjoint = vec![0.0; tape.len()];
        adjoint[out.0] = 1.0;
        // Parents always precede their children on the tape, so one reverse
        // sweep from the output is enough.
        for i in (0..=out.0).rev() {
            let node = tape[i];
            let a = adjoint[i];
            for k in 0..usize::from(node.arity) {
                adjoint[node.parents[k]] += a * node.partials[k];
            }
        }
        inputs.iter().map(|v| adjoint[v.0]).collect()
    });
    TAPE.with(|t| t.borrow_mut().clear());
    (value, grad)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::autodiff::{gradient, Dual};

    /// `rev` and `fwd` are the same generic function, instantiated for
    /// [`Var`] and [`Dual`](Dual).
    fn assert_matches_forward(
        rev: impl Fn(&[Var]) -> Var,
        fwd: impl Fn(&[Dual]) -> Dual,
        x: &[f64],
    ) {
        let (v_rev, g_rev) = reverse_gradient(rev, x);
        let (v_fwd, g_fwd) = gradient(fwd, x);
        assert!((v_rev - v_fwd).abs() < 1e-12);
        for (a, b) in g_rev.iter().zip(&g_fwd) {
            assert!((a - b).abs() < 1e-10 * (1.0 + b.abs()), "{a} vs {b}");
        }
    }

    #[test]
    fn every_operation_matches_forward_mode() {
        fn f<T: Real>(v: &[T]) -> T {
            let (x, y) = (v[0], v[1]);
            let a = (x * y - x / y + (-y)) * x.powi(3);
            let b = (x * x).exp() + (y + T::constant(2.0)).sqrt();
            let c = (T::constant(0.5) * x).acos() * y.powi(-2);
            a + b - c
        }
        assert_matches_forward(f, f, &[0.7, 1.3]);
    }

    #[test]
    fn shared_subexpressions_accumulate_adjoints() {
        // f = s + s·s with s = x·y: ∂f/∂x = y (1 + 2s) at each node use.
        fn f<T: Real>(v: &[T]) -> T {
            let s = v[0] * v[1];
            s + s * s
        }
        let (x, y) = (0.6, 1.1);
        let (_, g) = reverse_gradient(f, &[x, y]);
        let s = x * y;
        assert!((g[0] - y * (1.0 + 2.0 * s)).abs() < 1e-12);
        assert!((g[1] - x * (1.0 + 2.0 * s)).abs() < 1e-12);
    }

    #[test]
    fn many_inputs_in_one_pass() {
        // Σ_i (x_i · x_{i+1})^2: each gradient entry touches two neighbours.
        fn f<T: Real>(v: &[T]) -> T {
            let mut total = T::constant(0.0);
            for w in v.windows(2) {
                let p = w[0] * w[1];
                total = total + p * p;
            }
            total
        }
        let x: Vec<f64> = (0..200).map(|i| 0.01 * f64::from(i) - 1.0).collect();
        assert_matches_forward(f, f, &x);
    }

    #[test]
    fn tape_is_cleared_between_calls() {
        fn f<T: Real>(v: &[T]) -> T {
            v[0] * v[0]
        }
        let _ = reverse_gradient(f, &[1.0]);
        let _ = reverse_gradient(f, &[2.0]);
        TAPE.with(|t| assert!(t.borrow().is_empty()));
    }

    #[test]
    fn constant_only_function_has_zero_gradient() {
        fn f<T: Real>(_: &[T]) -> T {
            T::constant(3.0)
        }
        let (v, g) = reverse_gradient(f, &[1.0, 2.0]);
        assert_eq!(v, 3.0);
        assert_eq!(g, vec![0.0, 0.0]);
    }
}
