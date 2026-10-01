//! Minimal dense linear algebra for the Roothaan–Hall equations: a
//! from-scratch cyclic Jacobi eigensolver for real symmetric matrices.
//!
//! spec.txt §4 routes HF through `tpt-math`'s eigensolvers, but
//! `tpt-math-linalg-dense` exposes no eigendecomposition and `nalgebra`
//! (the usual Rust backend) is Apache-2.0-only — banned by the workspace
//! license policy (spec.txt §2). The Jacobi method is the classic
//! from-scratch answer: unconditionally stable for symmetric matrices,
//! accurate, and perfectly adequate for the small dense matrices of
//! minimal-basis HF.

use tpt_chem_core::num;

/// Error type for the eigensolver.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LinError {
    /// The sweep limit was hit before the off-diagonal norm converged.
    NoConvergence,
    /// A dimension mismatch between operands.
    DimensionMismatch,
}

impl core::fmt::Display for LinError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            LinError::NoConvergence => write!(f, "Jacobi eigensolver did not converge"),
            LinError::DimensionMismatch => write!(f, "dimension mismatch"),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for LinError {}

/// Eigen-decompose a symmetric row-major `n × n` matrix.
///
/// Returns `(eigenvalues ascending, eigenvectors as columns in a row-major
/// `n × n` matrix)` such that `A = V diag(w) Vᵀ`.
///
/// # Errors
/// [`LinError::NoConvergence`] if not converged within `max_sweeps`
/// full sweeps (50 is far beyond what symmetric matrices need).
pub fn jacobi_eigh(
    a: &[f64],
    n: usize,
    tol: f64,
    max_sweeps: usize,
) -> Result<(Vec<f64>, Vec<f64>), LinError> {
    if a.len() != n * n {
        return Err(LinError::DimensionMismatch);
    }
    let mut m = a.to_vec();
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    for _sweep in 0..max_sweeps {
        // Off-diagonal Frobenius norm².
        let mut off = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                off += m[i * n + j] * m[i * n + j];
            }
        }
        if off.sqrt() < tol {
            break;
        }
        for p in 0..n {
            for q in (p + 1)..n {
                let apq = m[p * n + q];
                if apq.abs() < tol * 1e-3 {
                    continue;
                }
                // Jacobi rotation angles (numerically stable form).
                let theta = (m[q * n + q] - m[p * n + p]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + num::sqrt(theta * theta + 1.0));
                let c = 1.0 / num::sqrt(t * t + 1.0);
                let s = t * c;
                // Apply rotation to rows/columns p, q of m and columns of v.
                for k in 0..n {
                    let mkp = m[k * n + p];
                    let mkq = m[k * n + q];
                    m[k * n + p] = c * mkp - s * mkq;
                    m[k * n + q] = s * mkp + c * mkq;
                }
                for k in 0..n {
                    let mpk = m[p * n + k];
                    let mqk = m[q * n + k];
                    m[p * n + k] = c * mpk - s * mqk;
                    m[q * n + k] = s * mpk + c * mqk;
                }
                for k in 0..n {
                    let vkp = v[k * n + p];
                    let vkq = v[k * n + q];
                    v[k * n + p] = c * vkp - s * vkq;
                    v[k * n + q] = s * vkp + c * vkq;
                }
            }
        }
    }
    // Final convergence check.
    let mut off = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            off += m[i * n + j] * m[i * n + j];
        }
    }
    if off.sqrt() > tol {
        return Err(LinError::NoConvergence);
    }
    // Extract ascending eigenvalues, permuting eigenvector columns.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| {
        m[i * n + i]
            .partial_cmp(&m[j * n + j])
            .unwrap_or(core::cmp::Ordering::Equal)
    });
    let w: Vec<f64> = order.iter().map(|&i| m[i * n + i]).collect();
    let mut vecs = vec![0.0; n * n];
    for (new_col, &old_col) in order.iter().enumerate() {
        for r in 0..n {
            vecs[r * n + new_col] = v[r * n + old_col];
        }
    }
    Ok((w, vecs))
}

/// Matrix product of row-major `a (n×n)` and `b (n×n)`.
pub fn mat_mul(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut out = vec![0.0; n * n];
    for i in 0..n {
        for k in 0..n {
            let aik = a[i * n + k];
            if aik == 0.0 {
                continue;
            }
            for j in 0..n {
                out[i * n + j] += aik * b[k * n + j];
            }
        }
    }
    out
}

/// Transpose of a row-major `n × n` matrix.
pub fn transpose(a: &[f64], n: usize) -> Vec<f64> {
    let mut out = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            out[j * n + i] = a[i * n + j];
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn apply(a: &[f64], n: usize, x: &[f64]) -> Vec<f64> {
        (0..n)
            .map(|i| (0..n).map(|j| a[i * n + j] * x[j]).sum())
            .collect()
    }

    #[test]
    fn diagonal_matrix_passes_through() {
        let a = vec![3.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 2.0];
        let (w, v) = jacobi_eigh(&a, 3, 1e-12, 50).unwrap();
        assert_eq!(w, vec![1.0, 2.0, 3.0]);
        for (k, &ev) in w.iter().enumerate() {
            let x: Vec<f64> = (0..3).map(|r| v[r * 3 + k]).collect();
            let ax = apply(&a, 3, &x);
            for r in 0..3 {
                assert!((ax[r] - ev * x[r]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn general_symmetric_3x3() {
        let a = vec![2.0, -1.0, 0.0, -1.0, 2.0, -1.0, 0.0, -1.0, 2.0];
        let (w, v) = jacobi_eigh(&a, 3, 1e-12, 50).unwrap();
        // Analytic eigenvalues of the tridiagonal Toeplitz: 2 − √2, 2, 2 + √2.
        assert!((w[0] - (2.0 - core::f64::consts::SQRT_2)).abs() < 1e-10);
        assert!((w[1] - 2.0).abs() < 1e-10);
        assert!((w[2] - (2.0 + core::f64::consts::SQRT_2)).abs() < 1e-10);
        for (k, &ev) in w.iter().enumerate() {
            let x: Vec<f64> = (0..3).map(|r| v[r * 3 + k]).collect();
            let ax = apply(&a, 3, &x);
            for r in 0..3 {
                assert!((ax[r] - ev * x[r]).abs() < 1e-9);
            }
        }
    }

    #[test]
    fn reconstruction_identity() {
        // A = V diag(w) Vᵀ.
        let a = vec![4.0, 1.0, 0.5, 1.0, 3.0, 0.2, 0.5, 0.2, 1.0];
        let n = 3;
        let (w, v) = jacobi_eigh(&a, n, 1e-13, 50).unwrap();
        let mut diag = vec![0.0; n * n];
        for i in 0..n {
            diag[i * n + i] = w[i];
        }
        let vt = transpose(&v, n);
        let recomposed = mat_mul(&mat_mul(&v, &diag, n), &vt, n);
        for i in 0..n * n {
            assert!((recomposed[i] - a[i]).abs() < 1e-10);
        }
    }

    #[test]
    fn dimension_mismatch() {
        assert_eq!(
            jacobi_eigh(&[1.0, 2.0], 2, 1e-12, 5),
            Err(LinError::DimensionMismatch)
        );
    }
}
