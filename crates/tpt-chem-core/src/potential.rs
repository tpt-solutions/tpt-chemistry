//! A common energy/force interface shared by classical force fields,
//! quantum-chemistry back-ends and (future) machine-learned potentials,
//! plus a generic BFGS minimiser that works with any of them.
//!
//! Units are the MD units used across the workspace: positions in Å,
//! energies in kJ·mol⁻¹, forces in kJ·mol⁻¹·Å⁻¹ (force = −∇E).

use alloc::vec;
use alloc::vec::Vec;

use crate::vec3::Vec3;

/// Anything that maps atomic positions to an energy and forces.
pub trait Potential {
    /// Evaluation failure (e.g. SCF non-convergence).
    type Error;

    /// Number of atoms the potential expects.
    fn n_atoms(&self) -> usize;

    /// Energy (kJ·mol⁻¹) and forces (kJ·mol⁻¹·Å⁻¹) at `positions` (Å).
    ///
    /// # Errors
    /// Implementation-specific.
    fn energy_forces(&mut self, positions: &[Vec3]) -> Result<(f64, Vec<Vec3>), Self::Error>;
}

/// Settings for [`minimize_bfgs`].
#[derive(Clone, Copy, Debug)]
pub struct MinimizeSettings {
    /// Maximum accepted steps.
    pub max_steps: usize,
    /// Convergence threshold on the largest force component
    /// (kJ·mol⁻¹·Å⁻¹).
    pub force_tol: f64,
    /// Largest allowed coordinate displacement per step (Å).
    pub max_step: f64,
}

impl Default for MinimizeSettings {
    fn default() -> Self {
        MinimizeSettings {
            max_steps: 200,
            force_tol: 1.0,
            max_step: 0.2,
        }
    }
}

/// Outcome of [`minimize_bfgs`].
#[derive(Clone, Copy, Debug)]
pub struct MinimizeReport {
    /// Accepted steps.
    pub steps: usize,
    /// Final energy (kJ·mol⁻¹).
    pub energy: f64,
    /// Final largest force component (kJ·mol⁻¹·Å⁻¹).
    pub max_force: f64,
    /// Whether `force_tol` was reached.
    pub converged: bool,
}

fn abs(x: f64) -> f64 {
    if x < 0.0 {
        -x
    } else {
        x
    }
}

fn flat(v: &[Vec3]) -> Vec<f64> {
    v.iter().flat_map(|p| [p.x, p.y, p.z]).collect()
}

fn max_abs(v: &[f64]) -> f64 {
    v.iter()
        .fold(0.0, |m, &x| if abs(x) > m { abs(x) } else { m })
}

fn unflat(x: &[f64]) -> Vec<Vec3> {
    x.chunks(3).map(|c| Vec3::new(c[0], c[1], c[2])).collect()
}

/// Minimise a [`Potential`] with BFGS (dense inverse Hessian, backtracking
/// line search, per-step displacement cap). `positions` is updated in
/// place.
///
/// # Errors
/// Propagates the potential's own error.
pub fn minimize_bfgs<P: Potential>(
    pot: &mut P,
    positions: &mut [Vec3],
    settings: MinimizeSettings,
) -> Result<MinimizeReport, P::Error> {
    let n = 3 * positions.len();
    let ident = || {
        let mut h = vec![0.0; n * n];
        for i in 0..n {
            h[i * n + i] = 1.0;
        }
        h
    };
    let mut x = flat(positions);
    let (mut e, f0) = pot.energy_forces(positions)?;
    // Gradient = -force.
    let mut g: Vec<f64> = flat(&f0).iter().map(|v| -v).collect();
    // Initial inverse Hessian scale: a stiff-ish 1/(500 kJ/mol/A^2).
    let h0 = 1.0 / 500.0;
    let scaled_ident = || {
        let mut h = ident();
        for v in h.iter_mut() {
            *v *= h0;
        }
        h
    };
    let mut hinv = scaled_ident();
    let mut steps = 0;
    while max_abs(&g) > settings.force_tol && steps < settings.max_steps {
        let mut dir: Vec<f64> = (0..n)
            .map(|i| -(0..n).map(|j| hinv[i * n + j] * g[j]).sum::<f64>())
            .collect();
        if dir.iter().zip(&g).map(|(d, gg)| d * gg).sum::<f64>() >= 0.0 {
            hinv = scaled_ident();
            dir = g.iter().map(|v| -v * h0).collect();
        }
        let dmax = max_abs(&dir);
        let mut scale = if dmax > settings.max_step {
            settings.max_step / dmax
        } else {
            1.0
        };
        let mut accepted = None;
        for _ in 0..12 {
            let xn: Vec<f64> = x.iter().zip(&dir).map(|(a, d)| a + scale * d).collect();
            let (en, fnew) = pot.energy_forces(&unflat(&xn))?;
            if en < e + 1e-12 {
                accepted = Some((xn, en, fnew));
                break;
            }
            scale *= 0.5;
        }
        let Some((xn, en, fnew)) = accepted else {
            break;
        };
        let gn: Vec<f64> = flat(&fnew).iter().map(|v| -v).collect();
        let s: Vec<f64> = xn.iter().zip(&x).map(|(a, b)| a - b).collect();
        let y: Vec<f64> = gn.iter().zip(&g).map(|(a, b)| a - b).collect();
        let sy: f64 = s.iter().zip(&y).map(|(a, b)| a * b).sum();
        if sy > 1e-12 {
            let hy: Vec<f64> = (0..n)
                .map(|i| (0..n).map(|j| hinv[i * n + j] * y[j]).sum())
                .collect();
            let yhy: f64 = y.iter().zip(&hy).map(|(a, b)| a * b).sum();
            for i in 0..n {
                for j in 0..n {
                    hinv[i * n + j] +=
                        (1.0 + yhy / sy) * s[i] * s[j] / sy - (hy[i] * s[j] + s[i] * hy[j]) / sy;
                }
            }
        }
        x = xn;
        e = en;
        g = gn;
        steps += 1;
    }
    for (p, c) in positions.iter_mut().zip(x.chunks(3)) {
        *p = Vec3::new(c[0], c[1], c[2]);
    }
    let max_force = max_abs(&g);
    Ok(MinimizeReport {
        steps,
        energy: e,
        max_force,
        converged: max_force <= settings.force_tol,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::convert::Infallible;

    /// Two atoms joined by a harmonic spring (k, r0) plus a weak z-tether.
    struct Spring {
        k: f64,
        r0: f64,
    }

    impl Potential for Spring {
        type Error = Infallible;
        fn n_atoms(&self) -> usize {
            2
        }
        fn energy_forces(&mut self, p: &[Vec3]) -> Result<(f64, Vec<Vec3>), Infallible> {
            let d = p[1] - p[0];
            let r = d.norm();
            let e = 0.5 * self.k * (r - self.r0) * (r - self.r0);
            let f1 = d * (-self.k * (r - self.r0) / r);
            Ok((e, vec![-f1, f1]))
        }
    }

    #[test]
    fn bfgs_relaxes_a_spring() {
        let mut pot = Spring { k: 400.0, r0: 1.2 };
        let mut pos = [Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.3, 0.9, 1.7)];
        let rep = minimize_bfgs(&mut pot, &mut pos, MinimizeSettings::default()).unwrap();
        assert!(rep.converged, "{rep:?}");
        assert!(((pos[1] - pos[0]).norm() - 1.2).abs() < 1e-3);
        assert!(rep.energy < 1e-4);
    }
}
