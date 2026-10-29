//! Machine-learned interatomic potential: a radial Behler–Parrinello-style
//! model with a pure-Rust multilayer perceptron.
//!
//! Each atom gets a vector of radial symmetry functions
//! `G_k(i) = Σ_j exp(−η_k (r_ij − R_k)²) · f_c(r_ij)` with a cosine cutoff
//! `f_c`. A shared tanh network maps each atom's `G(i)` to an atomic energy,
//! and the total energy is the sum over atoms.
//!
//! Forces come from analytic backpropagation through the network and the
//! symmetry functions, so they are the exact gradient of the energy. Under
//! [`VelocityVerlet::step_with`](crate::integrator::VelocityVerlet::step_with)
//! that makes energy conserved to integrator accuracy, which the tests check.
//! [`RadialMlp::train_energies`] fits the network to labelled total energies
//! by full-batch gradient descent.
//!
//! Units follow the crate convention: Å for distances, kJ·mol⁻¹ for energies.
//! The network is untrained when constructed. Fit it before use.

use tpt_chem_core::rng::Rng;
use tpt_chem_core::vec3::Vec3;

use crate::box3::Box3;
use crate::system::System;

/// Radial Gaussian symmetry functions with a cosine cutoff.
#[derive(Clone, Debug, PartialEq)]
pub struct RadialSymmetry {
    /// Cutoff radius (Å); pairs beyond it contribute nothing.
    pub r_cut: f64,
    /// Gaussian centres `R_k` (Å).
    pub centers: Vec<f64>,
    /// Gaussian widths `η_k` (Å⁻²).
    pub etas: Vec<f64>,
}

impl RadialSymmetry {
    /// `n` Gaussians with centres spaced evenly from `r_min` to `r_cut`, all
    /// with width `eta`.
    ///
    /// # Panics
    /// Panics when `n` is zero.
    pub fn uniform(r_min: f64, r_cut: f64, n: usize, eta: f64) -> Self {
        assert!(n > 0, "need at least one symmetry function");
        let centers = (0..n)
            .map(|k| {
                if n == 1 {
                    r_min
                } else {
                    r_min + (r_cut - r_min) * k as f64 / (n - 1) as f64
                }
            })
            .collect();
        RadialSymmetry {
            r_cut,
            centers,
            etas: vec![eta; n],
        }
    }

    /// Number of symmetry functions.
    pub fn len(&self) -> usize {
        self.centers.len()
    }

    /// Whether there are no symmetry functions.
    pub fn is_empty(&self) -> bool {
        self.centers.is_empty()
    }

    /// Values `g_k(r)` and derivatives `dg_k/dr` at distance `r`.
    fn radial(&self, r: f64) -> (Vec<f64>, Vec<f64>) {
        let k = self.len();
        if r >= self.r_cut {
            return (vec![0.0; k], vec![0.0; k]);
        }
        let phase = std::f64::consts::PI * r / self.r_cut;
        let fc = 0.5 * (phase.cos() + 1.0);
        let dfc = -0.5 * std::f64::consts::PI / self.r_cut * phase.sin();
        let mut g = Vec::with_capacity(k);
        let mut dg = Vec::with_capacity(k);
        for (&rs, &eta) in self.centers.iter().zip(&self.etas) {
            let d = r - rs;
            let gauss = (-eta * d * d).exp();
            g.push(gauss * fc);
            dg.push(gauss * (-2.0 * eta * d * fc + dfc));
        }
        (g, dg)
    }
}

/// One dense layer: `out = W·in + b`, with `W` stored row-major (`n_out × n_in`).
#[derive(Clone, Debug)]
struct Dense {
    n_in: usize,
    w: Vec<f64>,
    b: Vec<f64>,
}

/// Parameter gradients, one `(dW, db)` pair per layer.
type Grads = Vec<(Vec<f64>, Vec<f64>)>;

/// A labelled configuration for training: positions (Å), optional periodic
/// box, and the reference total energy (kJ·mol⁻¹).
#[derive(Clone, Debug)]
pub struct LabelledConfig {
    /// Atomic positions (Å).
    pub pos: Vec<Vec3>,
    /// Periodic box, or `None` for an isolated cluster.
    pub box_: Option<Box3>,
    /// Reference total energy (kJ·mol⁻¹).
    pub energy: f64,
}

/// Radial symmetry functions followed by a tanh multilayer perceptron.
#[derive(Clone, Debug)]
pub struct RadialMlp {
    sym: RadialSymmetry,
    layers: Vec<Dense>,
}

/// One interacting pair with the geometry needed for the force chain rule.
struct PairTerm {
    i: usize,
    j: usize,
    /// Unit vector from `j` to `i`.
    unit: Vec3,
    dg: Vec<f64>,
}

impl RadialMlp {
    /// Network with `hidden` tanh layer widths, Xavier-uniform weights from
    /// the seeded generator `seed`, and a scalar linear output.
    pub fn new(sym: RadialSymmetry, hidden: &[usize], seed: u64) -> Self {
        let mut sizes = vec![sym.len()];
        sizes.extend_from_slice(hidden);
        sizes.push(1);
        let mut rng = Rng::new(seed);
        let layers = sizes
            .windows(2)
            .map(|w| {
                let (n_in, n_out) = (w[0], w[1]);
                let scale = (6.0 / (n_in + n_out) as f64).sqrt();
                Dense {
                    n_in,
                    w: (0..n_in * n_out)
                        .map(|_| rng.range(-scale, scale))
                        .collect(),
                    b: vec![0.0; n_out],
                }
            })
            .collect();
        RadialMlp { sym, layers }
    }

    /// The symmetry-function set.
    pub fn symmetry(&self) -> &RadialSymmetry {
        &self.sym
    }

    /// Forward pass; returns every layer's activations, starting with the input.
    fn forward(&self, x: &[f64]) -> Vec<Vec<f64>> {
        let last = self.layers.len() - 1;
        let mut acts = vec![x.to_vec()];
        for (l, layer) in self.layers.iter().enumerate() {
            let input = &acts[l];
            let mut out = layer.b.clone();
            for (o, value) in out.iter_mut().enumerate() {
                let row = &layer.w[o * layer.n_in..(o + 1) * layer.n_in];
                *value += row.iter().zip(input).map(|(w, a)| w * a).sum::<f64>();
                if l < last {
                    *value = value.tanh();
                }
            }
            acts.push(out);
        }
        acts
    }

    /// Backward pass for upstream gradient `upstream` on the scalar output.
    ///
    /// Returns `d(output)/d(input)`. When `grads` is given, the parameter
    /// gradients are added into it.
    fn backward(
        &self,
        acts: &[Vec<f64>],
        upstream: f64,
        mut grads: Option<&mut [(Vec<f64>, Vec<f64>)]>,
    ) -> Vec<f64> {
        let last = self.layers.len() - 1;
        let mut delta = vec![upstream];
        for (l, layer) in self.layers.iter().enumerate().rev() {
            let out = &acts[l + 1];
            // Through the tanh (hidden layers only; the output is linear).
            let dz: Vec<f64> = if l < last {
                delta
                    .iter()
                    .zip(out)
                    .map(|(d, a)| d * (1.0 - a * a))
                    .collect()
            } else {
                delta.clone()
            };
            if let Some(g) = grads.as_deref_mut() {
                let (dw, db) = &mut g[l];
                for (o, &dzo) in dz.iter().enumerate() {
                    db[o] += dzo;
                    for (i, &a) in acts[l].iter().enumerate() {
                        dw[o * layer.n_in + i] += dzo * a;
                    }
                }
            }
            let mut prev = vec![0.0; layer.n_in];
            for (o, &dzo) in dz.iter().enumerate() {
                let row = &layer.w[o * layer.n_in..(o + 1) * layer.n_in];
                for (p, w) in prev.iter_mut().zip(row) {
                    *p += w * dzo;
                }
            }
            delta = prev;
        }
        delta
    }

    fn zero_grads(&self) -> Grads {
        self.layers
            .iter()
            .map(|l| (vec![0.0; l.w.len()], vec![0.0; l.b.len()]))
            .collect()
    }

    /// Symmetry vectors per atom and the list of interacting pairs.
    fn descriptors(&self, pos: &[Vec3], box_: Option<Box3>) -> (Vec<Vec<f64>>, Vec<PairTerm>) {
        let n = pos.len();
        let k = self.sym.len();
        let mut g = vec![vec![0.0; k]; n];
        let mut pairs = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                let raw = pos[i] - pos[j];
                let d = match box_ {
                    Some(b) => b.min_image(raw),
                    None => raw,
                };
                let r = d.norm();
                if r <= 0.0 || r >= self.sym.r_cut {
                    continue;
                }
                let (gv, dg) = self.sym.radial(r);
                for t in 0..k {
                    g[i][t] += gv[t];
                    g[j][t] += gv[t];
                }
                pairs.push(PairTerm {
                    i,
                    j,
                    unit: d * (1.0 / r),
                    dg,
                });
            }
        }
        (g, pairs)
    }

    /// Total energy (kJ·mol⁻¹) for positions `pos` in `box_`.
    pub fn energy(&self, pos: &[Vec3], box_: Option<Box3>) -> f64 {
        let (g, _) = self.descriptors(pos, box_);
        g.iter()
            .map(|gi| self.forward(gi).last().expect("network has an output")[0])
            .sum()
    }

    /// Total energy (kJ·mol⁻¹) and forces (kJ·mol⁻¹·Å⁻¹) for positions `pos`.
    pub fn energy_forces(&self, pos: &[Vec3], box_: Option<Box3>) -> (Vec<Vec3>, f64) {
        let n = pos.len();
        let (g, pairs) = self.descriptors(pos, box_);
        let mut energy = 0.0;
        // dE/dG_i for every atom.
        let mut dg_atom = Vec::with_capacity(n);
        for gi in &g {
            let acts = self.forward(gi);
            energy += acts.last().expect("network has an output")[0];
            dg_atom.push(self.backward(&acts, 1.0, None));
        }
        let mut force = vec![Vec3::ZERO; n];
        for p in &pairs {
            // G_k(i) and G_k(j) both depend on r_ij.
            let de_dr: f64 = (0..self.sym.len())
                .map(|t| (dg_atom[p.i][t] + dg_atom[p.j][t]) * p.dg[t])
                .sum();
            let f = p.unit * (-de_dr);
            force[p.i] += f;
            force[p.j] -= f;
        }
        (force, energy)
    }

    /// Energy and forces for a [`System`], including its bonded topology.
    ///
    /// This is the force closure for
    /// [`VelocityVerlet::step_with`](crate::integrator::VelocityVerlet::step_with).
    pub fn evaluate(&self, sys: &System) -> (Vec<Vec3>, f64) {
        sys.with_bonded(self.energy_forces(&sys.pos, sys.box_))
    }

    /// Mean squared error of total energies over `data` (kJ²·mol⁻²).
    pub fn mse(&self, data: &[LabelledConfig]) -> f64 {
        if data.is_empty() {
            return 0.0;
        }
        let sum: f64 = data
            .iter()
            .map(|c| {
                let e = self.energy(&c.pos, c.box_);
                (e - c.energy) * (e - c.energy)
            })
            .sum();
        sum / data.len() as f64
    }

    /// Fit the network to labelled total energies by full-batch gradient
    /// descent on the mean squared error. Returns the MSE after training.
    pub fn train_energies(&mut self, data: &[LabelledConfig], epochs: usize, rate: f64) -> f64 {
        if data.is_empty() {
            return 0.0;
        }
        let scale = 2.0 / data.len() as f64;
        for _ in 0..epochs {
            let mut grads = self.zero_grads();
            for c in data {
                let (g, _) = self.descriptors(&c.pos, c.box_);
                let acts: Vec<Vec<Vec<f64>>> = g.iter().map(|gi| self.forward(gi)).collect();
                let energy: f64 = acts
                    .iter()
                    .map(|a| a.last().expect("network has an output")[0])
                    .sum();
                // dMSE/dE for this configuration, shared by all its atoms.
                let upstream = scale * (energy - c.energy);
                for a in &acts {
                    self.backward(a, upstream, Some(&mut grads));
                }
            }
            for (layer, (dw, db)) in self.layers.iter_mut().zip(&grads) {
                for (w, d) in layer.w.iter_mut().zip(dw) {
                    *w -= rate * d;
                }
                for (b, d) in layer.b.iter_mut().zip(db) {
                    *b -= rate * d;
                }
            }
        }
        self.mse(data)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forces::ForceModel;
    use crate::integrator::VelocityVerlet;
    use tpt_chem_core::forcefield::LennardJones;

    fn sym() -> RadialSymmetry {
        RadialSymmetry::uniform(0.8, 6.0, 6, 0.6)
    }

    fn cluster(seed: u64) -> Vec<Vec3> {
        let mut rng = Rng::new(seed);
        (0..5)
            .map(|_| {
                Vec3::new(
                    rng.range(-2.0, 2.0),
                    rng.range(-2.0, 2.0),
                    rng.range(-2.0, 2.0),
                )
            })
            .collect()
    }

    #[test]
    fn analytic_forces_match_finite_differences() {
        let mlp = RadialMlp::new(sym(), &[7, 5], 3);
        let pos = cluster(11);
        let (forces, _) = mlp.energy_forces(&pos, None);
        let h = 1e-6;
        for i in 0..pos.len() {
            for axis in 0..3 {
                let mut plus = pos.clone();
                let mut minus = pos.clone();
                let base = pos[i].get(axis);
                plus[i].set(axis, base + h);
                minus[i].set(axis, base - h);
                let fd = -(mlp.energy(&plus, None) - mlp.energy(&minus, None)) / (2.0 * h);
                let analytic = forces[i].get(axis);
                assert!(
                    (analytic - fd).abs() < 1e-6 * (1.0 + fd.abs()),
                    "atom {i} axis {axis}: analytic {analytic} vs fd {fd}"
                );
            }
        }
    }

    #[test]
    fn weight_gradients_match_finite_differences() {
        let mut mlp = RadialMlp::new(sym(), &[4], 5);
        let pos = cluster(21);
        let (g, _) = mlp.descriptors(&pos, None);
        let acts: Vec<Vec<Vec<f64>>> = g.iter().map(|gi| mlp.forward(gi)).collect();
        let mut grads = mlp.zero_grads();
        for a in &acts {
            mlp.backward(a, 1.0, Some(&mut grads));
        }
        // Check one weight in each layer against a central difference of the energy.
        let h = 1e-6;
        for (l, (dw, _)) in grads.iter().enumerate() {
            let idx = mlp.layers[l].w.len() / 2;
            let base = mlp.layers[l].w[idx];
            mlp.layers[l].w[idx] = base + h;
            let up = mlp.energy(&pos, None);
            mlp.layers[l].w[idx] = base - h;
            let down = mlp.energy(&pos, None);
            mlp.layers[l].w[idx] = base;
            let fd = (up - down) / (2.0 * h);
            let analytic = dw[idx];
            assert!(
                (analytic - fd).abs() < 1e-6 * (1.0 + fd.abs()),
                "layer {l}: analytic {analytic} vs fd {fd}"
            );
        }
    }

    #[test]
    fn energy_is_zero_beyond_the_cutoff() {
        // Two atoms farther apart than r_cut see no symmetry functions, so the
        // energy equals twice the network's output on an all-zero descriptor.
        let mlp = RadialMlp::new(sym(), &[3], 9);
        let far = [Vec3::ZERO, Vec3::new(10.0, 0.0, 0.0)];
        let isolated = mlp.forward(&vec![0.0; mlp.sym.len()]).last().unwrap()[0];
        assert!((mlp.energy(&far, None) - 2.0 * isolated).abs() < 1e-12);
    }

    #[test]
    fn nve_conserves_energy_with_the_learned_forces() {
        let mlp = RadialMlp::new(RadialSymmetry::uniform(0.8, 5.0, 8, 1.0), &[8], 4);
        let mut sys = System::new();
        let lj = LennardJones {
            sigma: 3.4,
            epsilon: 0.997,
        };
        for p in cluster(33) {
            sys.add_atom(lj, 0.0, 39.95, p * 1.2);
        }
        let integrator = VelocityVerlet::new(0.1);
        integrator.step_with(&mut sys, |s| mlp.evaluate(s));
        let e0 = sys.total_energy();
        let mut worst: f64 = 0.0;
        for _ in 0..300 {
            integrator.step_with(&mut sys, |s| mlp.evaluate(s));
            worst = worst.max((sys.total_energy() - e0).abs());
        }
        let scale = e0.abs().max(1.0);
        assert!(worst / scale < 1e-3, "drift {worst} on energy {e0}");
    }

    #[test]
    fn training_reduces_error_on_lennard_jones_energies() {
        // Reference: classical LJ energies of small random clusters.
        let lj = LennardJones {
            sigma: 3.4,
            epsilon: 0.997,
        };
        let model = ForceModel::AllPairs { cutoff: None };
        let mut rng = Rng::new(77);
        let data: Vec<LabelledConfig> = (0..24)
            .map(|_| {
                let pos: Vec<Vec3> = (0..4)
                    .map(|k| {
                        // Near-equilibrium spacing: keeps the LJ wall from dominating the targets.
                        let base = Vec3::new(k as f64 * 3.8, 0.0, 0.0);
                        base + Vec3::new(
                            rng.range(-0.3, 0.3),
                            rng.range(-0.3, 0.3),
                            rng.range(-0.3, 0.3),
                        )
                    })
                    .collect();
                let mut sys = System::new();
                for &p in &pos {
                    sys.add_atom(lj, 0.0, 39.95, p);
                }
                let energy = model.evaluate(&sys).1;
                LabelledConfig {
                    pos,
                    box_: None,
                    energy,
                }
            })
            .collect();

        let mut mlp = RadialMlp::new(RadialSymmetry::uniform(2.5, 6.5, 8, 1.0), &[10], 1);
        let before = mlp.mse(&data);
        let after = mlp.train_energies(&data, 3000, 0.02);
        assert!(after.is_finite());
        assert!(
            after < 0.01 * before,
            "MSE should drop by 100×: before {before}, after {after}"
        );
    }
}
