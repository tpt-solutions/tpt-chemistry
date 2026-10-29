//! An AMBER-style force field: harmonic bonds and angles, periodic torsions,
//! and Lennard-Jones plus Coulomb nonbonded terms with 1-4 scaling.
//!
//! The functional forms and units follow OpenMM's conventions, so parameters
//! from OpenMM's AMBER files apply directly (Å, kJ/mol, elementary charges).
//! Atom types are not assigned here: the caller supplies per-atom charges and
//! Lennard-Jones parameters, and the bonded terms as explicit lists.
//!
//! - bond: `E = ½ k (r − r₀)²`
//! - angle: `E = ½ k (θ − θ₀)²`
//! - torsion: `E = k (1 + cos(n φ − φ₀))`, one term per listed quadruple
//!   (impropers are ordinary torsion entries)
//! - nonbonded: Lennard-Jones `4ε[(σ/r)¹² − (σ/r)⁶]` with arithmetic σ and
//!   geometric ε, and Coulomb. Pairs one or two bonds apart are excluded.
//!   Pairs three bonds apart (1-4) have their Coulomb term divided by `scee`
//!   and their LJ well depth divided by `scnb`.
//!
//! Energies are written against [`Real`], so forces come from automatic
//! differentiation with no hand derivatives. [`AmberModel::energy_forces`] uses
//! reverse mode, which gives every force in one backward pass.
//! [`AmberModel::energy_forces_forward`] uses forward mode, one pass per
//! coordinate, and is kept as an independent check.

use std::collections::VecDeque;

use tpt_chem_core::autodiff::{gradient, Real};
use tpt_chem_core::autodiff_reverse::reverse_gradient;
use tpt_chem_core::units::COULOMB_PREFACTOR_KJ_ANG;
use tpt_chem_core::vec3::Vec3;

/// Harmonic bond between atoms `i` and `j`.
#[derive(Clone, Debug, PartialEq)]
pub struct BondTerm {
    /// First atom.
    pub i: usize,
    /// Second atom.
    pub j: usize,
    /// Equilibrium length (Å).
    pub r0: f64,
    /// Force constant (kJ·mol⁻¹·Å⁻²).
    pub k: f64,
}

/// Harmonic angle `i–j–k` (vertex `j`).
#[derive(Clone, Debug, PartialEq)]
pub struct AngleTerm {
    /// First atom.
    pub i: usize,
    /// Vertex atom.
    pub j: usize,
    /// Third atom.
    pub k: usize,
    /// Equilibrium angle (rad).
    pub theta0: f64,
    /// Force constant (kJ·mol⁻¹·rad⁻²).
    pub k_theta: f64,
}

/// Ryckaert–Bellemans torsion over four atoms in order.
///
/// `E = Σₙ Cₙ cosⁿ ψ` with `ψ = φ − π`, the form OPLS-AA uses. Since
/// `cos ψ = −cos φ`, only the dihedral cosine is needed.
#[derive(Clone, Debug, PartialEq)]
pub struct RbTorsionTerm {
    /// Atoms of the quadruple, in order.
    pub atoms: [usize; 4],
    /// Coefficients `C₀ … C₅` (kJ/mol).
    pub c: [f64; 6],
}

/// Periodic torsion over four atoms in order.
#[derive(Clone, Debug, PartialEq)]
pub struct TorsionTerm {
    /// Atoms of the quadruple, in order.
    pub atoms: [usize; 4],
    /// Periodicity `n` (≥ 1).
    pub periodicity: u32,
    /// Phase `φ₀` (rad).
    pub phase: f64,
    /// Barrier height `k` (kJ/mol).
    pub k: f64,
}

/// Relation of a pair of atoms in the bond graph.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PairKind {
    /// One or two bonds apart: no nonbonded interaction.
    Excluded,
    /// Three bonds apart (shortest path): scaled.
    OneFour,
    /// Four or more bonds apart (or unconnected): full.
    Full,
}

/// An AMBER-style system: atoms, bonded terms and scaling factors.
#[derive(Clone, Debug, PartialEq)]
pub struct AmberModel {
    /// Partial charges (e).
    pub charges: Vec<f64>,
    /// Lennard-Jones σ (Å).
    pub sigma: Vec<f64>,
    /// Lennard-Jones ε (kJ/mol).
    pub epsilon: Vec<f64>,
    /// Harmonic bonds.
    pub bonds: Vec<BondTerm>,
    /// Harmonic angles.
    pub angles: Vec<AngleTerm>,
    /// Periodic torsions.
    pub torsions: Vec<TorsionTerm>,
    /// Ryckaert–Bellemans torsions (OPLS-AA form).
    pub rb_torsions: Vec<RbTorsionTerm>,
    /// Coulomb divisor for 1-4 pairs (AMBER: 1.2).
    pub scee: f64,
    /// Lennard-Jones well-depth divisor for 1-4 pairs (AMBER: 2.0).
    pub scnb: f64,
}

type V3<T> = [T; 3];

fn sub<T: Real>(a: V3<T>, b: V3<T>) -> V3<T> {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot<T: Real>(a: V3<T>, b: V3<T>) -> T {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross<T: Real>(a: V3<T>, b: V3<T>) -> V3<T> {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm<T: Real>(a: V3<T>) -> T {
    dot(a, a).sqrt()
}

fn flat_coordinates(positions: &[Vec3]) -> Vec<f64> {
    positions.iter().flat_map(|p| [p.x, p.y, p.z]).collect()
}

fn forces_from_gradient(n_atoms: usize, grad: &[f64]) -> Vec<Vec3> {
    (0..n_atoms)
        .map(|i| Vec3::new(-grad[3 * i], -grad[3 * i + 1], -grad[3 * i + 2]))
        .collect()
}

fn point<T: Real>(x: &[T], i: usize) -> V3<T> {
    [x[3 * i], x[3 * i + 1], x[3 * i + 2]]
}

/// `(cos φ, sin φ)` of the dihedral through atoms `atoms` (IUPAC sign).
///
/// Built from unnormalized components, so no `atan2` is needed:
/// `cos φ ∝ n₁·n₂` and `sin φ ∝ |b₂| b₁·n₂`, with `n₁ = b₁×b₂` and `n₂ = b₂×b₃`.
fn dihedral_cos_sin<T: Real>(x: &[T], atoms: [usize; 4]) -> (T, T) {
    let p = atoms.map(|i| point(x, i));
    let b1 = sub(p[1], p[0]);
    let b2 = sub(p[2], p[1]);
    let b3 = sub(p[3], p[2]);
    let n1 = cross(b1, b2);
    let n2 = cross(b2, b3);
    let cx = dot(n1, n2);
    let sy = norm(b2) * dot(b1, n2);
    let r = (cx * cx + sy * sy).sqrt();
    (cx / r, sy / r)
}

impl AmberModel {
    /// Classify every atom pair from the bond graph (shortest path ≤ 3 bonds).
    fn pair_kinds(&self) -> Vec<Vec<PairKind>> {
        let n = self.charges.len();
        let mut adjacency = vec![Vec::new(); n];
        for b in &self.bonds {
            adjacency[b.i].push(b.j);
            adjacency[b.j].push(b.i);
        }
        let mut kinds = vec![vec![PairKind::Full; n]; n];
        for start in 0..n {
            let mut dist = vec![usize::MAX; n];
            dist[start] = 0;
            let mut queue = VecDeque::from([start]);
            while let Some(u) = queue.pop_front() {
                if dist[u] >= 3 {
                    continue;
                }
                for &v in &adjacency[u] {
                    if dist[v] == usize::MAX {
                        dist[v] = dist[u] + 1;
                        queue.push_back(v);
                    }
                }
            }
            for end in 0..n {
                kinds[start][end] = match dist[end] {
                    1 | 2 => PairKind::Excluded,
                    3 => PairKind::OneFour,
                    _ => PairKind::Full,
                };
            }
        }
        kinds
    }

    /// Energy (kJ/mol) for flat coordinates `x = [x₀, y₀, z₀, x₁, …]` (Å).
    ///
    /// Generic over [`Real`], so the same code gives the energy (`f64`) and
    /// its derivatives (`Dual`).
    pub fn energy<T: Real>(&self, x: &[T]) -> T {
        let kinds = self.pair_kinds();
        let mut e = T::constant(0.0);

        for b in &self.bonds {
            let d = norm(sub(point(x, b.i), point(x, b.j))) - T::constant(b.r0);
            e = e + T::constant(0.5 * b.k) * d * d;
        }

        for a in &self.angles {
            let u = sub(point(x, a.i), point(x, a.j));
            let v = sub(point(x, a.k), point(x, a.j));
            let cos_theta = dot(u, v) / (norm(u) * norm(v));
            let d = cos_theta.acos() - T::constant(a.theta0);
            e = e + T::constant(0.5 * a.k_theta) * d * d;
        }

        for rb in &self.rb_torsions {
            let [i, j, k, l] = rb.atoms;
            let (cos_phi, _) = dihedral_cos_sin(x, [i, j, k, l]);
            // cos ψ = cos(φ − π) = −cos φ
            let cos_psi = -cos_phi;
            let mut power = T::constant(1.0);
            for &coef in &rb.c {
                e = e + T::constant(coef) * power;
                power = power * cos_psi;
            }
        }

        for t in &self.torsions {
            let [i, j, k, l] = t.atoms;
            let (c, s) = dihedral_cos_sin(x, [i, j, k, l]);
            // cos(nφ), sin(nφ) by the recurrence in (cos φ, sin φ).
            let (mut cn, mut sn) = (c, s);
            for _ in 1..t.periodicity {
                let next_c = cn * c - sn * s;
                let next_s = sn * c + cn * s;
                cn = next_c;
                sn = next_s;
            }
            // cos(nφ − φ₀) = cos nφ cos φ₀ + sin nφ sin φ₀
            let term = cn * T::constant(t.phase.cos()) + sn * T::constant(t.phase.sin());
            e = e + T::constant(t.k) * (T::constant(1.0) + term);
        }

        for (i, kind_row) in kinds.iter().enumerate() {
            for (j, &kind) in kind_row.iter().enumerate().skip(i + 1) {
                let scale = match kind {
                    PairKind::Excluded => continue,
                    PairKind::OneFour => Some((self.scee, self.scnb)),
                    PairKind::Full => None,
                };
                let d = norm(sub(point(x, i), point(x, j)));
                let sigma = T::constant(0.5 * (self.sigma[i] + self.sigma[j]));
                let mut eps = (self.epsilon[i] * self.epsilon[j]).sqrt();
                let mut qq = self.charges[i] * self.charges[j] * COULOMB_PREFACTOR_KJ_ANG;
                if let Some((scee, scnb)) = scale {
                    eps /= scnb;
                    qq /= scee;
                }
                let s6 = (sigma / d).powi(6);
                let lj = T::constant(4.0 * eps) * (s6 * s6 - s6);
                let coulomb = T::constant(qq) / d;
                e = e + lj + coulomb;
            }
        }
        e
    }

    /// Energy (kJ/mol) and forces (kJ/mol/Å) for atomic positions in Å.
    ///
    /// Forces are exact derivatives of [`AmberModel::energy`], from one
    /// reverse-mode pass for all Cartesian coordinates.
    pub fn energy_forces(&self, positions: &[Vec3]) -> (f64, Vec<Vec3>) {
        let flat = flat_coordinates(positions);
        let (energy, grad) = reverse_gradient(|x| self.energy(x), &flat);
        (energy, forces_from_gradient(positions.len(), &grad))
    }

    /// Same as [`AmberModel::energy_forces`], but with forward mode (one pass per
    /// coordinate). Kept as an independent check of the reverse-mode path.
    pub fn energy_forces_forward(&self, positions: &[Vec3]) -> (f64, Vec<Vec3>) {
        let flat = flat_coordinates(positions);
        let (energy, grad) = gradient(|x| self.energy(x), &flat);
        (energy, forces_from_gradient(positions.len(), &grad))
    }
}
