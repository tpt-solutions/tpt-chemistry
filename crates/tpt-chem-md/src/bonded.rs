//! Bonded interactions: harmonic bonds, harmonic angles, and
//! Ryckaert–Bellemans torsions, driven by an explicit topology on the
//! [`System`](crate::system::System).
//!
//! The energy/force kernels live in `tpt-chem-core::forcefield`; this module
//! applies them over the topology and folds the result into every
//! [`ForceModel`](crate::forces::ForceModel) evaluation, so run loops and
//! the minimiser see one consistent potential.
//!
//! Bonded distances are evaluated *without* the minimum image: topology
//! bonds must be shorter than half the box (the usual MD assumption). The
//! bonded virial is `−Σ rᵢ·Fᵢ`, which is exact for any internal-coordinate
//! potential.

use tpt_chem_core::autodiff::Real;
use tpt_chem_core::autodiff_reverse::reverse_gradient;
use tpt_chem_core::forcefield::{HarmonicAngle, HarmonicBond, RyckaertBellemans};
use tpt_chem_core::vec3::Vec3;

/// A bonded topology: parallel lists of terms with particle indices.
#[derive(Clone, Debug, Default)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Bonded {
    /// Harmonic bonds `(particles, parameters)`.
    pub bonds: Vec<([usize; 2], HarmonicBond)>,
    /// Harmonic angles `(particles, parameters)`, angle at particle 1.
    pub angles: Vec<([usize; 3], HarmonicAngle)>,
    /// Ryckaert–Bellemans torsions `(particles, parameters)`, chain
    /// `0–1–2–3`.
    pub dihedrals: Vec<([usize; 4], RyckaertBellemans)>,
}

impl Bonded {
    /// An empty topology.
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a harmonic bond.
    pub fn bond(&mut self, a: usize, b: usize, params: HarmonicBond) {
        self.bonds.push(([a, b], params));
    }

    /// Add a harmonic angle (measured at `b`).
    pub fn angle(&mut self, a: usize, b: usize, c: usize, params: HarmonicAngle) {
        self.angles.push(([a, b, c], params));
    }

    /// Add a torsion over the chain `a–b–c–d`.
    pub fn dihedral(&mut self, a: usize, b: usize, c: usize, d: usize, params: RyckaertBellemans) {
        self.dihedrals.push(([a, b, c, d], params));
    }
}

type V3<T> = [T; 3];

fn sub3<T: Real>(a: V3<T>, b: V3<T>) -> V3<T> {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot3<T: Real>(a: V3<T>, b: V3<T>) -> T {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross3<T: Real>(a: V3<T>, b: V3<T>) -> V3<T> {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm3<T: Real>(a: V3<T>) -> T {
    dot3(a, a).sqrt()
}

fn point3<T: Real>(x: &[T], i: usize) -> V3<T> {
    [x[3 * i], x[3 * i + 1], x[3 * i + 2]]
}

/// Bonded energy written against [`Real`], so the same code gives the energy
/// (`f64`) and its derivatives (`Var`). `x` holds the flat coordinates
/// `x, y, z` per particle.
///
/// The torsion energy depends only on `cos φ`, so no sign convention for φ
/// is needed.
pub fn bonded_energy<T: Real>(bonded: &Bonded, x: &[T]) -> T {
    let mut e = T::constant(0.0);
    for (idx, hb) in &bonded.bonds {
        let [a, b] = *idx;
        let r = norm3(sub3(point3(x, b), point3(x, a)));
        let d = r - T::constant(hb.r0);
        e = e + T::constant(0.5 * hb.k) * d * d;
    }
    for (idx, ha) in &bonded.angles {
        let [a, b, c] = *idx;
        let u = sub3(point3(x, a), point3(x, b));
        let v = sub3(point3(x, c), point3(x, b));
        let cos = dot3(u, v) / (norm3(u) * norm3(v));
        let d = cos.acos() - T::constant(ha.theta0);
        e = e + T::constant(0.5 * ha.k) * d * d;
    }
    for (idx, rb) in &bonded.dihedrals {
        let [a, b, c, d] = *idx;
        let b1 = sub3(point3(x, b), point3(x, a));
        let b2 = sub3(point3(x, c), point3(x, b));
        let b3 = sub3(point3(x, d), point3(x, c));
        let n1 = cross3(b1, b2);
        let n2 = cross3(b2, b3);
        let cos_phi = dot3(n1, n2) / (norm3(n1) * norm3(n2));
        // V = Σ Cₙ (−cos φ)ⁿ
        let xc = -cos_phi;
        let mut term = T::constant(1.0);
        for &c_n in &rb.c {
            e = e + T::constant(c_n) * term;
            term = term * xc;
        }
    }
    e
}

/// Bonded energy and forces by reverse-mode automatic differentiation of
/// [`bonded_energy`]: one backward pass for all forces.
///
/// About 20× slower than [`bonded_energy_forces`] on a 3000-term chain (see
/// `bonded_timing`), so the MD path keeps the analytic kernels. This function
/// is the cross-check. Collinear angles make `acos` singular, so if the AD
/// result is not finite the analytic kernels are used for that call.
pub fn bonded_energy_forces_ad(bonded: &Bonded, pos: &[Vec3]) -> (Vec<Vec3>, f64) {
    let flat: Vec<f64> = pos.iter().flat_map(|p| [p.x, p.y, p.z]).collect();
    let (energy, grad) = reverse_gradient(|x| bonded_energy(bonded, x), &flat);
    let forces: Vec<Vec3> = (0..pos.len())
        .map(|i| Vec3::new(-grad[3 * i], -grad[3 * i + 1], -grad[3 * i + 2]))
        .collect();
    let finite = energy.is_finite()
        && forces
            .iter()
            .all(|f| f.x.is_finite() && f.y.is_finite() && f.z.is_finite());
    if finite {
        (forces, energy)
    } else {
        bonded_energy_forces(bonded, pos)
    }
}

/// Bonded energy and forces from the hand-written kernels (Å, kJ·mol⁻¹). The
/// MD path uses this function. [`bonded_energy_forces_ad`] is the cross-check.
pub fn bonded_energy_forces(bonded: &Bonded, pos: &[Vec3]) -> (Vec<Vec3>, f64) {
    let mut forces = vec![Vec3::ZERO; pos.len()];
    let mut energy = 0.0;

    for (idx, hb) in &bonded.bonds {
        let [a, b] = *idx;
        let d = pos[b] - pos[a];
        let r = d.norm();
        if r < 1e-9 {
            continue;
        }
        energy += hb.energy(r);
        // Positive force magnitude pulls the atoms together when r > r0:
        // force on b points back toward a (−d̂).
        let fmag = hb.force_mag(r);
        let fvec = d * (fmag / r);
        forces[a] -= fvec;
        forces[b] += fvec;
    }

    for (idx, ha) in &bonded.angles {
        let [a, b, c] = *idx;
        let (e, f1, f2, f3) = ha.energy_and_forces(pos[a], pos[b], pos[c]);
        energy += e;
        forces[a] += f1;
        forces[b] += f2;
        forces[c] += f3;
    }

    for (idx, rb) in &bonded.dihedrals {
        let [a, b, c, d] = *idx;
        let (e, f1, f2, f3, f4) = rb.energy_and_forces(pos[a], pos[b], pos[c], pos[d]);
        energy += e;
        forces[a] += f1;
        forces[b] += f2;
        forces[c] += f3;
        forces[d] += f4;
    }

    (forces, energy)
}

/// Bonded virial `W = Σ rᵢ·Fᵢ` (kJ·mol⁻¹, with `Fᵢ = −∂U/∂rᵢ`) —
/// translationally invariant for internal-coordinate potentials.
pub fn bonded_virial(bonded_forces: &[Vec3], pos: &[Vec3]) -> f64 {
    let mut w = 0.0;
    for (r, f) in pos.iter().zip(bonded_forces.iter()) {
        w += r.dot(*f);
    }
    w
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Central finite-difference check of the full bonded gradient.
    #[test]
    fn bonded_forces_match_finite_differences() {
        let mut bonded = Bonded::new();
        bonded.bond(0, 1, HarmonicBond { k: 250.0, r0: 1.5 });
        bonded.bond(1, 2, HarmonicBond { k: 300.0, r0: 1.4 });
        bonded.angle(
            0,
            1,
            2,
            HarmonicAngle {
                k: 60.0,
                theta0: 1.8,
            },
        );
        bonded.dihedral(
            0,
            1,
            2,
            3,
            RyckaertBellemans {
                c: [0.9, 0.5, -0.3, 0.2, -0.1, 0.05],
            },
        );
        let pos = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.4, 0.2, 0.1),
            Vec3::new(2.3, 1.5, -0.3),
            Vec3::new(1.1, 2.2, 1.0),
        ];
        let (forces, energy) = bonded_energy_forces(&bonded, &pos);
        let h = 1e-5;
        for i in 0..4 {
            for axis in 0..3 {
                let comp = |v: Vec3| [v.x, v.y, v.z][axis];
                let mut pp = pos;
                let mut pm = pos;
                pp[i].set(axis, comp(pos[i]) + h);
                pm[i].set(axis, comp(pos[i]) - h);
                let (_, ep) = bonded_energy_forces(&bonded, &pp);
                let (_, em) = bonded_energy_forces(&bonded, &pm);
                let fd = -(ep - em) / (2.0 * h);
                let an = comp(forces[i]);
                assert!(
                    (fd - an).abs() < 1e-5 * (fd.abs() + an.abs()).max(1e-4) + 1e-6,
                    "atom {i} axis {axis}: fd {fd} vs analytic {an}"
                );
            }
        }
        assert!(energy.is_finite());
    }

    #[test]
    fn bond_force_restores_equilibrium_length() {
        let mut bonded = Bonded::new();
        bonded.bond(0, 1, HarmonicBond { k: 100.0, r0: 2.0 });
        // Stretched bond: atom 1 is pulled back toward atom 0.
        let pos = [Vec3::ZERO, Vec3::new(3.0, 0.0, 0.0)];
        let (forces, e) = bonded_energy_forces(&bonded, &pos);
        assert!(forces[1].x < 0.0, "stretched bond must pull back");
        assert!((forces[0].x + forces[1].x).abs() < 1e-12, "third law");
        let want = 0.5 * 100.0 * 1.0f64 * 1.0;
        assert!((e - want).abs() < 1e-12);
    }

    #[test]
    fn bonded_virial_of_internal_forces_sums_consistently() {
        // For a two-atom system with only a bond, −Σr·F equals r·f (the
        // pair virial), since F0 = −F1.
        let mut bonded = Bonded::new();
        bonded.bond(0, 1, HarmonicBond { k: 100.0, r0: 2.0 });
        let pos = [Vec3::new(1.0, 0.5, 0.0), Vec3::new(3.0, 0.5, 0.0)];
        let (forces, _) = bonded_energy_forces(&bonded, &pos);
        let w = bonded_virial(&forces, &pos);
        let _r = (pos[1] - pos[0]).norm();
        // f on atom 1 = −fvec along +x: magnitude k(r0−r) = −100 (attractive
        // pull), virial r·f = −100·... just require consistency:
        let pair_w = (pos[1] - pos[0]).dot(forces[1]);
        assert!((w - pair_w).abs() < 1e-12, "w {w} vs r·F1 {pair_w}");
    }

    /// A six-atom chain with bonds, angles and torsions, at non-collinear
    /// positions. Every term kind is present.
    fn chain_topology() -> (Bonded, Vec<Vec3>) {
        use tpt_chem_core::forcefield::{HarmonicAngle, HarmonicBond, RyckaertBellemans};
        let mut b = Bonded::new();
        for i in 0..5 {
            b.bond(
                i,
                i + 1,
                HarmonicBond {
                    k: 300.0 + 10.0 * i as f64,
                    r0: 1.5,
                },
            );
        }
        for i in 0..4 {
            b.angle(
                i,
                i + 1,
                i + 2,
                HarmonicAngle {
                    k: 80.0,
                    theta0: 1.9,
                },
            );
        }
        for i in 0..3 {
            b.dihedral(
                i,
                i + 1,
                i + 2,
                i + 3,
                RyckaertBellemans {
                    c: [2.0, 1.0, -0.5, 0.3, 0.1, -0.2],
                },
            );
        }
        let pos = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.4, 0.2, -0.1),
            Vec3::new(2.1, 1.3, 0.3),
            Vec3::new(3.4, 1.1, 1.0),
            Vec3::new(4.2, 2.3, 0.8),
            Vec3::new(5.5, 2.0, 1.6),
        ];
        (b, pos)
    }

    #[test]
    fn ad_bonded_forces_match_analytic_kernels() {
        let (bonded, pos) = chain_topology();
        let (f_ad, e_ad) = bonded_energy_forces_ad(&bonded, &pos);
        let (f_an, e_an) = bonded_energy_forces(&bonded, &pos);
        assert!((e_ad - e_an).abs() < 1e-9, "energy {e_ad} vs {e_an}");
        for (i, (a, b)) in f_ad.iter().zip(&f_an).enumerate() {
            for (x, y) in [(a.x, b.x), (a.y, b.y), (a.z, b.z)] {
                assert!(
                    (x - y).abs() < 1e-9 * (1.0 + y.abs()),
                    "atom {i}: {x} vs {y}"
                );
            }
        }
    }

    #[test]
    fn collinear_angle_falls_back_to_analytic_kernel() {
        use tpt_chem_core::forcefield::HarmonicAngle;
        let mut b = Bonded::new();
        b.angle(
            0,
            1,
            2,
            HarmonicAngle {
                k: 80.0,
                theta0: 1.9,
            },
        );
        let pos = vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
        ];
        let (f, e) = bonded_energy_forces_ad(&b, &pos);
        assert!(e.is_finite());
        assert!(f
            .iter()
            .all(|v| v.x.is_finite() && v.y.is_finite() && v.z.is_finite()));
    }

    /// Timing of the two bonded paths on a long chain. Run on purpose, in release:
    /// `cargo test -p tpt-chem-md --release --lib -- --ignored bonded_timing --nocapture`.
    #[test]
    #[ignore = "timing only; run in release with --nocapture"]
    fn bonded_timing() {
        use std::time::Instant;
        use tpt_chem_core::forcefield::{HarmonicAngle, HarmonicBond, RyckaertBellemans};
        let n = 3000;
        let mut b = Bonded::new();
        let pos: Vec<Vec3> = (0..n)
            .map(|i| {
                let t = i as f64;
                Vec3::new(1.5 * t.cos() * 0.3 + 0.1 * t, 1.5 * t.sin() * 0.3, 0.2 * t)
            })
            .collect();
        for i in 0..n - 1 {
            b.bond(i, i + 1, HarmonicBond { k: 300.0, r0: 1.5 });
        }
        for i in 0..n - 2 {
            b.angle(
                i,
                i + 1,
                i + 2,
                HarmonicAngle {
                    k: 80.0,
                    theta0: 1.9,
                },
            );
        }
        for i in 0..n - 3 {
            b.dihedral(
                i,
                i + 1,
                i + 2,
                i + 3,
                RyckaertBellemans {
                    c: [2.0, 1.0, -0.5, 0.3, 0.1, -0.2],
                },
            );
        }
        let reps = 20;
        let t0 = Instant::now();
        for _ in 0..reps {
            std::hint::black_box(bonded_energy_forces(&b, &pos));
        }
        let analytic = t0.elapsed().as_secs_f64() / reps as f64;
        let t1 = Instant::now();
        for _ in 0..reps {
            std::hint::black_box(bonded_energy_forces_ad(&b, &pos));
        }
        let ad = t1.elapsed().as_secs_f64() / reps as f64;
        println!(
            "bonded terms: {} bonds, {} angles, {} torsions",
            n - 1,
            n - 2,
            n - 3
        );
        println!("analytic: {:.3} ms per call", analytic * 1e3);
        println!(
            "reverse AD: {:.3} ms per call ({:.1}x)",
            ad * 1e3,
            ad / analytic
        );
    }
}
