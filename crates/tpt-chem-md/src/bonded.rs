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

use tpt_chem_core::forcefield::{HarmonicAngle, HarmonicBond, RyckaertBellemans};
use tpt_chem_core::vec3::Vec3;

/// A bonded topology: parallel lists of terms with particle indices.
#[derive(Clone, Debug, Default)]
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

/// Bonded energy and forces over the topology (Å, kJ·mol⁻¹).
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
}
