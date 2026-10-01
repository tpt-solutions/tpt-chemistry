//! Force-field definitions: Lennard-Jones, Coulomb, harmonic bonds/angles,
//! and Ryckaert–Bellemans dihedrals (spec.txt §4, `tpt-chem-core`).
//!
//! All quantities are in the MD unit system: positions in Å, energies in
//! kJ·mol⁻¹, charges in elementary units e. Radial force magnitudes are
//! signed along `r̂` from the interaction center: **positive = repulsive**
//! (pushing the pair apart).
//!
//! ```
//! use tpt_chem_core::forcefield::LennardJones;
//!
//! // Argon-ish LJ pair in Å / kJ·mol⁻¹.
//! let lj = LennardJones { sigma: 3.4, epsilon: 0.997 };
//! let (e, f) = lj.energy_and_force(lj.r_min());
//! assert!(e < 0.0);          // well depth at the minimum
//! assert!(f.abs() < 1e-12);  // zero force at the minimum
//! ```

use alloc::string::String;
use alloc::vec::Vec;

use crate::num;
use crate::vec3::Vec3;

/// Coulomb pair interaction `q₁q₂/(4πε₀r)` between two point charges.
///
/// Charges are in units of e, distance in Å, energy in kJ·mol⁻¹.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Coulomb {
    /// Charge of the first particle (e).
    pub qa: f64,
    /// Charge of the second particle (e).
    pub qb: f64,
}

impl Coulomb {
    /// Pair energy at distance `r` (Å) in kJ·mol⁻¹.
    pub fn energy(&self, r: f64) -> f64 {
        crate::units::COULOMB_PREFACTOR_KJ_ANG * self.qa * self.qb / r
    }

    /// Radial force magnitude at distance `r` (positive = repulsive).
    pub fn force_mag(&self, r: f64) -> f64 {
        crate::units::COULOMB_PREFACTOR_KJ_ANG * self.qa * self.qb / (r * r)
    }

    /// Energy and force magnitude in one call.
    pub fn energy_and_force(&self, r: f64) -> (f64, f64) {
        (self.energy(r), self.force_mag(r))
    }
}

/// 12-6 Lennard-Jones pair potential.
///
/// `V(r) = 4ε[(σ/r)¹² − (σ/r)⁶]`, σ in Å, ε in kJ·mol⁻¹.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LennardJones {
    /// Collision diameter σ (Å).
    pub sigma: f64,
    /// Well depth ε (kJ·mol⁻¹).
    pub epsilon: f64,
}

impl LennardJones {
    /// Radius of the potential minimum, `2^(1/6) σ`.
    pub fn r_min(&self) -> f64 {
        num::powf(2.0, 1.0 / 6.0) * self.sigma
    }

    /// Pair energy at distance `r` (kJ·mol⁻¹).
    pub fn energy(&self, r: f64) -> f64 {
        let sr6 = num::powi(self.sigma / r, 6);
        4.0 * self.epsilon * (sr6 * sr6 - sr6)
    }

    /// Radial force magnitude at distance `r` (positive = repulsive).
    pub fn force_mag(&self, r: f64) -> f64 {
        let sr6 = num::powi(self.sigma / r, 6);
        24.0 * self.epsilon / r * (2.0 * sr6 * sr6 - sr6)
    }

    /// Energy and force magnitude in one call.
    pub fn energy_and_force(&self, r: f64) -> (f64, f64) {
        (self.energy(r), self.force_mag(r))
    }

    /// Lorentz–Berthelot mixing: `σᵢⱼ = (σᵢ+σⱼ)/2`, `εᵢⱼ = √(εᵢεⱼ)`.
    pub fn mix(a: &LennardJones, b: &LennardJones) -> LennardJones {
        LennardJones {
            sigma: 0.5 * (a.sigma + b.sigma),
            epsilon: num::sqrt(a.epsilon * b.epsilon),
        }
    }
}

/// A combined pair potential: Lennard-Jones plus Coulomb, the standard
/// non-bonded interaction of a force field.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PairPotential {
    /// Dispersion/repulsion parameters.
    pub lj: LennardJones,
    /// Product of the two partial charges (e²) — `qa * qb`.
    pub qq: f64,
}

impl PairPotential {
    /// Build from two atom parameters (charges and LJ terms).
    pub fn new(lj: LennardJones, qa: f64, qb: f64) -> Self {
        PairPotential { lj, qq: qa * qb }
    }

    /// Lorentz–Berthelot-mixed LJ with the plain Coulomb product.
    pub fn mixed(lj_a: &LennardJones, qa: f64, lj_b: &LennardJones, qb: f64) -> Self {
        PairPotential {
            lj: LennardJones::mix(lj_a, lj_b),
            qq: qa * qb,
        }
    }

    /// Untruncated pair energy at distance `r` (kJ·mol⁻¹).
    pub fn energy(&self, r: f64) -> f64 {
        self.lj.energy(r) + crate::units::COULOMB_PREFACTOR_KJ_ANG * self.qq / r
    }

    /// Untruncated radial force magnitude (positive = repulsive).
    pub fn force_mag(&self, r: f64) -> f64 {
        self.lj.force_mag(r) + crate::units::COULOMB_PREFACTOR_KJ_ANG * self.qq / (r * r)
    }

    /// Energy and force magnitude in one call.
    pub fn energy_and_force(&self, r: f64) -> (f64, f64) {
        (self.energy(r), self.force_mag(r))
    }
}

/// Harmonic bond: `V(r) = ½k(r − r₀)²`.
///
/// `k` in kJ·mol⁻¹·Å⁻², `r₀` in Å.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicBond {
    /// Force constant (kJ·mol⁻¹·Å⁻²).
    pub k: f64,
    /// Equilibrium bond length (Å).
    pub r0: f64,
}

impl HarmonicBond {
    /// Bond energy at length `r`.
    pub fn energy(&self, r: f64) -> f64 {
        let d = r - self.r0;
        0.5 * self.k * d * d
    }

    /// Magnitude of the restoring force at length `r`, signed so positive
    /// pulls the atoms together when `r > r₀`.
    pub fn force_mag(&self, r: f64) -> f64 {
        self.k * (self.r0 - r)
    }
}

/// Harmonic angle: `V(θ) = ½k(θ − θ₀)²`.
///
/// `k` in kJ·mol⁻¹·rad⁻², `θ₀` in radians.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HarmonicAngle {
    /// Force constant (kJ·mol⁻¹·rad⁻²).
    pub k: f64,
    /// Equilibrium angle (rad).
    pub theta0: f64,
}

impl HarmonicAngle {
    /// Angle energy at angle `theta` (rad).
    pub fn energy(&self, theta: f64) -> f64 {
        let d = theta - self.theta0;
        0.5 * self.k * d * d
    }

    /// Angle energy and the gradients ∂V/∂p₁, ∂V/∂p₂, ∂V/∂p₃ for the atom
    /// triple (p₁, p₂, p₃) with the angle measured at p₂.
    pub fn energy_and_forces(&self, p1: Vec3, p2: Vec3, p3: Vec3) -> (f64, Vec3, Vec3, Vec3) {
        let theta = crate::coords::angle(p1, p2, p3);
        let dv_dtheta = self.k * (theta - self.theta0);

        // u, v: arms from the vertex p2.
        let u = p1 - p2;
        let v = p3 - p2;
        let un = u.norm();
        let vn = v.norm();
        let sin_t = num::sin(theta);
        if sin_t < 1e-12 {
            // Collinear: gradient is numerically singular; zero force is the
            // physically expected limit for a finite dV/dθ.
            return (self.energy(theta), Vec3::ZERO, Vec3::ZERO, Vec3::ZERO);
        }
        let dtheta_du = -(v / (un * vn) - u * (u.dot(v) / (un * un * un * vn))) / sin_t;
        let dtheta_dv = -(u / (un * vn) - v * (u.dot(v) / (un * vn * vn * vn))) / sin_t;
        let f1 = -dv_dtheta * dtheta_du;
        let f3 = -dv_dtheta * dtheta_dv;
        (self.energy(theta), f1, -f1 - f3, f3)
    }
}

/// Ryckaert–Bellemans torsion potential:
/// `V(ψ) = Σₙ Cₙ cosⁿ ψ` with `ψ = φ − 180°` (φ the IUPAC dihedral), so the
/// trans state (φ = 180°) sits at `Σ Cₙ`.
///
/// Coefficients C₀..C₅ in kJ·mol⁻¹.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RyckaertBellemans {
    /// Expansion coefficients C₀..C₅ (kJ·mol⁻¹).
    pub c: [f64; 6],
}

impl RyckaertBellemans {
    /// Torsion energy at IUPAC dihedral `phi` (rad).
    pub fn energy(&self, phi: f64) -> f64 {
        // cos(ψ) = cos(φ − π) = −cos φ.
        let x = -num::cos(phi);
        let mut v = 0.0;
        let mut xpow = 1.0;
        for c in &self.c {
            v += c * xpow;
            xpow *= x;
        }
        v
    }

    /// dV/dφ at IUPAC dihedral `phi` (rad).
    pub fn denergy_dphi(&self, phi: f64) -> f64 {
        // V(ψ) = Σ cₙ xⁿ with x = cos ψ = −cos φ.
        // dV/dφ = Σ n·cₙ·xⁿ⁻¹ · dx/dφ, and dx/dφ = sin φ.
        let x = -num::cos(phi);
        let mut dv = 0.0;
        let mut xpow = 1.0; // x^(n-1)
        for (n, c) in self.c.iter().enumerate().skip(1) {
            dv += f64::from(n as u32) * c * xpow;
            xpow *= x;
        }
        dv * num::sin(phi)
    }

    /// Torsion energy and the forces on the four atoms of the chain
    /// `p1–p2–p3–p4` whose dihedral is `phi`.
    ///
    /// Returns `(V, F1, F2, F3, F4)`; the forces sum to zero.
    #[allow(clippy::too_many_arguments)]
    pub fn energy_and_forces(
        &self,
        p1: Vec3,
        p2: Vec3,
        p3: Vec3,
        p4: Vec3,
    ) -> (f64, Vec3, Vec3, Vec3, Vec3) {
        let phi = crate::coords::dihedral(p1, p2, p3, p4);
        let dv = self.denergy_dphi(phi);

        // Dihedral gradient (validated against finite differences):
        // with r_ij = p2−p1, r_kj = p3−p2, r_lk = p4−p3,
        // m = r_ij×r_kj, n = r_kj×r_lk, a = r_ij·r_kj/|r_kj|²,
        // b = r_lk·r_kj/|r_kj|²:
        //   ∂φ/∂p1 =  (|r_kj|/|m|²) m
        //   ∂φ/∂p4 = −(|r_kj|/|n|²) n
        //   ∂φ/∂p2 = −(a+1)·∂φ/∂p1 + b·∂φ/∂p4
        //   ∂φ/∂p3 = −(∂φ/∂p1 + ∂φ/∂p2 + ∂φ/∂p4)   (translation invariance)
        // Forces are F = −(dV/dφ)·∂φ/∂p.
        let r_ij = p2 - p1;
        let r_kj = p3 - p2;
        let r_lk = p4 - p3;
        let m = r_ij.cross(r_kj);
        let n = r_kj.cross(r_lk);
        let m2 = m.norm_sq();
        let n2 = n.norm_sq();
        let kj2 = r_kj.norm_sq();
        if m2 < 1e-24 || n2 < 1e-24 || kj2 < 1e-24 {
            return (self.energy(phi), Vec3::ZERO, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO);
        }
        let kjn = r_kj.norm();
        let g_i = (kjn / m2) * m;
        let g_l = -(kjn / n2) * n;
        let a = r_ij.dot(r_kj) / kj2;
        let b = r_lk.dot(r_kj) / kj2;
        let g_j = -(a + 1.0) * g_i + b * g_l;
        let g_k = -(g_i + g_j + g_l);
        (
            self.energy(phi),
            -dv * g_i,
            -dv * g_j,
            -dv * g_k,
            -dv * g_l,
        )
    }
}

/// One force-field atom type: a name plus its non-bonded parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct AtomType {
    /// Type name, e.g. `"opls_135"` or `"C3"`.
    pub name: String,
    /// Lennard-Jones parameters.
    pub lj: LennardJones,
    /// Partial charge (e).
    pub charge: f64,
    /// Standard atomic mass (amu), usually the element mass.
    pub mass: f64,
}

/// A force field: atom types plus bonded parameter tables keyed by type
/// indices (normalized ascending within each key).
#[derive(Clone, Debug, Default)]
pub struct ForceField {
    /// Atom-type table.
    pub atom_types: Vec<AtomType>,
    /// Harmonic bond parameters by `(type_i, type_j)`.
    pub bonds: alloc::collections::BTreeMap<[u16; 2], HarmonicBond>,
    /// Harmonic angle parameters by `(type_i, type_j, type_k)`.
    pub angles: alloc::collections::BTreeMap<[u16; 3], HarmonicAngle>,
    /// RB dihedral parameters by `(type_i, type_j, type_k, type_l)`.
    pub dihedrals: alloc::collections::BTreeMap<[u16; 4], RyckaertBellemans>,
}

impl ForceField {
    /// An empty force field.
    pub fn new() -> Self {
        Self::default()
    }

    /// Register an atom type and return its index.
    pub fn add_atom_type(&mut self, name: &str, lj: LennardJones, charge: f64, mass: f64) -> u16 {
        self.atom_types.push(AtomType {
            name: name.into(),
            lj,
            charge,
            mass,
        });
        (self.atom_types.len() - 1) as u16
    }

    /// Set the harmonic bond parameter for a pair of atom types.
    pub fn set_bond(&mut self, t1: u16, t2: u16, p: HarmonicBond) {
        self.bonds.insert(norm2(t1, t2), p);
    }

    /// Set the harmonic angle parameter for an atom-type triple.
    pub fn set_angle(&mut self, t1: u16, t2: u16, t3: u16, p: HarmonicAngle) {
        self.angles.insert(norm3(t1, t2, t3), p);
    }

    /// Set the RB dihedral parameter for an atom-type quadruple.
    pub fn set_dihedral(&mut self, t1: u16, t2: u16, t3: u16, t4: u16, p: RyckaertBellemans) {
        self.dihedrals.insert(norm4(t1, t2, t3, t4), p);
    }

    /// Bond parameter for two atom types, if defined.
    pub fn bond(&self, t1: u16, t2: u16) -> Option<HarmonicBond> {
        self.bonds.get(&norm2(t1, t2)).copied()
    }

    /// Angle parameter for a triple, if defined.
    pub fn angle(&self, t1: u16, t2: u16, t3: u16) -> Option<HarmonicAngle> {
        self.angles.get(&norm3(t1, t2, t3)).copied()
    }

    /// Dihedral parameter for a quadruple, if defined.
    pub fn dihedral(&self, t1: u16, t2: u16, t3: u16, t4: u16) -> Option<RyckaertBellemans> {
        self.dihedrals.get(&norm4(t1, t2, t3, t4)).copied()
    }

    /// Mixed (Lorentz–Berthelot) pair potential for two atom types.
    ///
    /// # Panics
    /// If either type index is out of range.
    pub fn pair_potential(&self, t1: u16, t2: u16) -> PairPotential {
        let a = &self.atom_types[t1 as usize];
        let b = &self.atom_types[t2 as usize];
        PairPotential::mixed(&a.lj, a.charge, &b.lj, b.charge)
    }
}

fn norm2(a: u16, b: u16) -> [u16; 2] {
    if a <= b {
        [a, b]
    } else {
        [b, a]
    }
}

/// Canonical angle key: reversed rather than sorted, so the vertex position
/// survives (O–C–O and C–O–C are different angle types).
fn norm3(a: u16, b: u16, c: u16) -> [u16; 3] {
    if a <= c {
        [a, b, c]
    } else {
        [c, b, a]
    }
}

/// Canonical dihedral key: reversed rather than sorted, so the two central
/// atoms stay central.
fn norm4(a: u16, b: u16, c: u16, d: u16) -> [u16; 4] {
    if (a, b, c, d) <= (d, c, b, a) {
        [a, b, c, d]
    } else {
        [d, c, b, a]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-10;

    #[test]
    fn lj_minimum_and_well_depth() {
        let lj = LennardJones { sigma: 3.4, epsilon: 0.997 };
        let rmin = lj.r_min();
        assert!((rmin - 3.4 * num::powf(2.0, 1.0 / 6.0)).abs() < TOL);
        let (e, f) = lj.energy_and_force(rmin);
        assert!((e + lj.epsilon).abs() < 1e-12);
        assert!(f.abs() < 1e-12);
        // Repulsive inside σ (V > 0), attractive near the well.
        assert!(lj.energy(lj.sigma * 0.9) > 0.0);
        assert!(lj.force_mag(lj.sigma * 0.9) > 0.0);
        assert!(lj.energy(rmin * 1.5) < 0.0);
        assert!(lj.force_mag(rmin * 1.5) < 0.0);
        // Finite-difference check: E'(r) = -F(r).
        let h = 1e-6;
        for r in [3.0, 3.6, 4.5, 5.5] {
            let d = (lj.energy(r + h) - lj.energy(r - h)) / (2.0 * h);
            assert!((d + lj.force_mag(r)).abs() < 1e-6, "r = {r}");
        }
    }

    #[test]
    fn coulomb_values() {
        let repulsive = Coulomb { qa: 1.0, qb: 1.0 };
        // Two unit charges 1 Å apart: 1389.35 kJ/mol.
        assert!((repulsive.energy(1.0) - 1389.354582655).abs() < 1e-6);
        assert!(repulsive.energy(2.0) > 0.0);
        let c = Coulomb { qa: 1.0, qb: -1.0 };
        assert!(c.energy(1.0) < 0.0);
        // dE/dr = -F
        let h = 1e-6;
        let d = (c.energy(2.0 + h) - c.energy(2.0 - h)) / (2.0 * h);
        assert!((d + c.force_mag(2.0)).abs() < 1e-5);
    }

    #[test]
    fn pair_potential_superposition() {
        let p = PairPotential::new(LennardJones { sigma: 3.0, epsilon: 1.0 }, 0.5, -0.2);
        let r = 3.3;
        assert!((p.energy(r) - (p.lj.energy(r) + Coulomb { qa: 0.5, qb: -0.2 }.energy(r))).abs() < 1e-12);
    }

    #[test]
    fn mixing_lorentz_berthelot() {
        let a = LennardJones { sigma: 3.0, epsilon: 0.5 };
        let b = LennardJones { sigma: 4.0, epsilon: 2.0 };
        let m = LennardJones::mix(&a, &b);
        assert!((m.sigma - 3.5).abs() < TOL);
        assert!((m.epsilon - 1.0).abs() < TOL);
    }

    #[test]
    fn harmonic_bond() {
        let b = HarmonicBond { k: 1500.0, r0: 1.1 };
        assert!((b.energy(1.1)).abs() < 1e-15);
        assert!((b.energy(1.2) - 0.5 * 1500.0 * 0.01).abs() < 1e-9);
        assert!((b.force_mag(1.2) + 150.0).abs() < 1e-9); // pulls inward
        assert!(b.force_mag(1.0) > 0.0); // pushes outward
    }

    #[test]
    fn harmonic_angle_gradients_by_finite_difference() {
        let ang = HarmonicAngle { k: 80.0, theta0: 109.47f64.to_radians() };
        let p2 = Vec3::new(0.1, -0.3, 0.2);
        let p1 = p2 + Vec3::new(1.0, 0.4, -0.2);
        let p3 = p2 + Vec3::new(-0.3, 0.8, 0.9);
        let (e, f1, f2, f3) = ang.energy_and_forces(p1, p2, p3);
        let forces = [f1, f2, f3];
        let points = [p1, p2, p3];
        let h = 1e-7;
        for i in 0..3 {
            for c in 0..3 {
                let mut pp = points;
                pp[i].set(c, points[i].get(c) + h);
                let e_plus = ang.energy_and_forces(pp[0], pp[1], pp[2]).0;
                pp[i].set(c, points[i].get(c) - h);
                let e_minus = ang.energy_and_forces(pp[0], pp[1], pp[2]).0;
                let fd = (e_plus - e_minus) / (2.0 * h);
                assert!(
                    (fd + forces[i].get(c)).abs() < 1e-5,
                    "atom {i} comp {c}: fd {fd} vs {}",
                    forces[i].get(c)
                );
            }
        }
        let _ = e;
    }

    #[test]
    fn rb_dihedral_values_and_gradient() {
        // A classic RB expansion of an OPLS-style torsion (C0..C5).
        let rb = RyckaertBellemans { c: [0.7, 1.5, 0.3, -0.6, 0.1, 0.0] };
        // energy(π) = Σ cₙ (trans state)
        let etrans = rb.energy(core::f64::consts::PI);
        assert!((etrans - 2.0).abs() < 1e-12);

        // Finite-difference dV/dφ.
        let h = 1e-7;
        for phi in [0.5, 1.0, 2.0, -1.3] {
            let fd = (rb.energy(phi + h) - rb.energy(phi - h)) / (2.0 * h);
            assert!((fd - rb.denergy_dphi(phi)).abs() < 1e-6, "phi = {phi}");
        }

        // Atomic-force check by finite differences.
        let p1 = Vec3::new(0.0, 0.0, 0.0);
        let p2 = Vec3::new(1.4, 0.0, 0.0);
        let p3 = p2 + Vec3::new(0.5, 1.3, 0.0);
        let p4 = p3 + Vec3::new(-0.6, 1.0, 0.9);
        let (_, f1, f2, f3, f4) = rb.energy_and_forces(p1, p2, p3, p4);
        let forces = [f1, f2, f3, f4];
        let points = [p1, p2, p3, p4];
        for i in 0..4 {
            for c in 0..3 {
                let mut pp = points;
                pp[i].set(c, points[i].get(c) + h);
                let e_plus = rb.energy_and_forces(pp[0], pp[1], pp[2], pp[3]).0;
                pp[i].set(c, points[i].get(c) - h);
                let e_minus = rb.energy_and_forces(pp[0], pp[1], pp[2], pp[3]).0;
                let fd = (e_plus - e_minus) / (2.0 * h);
                assert!(
                    (fd + forces[i].get(c)).abs() < 1e-4,
                    "atom {i} comp {c}: fd {fd} vs {}",
                    forces[i].get(c)
                );
            }
        }
        // Translational invariance: ΣF = 0.
        let sum = f1 + f2 + f3 + f4;
        assert!(sum.norm() < 1e-10);
    }

    #[test]
    fn force_field_container() {
        let mut ff = ForceField::new();
        let t_c = ff.add_atom_type(
            "CT",
            LennardJones { sigma: 3.5, epsilon: 0.276 },
            -0.18,
            12.011,
        );
        let t_o = ff.add_atom_type(
            "OT",
            LennardJones { sigma: 3.07, epsilon: 0.65 },
            -0.64,
            15.999,
        );
        let t_h = ff.add_atom_type(
            "HC",
            LennardJones { sigma: 2.5, epsilon: 0.126 },
            0.06,
            1.008,
        );
        ff.set_bond(t_c, t_o, HarmonicBond { k: 1600.0, r0: 1.43 });
        ff.set_angle(t_c, t_o, t_h, HarmonicAngle { k: 100.0, theta0: 109.5 });
        ff.set_dihedral(t_c, t_c, t_o, t_o, RyckaertBellemans { c: [0.6; 6] });
        assert!((ff.bond(t_o, t_c).unwrap().r0 - 1.43).abs() < TOL); // reversed key
        assert!((ff.angle(t_h, t_o, t_c).unwrap().k - 100.0).abs() < TOL); // reversed triple
        // Reversal normalization: O–O–C–C ≡ C–C–O–O.
        assert!(ff.dihedral(t_o, t_o, t_c, t_c).is_some());
        // C–O–H ≠ H–O–C is fine, but C–O–C was never defined.
        assert!(ff.angle(t_c, t_o, t_c).is_none());
        let p = ff.pair_potential(t_c, t_o);
        assert!((p.qq - 0.1152).abs() < 1e-9);
    }
}
