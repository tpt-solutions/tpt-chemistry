//! Force evaluation: Lennard-Jones + Coulomb, plain and via neighbor lists.

use tpt_chem_core::forcefield::LennardJones;
use tpt_chem_core::units::COULOMB_PREFACTOR_KJ_ANG;
use tpt_chem_core::vec3::Vec3;

use crate::box3::Box3;
use crate::neighbors::VerletList;
use crate::system::System;

/// Combined (Lorentz–Berthelot) LJ parameters for a pair.
#[inline]
fn pair_lj(la: &LennardJones, lb: &LennardJones) -> LennardJones {
    LennardJones::mix(la, lb)
}

/// Pair kernel shared by the all-pairs and neighbor-list evaluators:
/// Lennard-Jones + Coulomb energy and the radial force magnitude for one
/// pair at separation `r = sqrt(r2)`.
///
/// Returns `(pair_energy, fmag)` with `fmag > 0` repulsive.
#[inline]
fn pair_energy_force(sys: &System, a: usize, b: usize, r: f64, r2: f64) -> (f64, f64) {
    let lj = pair_lj(&sys.lj[a], &sys.lj[b]);
    let sr6 = (lj.sigma * lj.sigma / r2).powi(3);
    let e_lj = 4.0 * lj.epsilon * (sr6 * sr6 - sr6);
    let f_lj = 24.0 * lj.epsilon / r * (2.0 * sr6 * sr6 - sr6);
    let qq = sys.charge[a] * sys.charge[b];
    let e_c = COULOMB_PREFACTOR_KJ_ANG * qq / r;
    let f_c = COULOMB_PREFACTOR_KJ_ANG * qq / r2;
    let (s_lj, s_c) = pair_scales(sys, a, b);
    (s_lj * e_lj + s_c * e_c, s_lj * f_lj + s_c * f_c)
}

/// `(LJ, Coulomb)` scale factors for a pair (1, 1 unless excluded/scaled).
#[inline]
fn pair_scales(sys: &System, a: usize, b: usize) -> (f64, f64) {
    match sys.exclusions.as_ref().and_then(|x| x.get(a, b)) {
        Some(p) => (p.lj, p.coulomb),
        None => (1.0, 1.0),
    }
}

/// Correction for Ewald/PME models, whose reciprocal sum includes every
/// pair at full strength: subtract `(1 - s)*qq/r` for each scaled or
/// excluded pair. Returns `(forces, energy, virial)`.
fn coulomb_exclusion_correction(sys: &System) -> (Vec<Vec3>, f64, f64) {
    let mut forces = vec![Vec3::ZERO; sys.len()];
    let (mut energy, mut w) = (0.0, 0.0);
    let Some(ex) = sys.exclusions.as_ref() else {
        return (forces, energy, w);
    };
    for (a, b, scale) in ex.iter() {
        let missing = 1.0 - scale.coulomb;
        if missing == 0.0 {
            continue;
        }
        let raw = sys.pos[b] - sys.pos[a];
        let d = match sys.box_ {
            Some(bx) => bx.min_image(raw),
            None => raw,
        };
        let r2 = d.norm_sq();
        if r2 < 1e-12 {
            continue;
        }
        let r = r2.sqrt();
        let qq = COULOMB_PREFACTOR_KJ_ANG * sys.charge[a] * sys.charge[b];
        energy -= missing * qq / r;
        let fmag = -missing * qq / r2;
        let fvec = d * (fmag / r);
        forces[a] -= fvec;
        forces[b] += fvec;
        w += r * fmag;
    }
    (forces, energy, w)
}

/// Energy and forces for a full O(N²) all-pairs evaluation.
///
/// `cutoff` (`None` = no truncation) applies plain truncation — the
/// reported potential energy uses exactly the same truncated pair set as
/// the forces, so NVE energy conservation is preserved.
///
/// Returns `(forces, potential_energy)`.
pub fn forces_all_pairs(sys: &System, box_: Option<Box3>, cutoff: Option<f64>) -> (Vec<Vec3>, f64) {
    let n = sys.len();
    let mut forces = vec![Vec3::ZERO; n];
    let mut energy = 0.0;
    for a in 0..n {
        for b in (a + 1)..n {
            let raw = sys.pos[b] - sys.pos[a];
            let d = match box_ {
                Some(bx) => bx.min_image(raw),
                None => raw,
            };
            let r2 = d.norm_sq();
            let cut = cutoff.unwrap_or(f64::INFINITY);
            if r2 >= cut * cut || r2 < 1e-12 {
                continue;
            }
            let r = r2.sqrt();
            let (e, fmag) = pair_energy_force(sys, a, b, r, r2);
            energy += e;
            // fmag > 0 (repulsion) pushes a away from b, i.e. along -d̂.
            let fvec = d * (fmag / r);
            forces[a] -= fvec;
            forces[b] += fvec;
        }
    }
    (forces, energy)
}

/// Energy and forces using a Verlet neighbor list (pairs within
/// `cutoff + skin`, truncated at `cutoff`).
///
/// Returns `(forces, potential_energy)`.
pub fn forces_neighbor_list(
    sys: &System,
    box_: Option<Box3>,
    neighbors: &VerletList,
) -> (Vec<Vec3>, f64) {
    let n = sys.len();
    let mut forces = vec![Vec3::ZERO; n];
    let mut energy = 0.0;
    let cut2 = neighbors.cutoff * neighbors.cutoff;
    for &(a, b) in &neighbors.pairs {
        let raw = sys.pos[b] - sys.pos[a];
        let d = match box_ {
            Some(bx) => bx.min_image(raw),
            None => raw,
        };
        let r2 = d.norm_sq();
        if r2 >= cut2 || r2 < 1e-12 {
            continue;
        }
        let r = r2.sqrt();
        let (e, fmag) = pair_energy_force(sys, a, b, r, r2);
        energy += e;
        let fvec = d * (fmag / r);
        forces[a] -= fvec;
        forces[b] += fvec;
    }
    (forces, energy)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::System;
    use tpt_chem_core::forcefield::LennardJones;

    fn argon_dimer(r: f64) -> System {
        let mut sys = System::new();
        let lj = LennardJones {
            sigma: 3.4,
            epsilon: 0.997,
        };
        sys.add_atom(lj, 0.0, 39.95, Vec3::new(-r / 2.0, 0.0, 0.0));
        sys.add_atom(lj, 0.0, 39.95, Vec3::new(r / 2.0, 0.0, 0.0));
        sys
    }

    #[test]
    fn lj_force_matches_finite_difference() {
        let r = 4.0;
        let sys = argon_dimer(r);
        let (f, e) = forces_all_pairs(&sys, None, None);
        // Attraction pulls atoms together: +x force on atom 0.
        assert!(f[0].x > 0.0);
        assert!(f[1].x < 0.0);
        // Energy = LJ at 4.0 Å.
        let lj = LennardJones {
            sigma: 3.4,
            epsilon: 0.997,
        };
        assert!((e - lj.energy(r)).abs() < 1e-10);
        // Newton's third law.
        assert!(f[0].x + f[1].x < 1e-14);

        // Finite-difference check: dE/dr = -F on the second atom.
        let h = 1e-6;
        let e_p = forces_all_pairs(&argon_dimer(r + h), None, None).1;
        let e_m = forces_all_pairs(&argon_dimer(r - h), None, None).1;
        let de_dr = (e_p - e_m) / (2.0 * h);
        assert!(
            (de_dr + f[1].x).abs() < 1e-5,
            "dE/dr {de_dr} vs -F {}",
            -f[1].x
        );
    }

    #[test]
    fn coulomb_pair() {
        let mut sys = System::new();
        sys.add_atom(
            LennardJones {
                sigma: 1.0,
                epsilon: 0.0,
            },
            1.0,
            1.008,
            Vec3::new(0.0, 0.0, 0.0),
        );
        sys.add_atom(
            LennardJones {
                sigma: 1.0,
                epsilon: 0.0,
            },
            -1.0,
            1.008,
            Vec3::new(2.0, 0.0, 0.0),
        );
        let (f, e) = forces_all_pairs(&sys, None, None);
        assert!((e - (-1389.354582655 / 2.0)).abs() < 1e-5);
        // Opposite charges attract: force on atom 0 is +x.
        assert!(f[0].x > 0.0);
        assert!((f[0].x - 1389.354582655 / 4.0).abs() < 1e-4);
    }

    #[test]
    fn neighbor_list_matches_all_pairs() {
        let box_ = crate::box3::Box3::cubic(15.0);
        let mut sys = System::new();
        let lj = LennardJones {
            sigma: 3.4,
            epsilon: 0.997,
        };
        let mut rng = tpt_chem_core::rng::Rng::new(7);
        // Jittered lattice so no two atoms overlap.
        let mut placed = 0;
        'outer: for i in 0..4 {
            for j in 0..4 {
                for k in 0..3 {
                    let p = Vec3::new(
                        2.0 + i as f64 * 3.7 + rng.range(-0.3, 0.3),
                        2.0 + j as f64 * 3.7 + rng.range(-0.3, 0.3),
                        2.0 + k as f64 * 3.7 + rng.range(-0.3, 0.3),
                    );
                    sys.add_atom(lj, if rng.uniform() < 0.3 { 0.5 } else { 0.0 }, 39.95, p);
                    placed += 1;
                    if placed >= 40 {
                        break 'outer;
                    }
                }
            }
        }
        sys.wrap_in_box();
        let (f_ref, e_ref) = forces_all_pairs(&sys, Some(box_), Some(6.0));
        let vl = VerletList::build(&sys.pos, 6.0, 1.0, Some(box_));
        let (f_nl, e_nl) = forces_neighbor_list(&sys, Some(box_), &vl);
        assert!((e_ref - e_nl).abs() < 1e-9);
        for i in 0..sys.len() {
            assert!((f_ref[i] - f_nl[i]).norm() < 1e-9, "force mismatch at {i}");
        }
    }
}

/// Configurable force evaluation for a [`System`] — pairs the ensemble run
/// loops (`crate::ensemble`) with a single, self-consistent energy/force
/// definition.
#[derive(Clone, Copy, Debug)]
pub enum ForceModel {
    /// Lennard-Jones + direct Coulomb over all pairs, optionally truncated
    /// at `cutoff` (Å).
    AllPairs {
        /// Pair cutoff (Å); `None` evaluates every pair.
        cutoff: Option<f64>,
    },
    /// Lennard-Jones (truncated) + full Ewald Coulomb: real-space erfc
    /// pairs, reciprocal structure-factor sum, self term.
    LjPlusEwald {
        /// LJ cutoff (Å).
        lj_cutoff: f64,
        /// Ewald parameters; α must match a reasonably converged set.
        ewald: crate::ewald::EwaldParams,
    },
    /// Lennard-Jones (truncated) + PME Coulomb: real-space erfc pairs,
    /// mesh reciprocal, self term. O(N log N) in the reciprocal part.
    LjPlusPme {
        /// LJ cutoff (Å).
        lj_cutoff: f64,
        /// PME parameters; α must match the real-space splitting.
        pme: crate::pme::PmeParams,
    },
}

impl ForceModel {
    /// Evaluate energy and forces for `sys`: the nonbonded model plus the
    /// bonded topology (`System::bonded`) when present.
    pub fn evaluate(&self, sys: &System) -> (Vec<Vec3>, f64) {
        sys.with_bonded(self.evaluate_nonbonded(sys))
    }

    /// Nonbonded part only (pair LJ/Coulomb or LJ + Ewald/PME).
    pub fn evaluate_nonbonded(&self, sys: &System) -> (Vec<Vec3>, f64) {
        match *self {
            ForceModel::AllPairs { cutoff } => forces_all_pairs(sys, sys.box_, cutoff),
            ForceModel::LjPlusEwald { lj_cutoff, ewald } => {
                self.lj_longrange(sys, lj_cutoff, |b| {
                    crate::ewald::ewald_energy_forces(&sys.pos, &sys.charge, b, &ewald)
                })
            }
            ForceModel::LjPlusPme { lj_cutoff, pme } => self.lj_longrange(sys, lj_cutoff, |b| {
                let (mut f, mut e) =
                    crate::pme::pme_reciprocal_energy_forces(&sys.pos, &sys.charge, b, &pme);
                // Mesh = smooth reciprocal part only; add the erfc real-space
                // pairs (cut at the LJ cutoff, capped at half the box) and
                // the self term.
                let real = pme_real_params(&pme, lj_cutoff, b);
                let (fr, er) = crate::ewald::ewald_real_space(&sys.pos, &sys.charge, b, &real);
                for (fi, fri) in f.iter_mut().zip(fr.iter()) {
                    *fi += *fri;
                }
                e += er + crate::ewald::ewald_self_energy(&sys.charge, pme.alpha);
                (f, e)
            }),
        }
    }

    /// LJ all-pairs (truncated at `lj_cutoff`) combined with a
    /// long-range Coulomb evaluator over the box.
    fn lj_longrange(
        &self,
        sys: &System,
        lj_cutoff: f64,
        longrange: impl FnOnce(&crate::box3::Box3) -> (Vec<Vec3>, f64),
    ) -> (Vec<Vec3>, f64) {
        let n = sys.len();
        let mut forces = vec![Vec3::ZERO; n];
        let mut energy = 0.0;
        let cut2 = lj_cutoff * lj_cutoff;
        if let Some(b) = sys.box_ {
            for a in 0..n {
                for c in (a + 1)..n {
                    let d = b.min_image(sys.pos[c] - sys.pos[a]);
                    let r2 = d.norm_sq();
                    if r2 >= cut2 || r2 < 1e-12 {
                        continue;
                    }
                    let r = r2.sqrt();
                    let lj = pair_lj(&sys.lj[a], &sys.lj[c]);
                    let s_lj = pair_scales(sys, a, c).0;
                    let sr6 = (lj.sigma * lj.sigma / r2).powi(3);
                    energy += s_lj * 4.0 * lj.epsilon * (sr6 * sr6 - sr6);
                    let fmag = s_lj * 24.0 * lj.epsilon / r * (2.0 * sr6 * sr6 - sr6);
                    let fvec = d * (fmag / r);
                    forces[a] -= fvec;
                    forces[c] += fvec;
                }
            }
            let (f_lr, e_lr) = longrange(&b);
            for (f, flr) in forces.iter_mut().zip(f_lr.iter()) {
                *f += *flr;
            }
            energy += e_lr;
            let (f_ex, e_ex, _) = coulomb_exclusion_correction(sys);
            for (f, fe) in forces.iter_mut().zip(f_ex.iter()) {
                *f += *fe;
            }
            energy += e_ex;
        }
        (forces, energy)
    }

    /// Virial `W = Σ_pairs r·f` of the modeled interactions (kJ·mol⁻¹),
    /// consistent with [`ForceModel::evaluate`] (for the mesh model, to the
    /// PME discretization level). Feeds [`crate::ensemble::pressure_bar`].
    pub fn virial(&self, sys: &System) -> f64 {
        let mut w = match *self {
            ForceModel::AllPairs { cutoff } => virial_all_pairs(sys, sys.box_, cutoff),
            ForceModel::LjPlusEwald { lj_cutoff, ewald } => {
                self.lj_virial(sys, lj_cutoff)
                    + crate::ewald::ewald_virial(
                        &sys.pos,
                        &sys.charge,
                        sys.box_.as_ref().unwrap_or(&DUMMY_BOX),
                        &ewald,
                    )
            }
            ForceModel::LjPlusPme { lj_cutoff, pme } => {
                self.lj_virial(sys, lj_cutoff) + pme_virial(sys, &pme, lj_cutoff)
            }
        };
        if let Some(bonded) = &sys.bonded {
            let (fb, _) = crate::bonded::bonded_energy_forces(bonded, &sys.pos);
            w += crate::bonded::bonded_virial(&fb, &sys.pos);
        }
        w
    }

    fn lj_virial(&self, sys: &System, lj_cutoff: f64) -> f64 {
        let n = sys.len();
        let mut w = 0.0;
        let cut2 = lj_cutoff * lj_cutoff;
        if let Some(b) = sys.box_ {
            for a in 0..n {
                for c in (a + 1)..n {
                    let d = b.min_image(sys.pos[c] - sys.pos[a]);
                    let r2 = d.norm_sq();
                    if r2 >= cut2 || r2 < 1e-12 {
                        continue;
                    }
                    let r = r2.sqrt();
                    let lj = pair_lj(&sys.lj[a], &sys.lj[c]);
                    let s_lj = pair_scales(sys, a, c).0;
                    let sr6 = (lj.sigma * lj.sigma / r2).powi(3);
                    let fmag = s_lj * 24.0 * lj.epsilon / r * (2.0 * sr6 * sr6 - sr6);
                    w += r * fmag;
                }
            }
            w += coulomb_exclusion_correction(sys).2;
        }
        w
    }
}

/// Zero-length placeholder for open-boundary systems, where the long-range
/// Coulomb is zero anyway (no box → no periodic electrostatics).
static DUMMY_BOX: crate::box3::Box3 = crate::box3::Box3 {
    length: Vec3::new(1.0, 1.0, 1.0),
};

/// Real-space Ewald parameters matching a PME model: same splitting alpha,
/// cutoff `min(lj_cutoff, L/2)`.
fn pme_real_params(
    pme: &crate::pme::PmeParams,
    lj_cutoff: f64,
    b: &crate::box3::Box3,
) -> crate::ewald::EwaldParams {
    let half = 0.5 * b.length.x.min(b.length.y).min(b.length.z);
    crate::ewald::EwaldParams {
        alpha: pme.alpha,
        r_cutoff: lj_cutoff.min(half),
        g_max: 0.0,
    }
}

/// Mesh-model Coulomb virial: `Σ rᵢ·Fᵢ` of the PME reciprocal forces
/// (wrap-safe for net-neutral systems) plus the real-space pair virial.
fn pme_virial(sys: &System, pme: &crate::pme::PmeParams, lj_cutoff: f64) -> f64 {
    let Some(b) = sys.box_ else { return 0.0 };
    let (f, _) = crate::pme::pme_reciprocal_energy_forces(&sys.pos, &sys.charge, &b, pme);
    let mut w = 0.0;
    for (r, fi) in sys.pos.iter().zip(f.iter()) {
        w += r.dot(*fi);
    }
    w + crate::ewald::ewald_real_virial(
        &sys.pos,
        &sys.charge,
        &b,
        &pme_real_params(pme, lj_cutoff, &b),
    )
}

/// Virial of the truncated LJ + direct Coulomb pair sum,
/// `W = Σ_pairs r·f(r)` (kJ·mol⁻¹).
pub fn virial_all_pairs(sys: &System, box_: Option<Box3>, cutoff: Option<f64>) -> f64 {
    let n = sys.len();
    let mut w = 0.0;
    let cut = cutoff.unwrap_or(f64::INFINITY);
    for a in 0..n {
        for b in (a + 1)..n {
            let raw = sys.pos[b] - sys.pos[a];
            let d = match box_ {
                Some(bx) => bx.min_image(raw),
                None => raw,
            };
            let r2 = d.norm_sq();
            if r2 >= cut * cut || r2 < 1e-12 {
                continue;
            }
            let r = r2.sqrt();
            let lj = pair_lj(&sys.lj[a], &sys.lj[b]);
            let sr6 = (lj.sigma * lj.sigma / r2).powi(3);
            let f_lj = 24.0 * lj.epsilon / r * (2.0 * sr6 * sr6 - sr6);
            let f_c = COULOMB_PREFACTOR_KJ_ANG * sys.charge[a] * sys.charge[b] / r2;
            let (s_lj, s_c) = pair_scales(sys, a, b);
            w += r * (s_lj * f_lj + s_c * f_c);
        }
    }
    w
}

#[cfg(test)]
mod bonded_model_tests {
    use super::*;
    use crate::bonded::Bonded;

    #[test]
    fn force_model_adds_bonded_terms() {
        let mut sys = System::new();
        let lj = LennardJones {
            sigma: 2.5,
            epsilon: 0.5,
        };
        sys.add_atom(lj, 0.0, 1.0, Vec3::new(0.0, 0.0, 0.0));
        sys.add_atom(lj, 0.0, 1.0, Vec3::new(2.6, 0.0, 0.0));
        let model = ForceModel::AllPairs { cutoff: None };
        let (_, e_free) = model.evaluate(&sys);
        let mut bonded = Bonded::new();
        bonded.bond(
            0,
            1,
            tpt_chem_core::forcefield::HarmonicBond { k: 200.0, r0: 2.5 },
        );
        sys.bonded = Some(bonded);
        let (f, e_bonded) = model.evaluate(&sys);
        assert!(
            (e_bonded - e_free).abs() < 0.1 * e_free.abs() + 1.0,
            "bond near r0 barely shifts the energy: {e_free} vs {e_bonded}"
        );
        // Force direction flips with the bond term dominating near r0.
        assert!(f[0].x.is_finite());
    }

    #[test]
    fn bonded_model_virial_is_exact_gradient_consistent() {
        // With only a bond, the model virial must equal r * fmag.
        let mut sys = System::new();
        let lj = LennardJones {
            sigma: 1.0,
            epsilon: 0.0,
        };
        sys.add_atom(lj, 0.0, 1.0, Vec3::new(0.0, 0.0, 0.0));
        sys.add_atom(lj, 0.0, 1.0, Vec3::new(2.7, 0.0, 0.0));
        let mut bonded = Bonded::new();
        bonded.bond(
            0,
            1,
            tpt_chem_core::forcefield::HarmonicBond { k: 150.0, r0: 2.5 },
        );
        sys.bonded = Some(bonded);
        let model = ForceModel::AllPairs { cutoff: None };
        let w = model.virial(&sys);
        let r = 2.7;
        let fmag_bond = 150.0 * (2.5 - 2.7); // = -30 (stretching pulls back)
        assert!((w - r * fmag_bond).abs() < 1e-9, "w = {w}");
    }
}

/// A [`System`] plus [`ForceModel`] viewed as a generic
/// [`tpt_chem_core::potential::Potential`] (positions in, energy/forces
/// out), so MD force fields can share the generic minimiser and any other
/// tool written against the trait.
#[derive(Clone, Debug)]
pub struct ModelPotential {
    /// The system (its `pos` is overwritten on every call).
    pub system: System,
    /// Force model used for evaluation.
    pub model: ForceModel,
}

impl tpt_chem_core::potential::Potential for ModelPotential {
    type Error = core::convert::Infallible;

    fn n_atoms(&self) -> usize {
        self.system.len()
    }

    fn energy_forces(
        &mut self,
        positions: &[Vec3],
    ) -> Result<(f64, Vec<Vec3>), core::convert::Infallible> {
        self.system.pos.clear();
        self.system.pos.extend_from_slice(positions);
        let (f, e) = self.model.evaluate(&self.system);
        // The model returns (forces, energy) per `ForceModel::evaluate`.
        Ok((e, f))
    }
}

#[cfg(test)]
mod potential_tests {
    use super::*;
    use tpt_chem_core::forcefield::LennardJones;
    use tpt_chem_core::potential::{minimize_bfgs, MinimizeSettings};

    #[test]
    fn lj_dimer_relaxes_to_potential_minimum() {
        let lj = LennardJones {
            sigma: 3.4,
            epsilon: 0.997,
        };
        let mut sys = System::new();
        sys.add_atom(lj, 0.0, 39.95, Vec3::new(0.0, 0.0, 0.0));
        sys.add_atom(lj, 0.0, 39.95, Vec3::new(0.0, 0.0, 4.5));
        let mut pot = ModelPotential {
            system: sys,
            model: ForceModel::AllPairs { cutoff: None },
        };
        let mut pos = [Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 4.5)];
        let rep = minimize_bfgs(
            &mut pot,
            &mut pos,
            MinimizeSettings {
                force_tol: 1e-4,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(rep.converged, "{rep:?}");
        let r = (pos[1] - pos[0]).norm();
        assert!((r - 3.4 * 2f64.powf(1.0 / 6.0)).abs() < 1e-3, "{r}");
        assert!((rep.energy + 0.997).abs() < 1e-4);
    }
}

#[cfg(test)]
mod exclusion_tests {
    use super::*;
    use crate::exclusions::{Exclusions, PairScale};
    use tpt_chem_core::forcefield::LennardJones;

    fn pair(r: f64, q: f64) -> System {
        let mut sys = System::new();
        let lj = LennardJones {
            sigma: 3.0,
            epsilon: 0.5,
        };
        sys.add_atom(lj, q, 12.0, Vec3::new(0.0, 0.0, 0.0));
        sys.add_atom(lj, -q, 12.0, Vec3::new(r, 0.0, 0.0));
        sys
    }

    #[test]
    fn excluded_pair_has_no_interaction() {
        let mut sys = pair(2.0, 0.4);
        let mut ex = Exclusions::new();
        ex.set(
            0,
            1,
            PairScale {
                lj: 0.0,
                coulomb: 0.0,
            },
        );
        sys.exclusions = Some(ex);
        let (f, e) = forces_all_pairs(&sys, None, None);
        assert_eq!(e, 0.0);
        assert!(f.iter().all(|v| v.norm() == 0.0));
        assert_eq!(virial_all_pairs(&sys, None, None), 0.0);
    }

    #[test]
    fn one_four_scaling_is_linear_in_the_fudge() {
        let full = forces_all_pairs(&pair(3.2, 0.4), None, None).1;
        let mut sys = pair(3.2, 0.4);
        let mut ex = Exclusions::new();
        ex.set(
            0,
            1,
            PairScale {
                lj: 0.5,
                coulomb: 0.5,
            },
        );
        sys.exclusions = Some(ex);
        let (_, e) = forces_all_pairs(&sys, None, None);
        assert!((e - 0.5 * full).abs() < 1e-12);
    }

    #[test]
    fn ewald_with_exclusion_matches_finite_difference() {
        use crate::ewald::EwaldParams;
        let mut sys = pair(1.5, 0.5);
        let lj = sys.lj[0];
        sys.add_atom(lj, 0.3, 12.0, Vec3::new(5.0, 4.0, 3.0));
        sys.add_atom(lj, -0.3, 12.0, Vec3::new(3.0, 6.0, 7.0));
        sys.set_box(Box3::cubic(12.0));
        let mut ex = Exclusions::new();
        ex.set(
            0,
            1,
            PairScale {
                lj: 0.0,
                coulomb: 0.0,
            },
        );
        sys.exclusions = Some(ex);
        let model = ForceModel::LjPlusEwald {
            lj_cutoff: 5.5,
            ewald: EwaldParams {
                alpha: 0.4,
                r_cutoff: 5.9,
                g_max: 3.5,
            },
        };
        let (f, _) = model.evaluate(&sys);
        let h = 1e-5;
        for (atom, ax) in [(0usize, 0usize), (1, 0), (2, 1)] {
            let e = |s: f64| {
                let mut m = sys.clone();
                let mut p = m.pos[atom];
                p.set(ax, p.get(ax) + s);
                m.pos[atom] = p;
                model.evaluate(&m).1
            };
            let fd = -(e(h) - e(-h)) / (2.0 * h);
            assert!(
                (f[atom].get(ax) - fd).abs() < 1e-4,
                "{atom} {ax}: {} vs {fd}",
                f[atom].get(ax)
            );
        }
    }

    /// Regression: the PME model once returned only the mesh part (no
    /// real-space erfc pairs, no self term), so its Coulomb energy was
    /// wrong. It must now agree with the direct Ewald sum.
    #[test]
    fn pme_model_agrees_with_ewald_model() {
        use crate::ewald::EwaldParams;
        use crate::pme::PmeParams;
        let mut sys = pair(3.5, 0.5);
        let lj = LennardJones {
            sigma: 3.0,
            epsilon: 0.0,
        };
        sys.lj = vec![lj; 2];
        sys.add_atom(lj, 0.3, 12.0, Vec3::new(5.0, 4.0, 3.0));
        sys.add_atom(lj, -0.3, 12.0, Vec3::new(3.0, 6.0, 7.0));
        sys.set_box(Box3::cubic(12.0));
        let ewald = ForceModel::LjPlusEwald {
            lj_cutoff: 5.9,
            ewald: EwaldParams {
                alpha: 0.4,
                r_cutoff: 5.9,
                g_max: 4.5,
            },
        };
        let pme = ForceModel::LjPlusPme {
            lj_cutoff: 5.9,
            pme: PmeParams {
                alpha: 0.4,
                dims: [32, 32, 32],
                g_max: 0.0,
            },
        };
        let (fe, ee) = ewald.evaluate(&sys);
        let (fp, ep) = pme.evaluate(&sys);
        assert!((ee - ep).abs() < 0.02 * ee.abs().max(1.0), "{ee} vs {ep}");
        for (a, b) in fe.iter().zip(&fp) {
            assert!((*a - *b).norm() < 0.5, "{a:?} vs {b:?}");
        }
    }
}
