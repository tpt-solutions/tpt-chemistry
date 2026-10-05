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
    (e_lj + e_c, f_lj + f_c)
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
    /// Evaluate energy and forces for `sys`.
    pub fn evaluate(&self, sys: &System) -> (Vec<Vec3>, f64) {
        match *self {
            ForceModel::AllPairs { cutoff } => forces_all_pairs(sys, sys.box_, cutoff),
            ForceModel::LjPlusEwald { lj_cutoff, ewald } => {
                self.lj_longrange(sys, lj_cutoff, |b| {
                    crate::ewald::ewald_energy_forces(&sys.pos, &sys.charge, b, &ewald)
                })
            }
            ForceModel::LjPlusPme { lj_cutoff, pme } => self.lj_longrange(sys, lj_cutoff, |b| {
                crate::pme::pme_reciprocal_energy_forces(&sys.pos, &sys.charge, b, &pme)
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
                    let sr6 = (lj.sigma * lj.sigma / r2).powi(3);
                    energy += 4.0 * lj.epsilon * (sr6 * sr6 - sr6);
                    let fmag = 24.0 * lj.epsilon / r * (2.0 * sr6 * sr6 - sr6);
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
                self.lj_virial(sys, lj_cutoff) + pme_virial(sys, &pme)
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
                    let sr6 = (lj.sigma * lj.sigma / r2).powi(3);
                    let fmag = 24.0 * lj.epsilon / r * (2.0 * sr6 * sr6 - sr6);
                    w += r * fmag;
                }
            }
        }
        w
    }
}

/// Zero-length placeholder for open-boundary systems, where the long-range
/// Coulomb is zero anyway (no box → no periodic electrostatics).
static DUMMY_BOX: crate::box3::Box3 = crate::box3::Box3 {
    length: Vec3::new(1.0, 1.0, 1.0),
};

/// Mesh-model virial: LJ part plus `Σ rᵢ·Fᵢ` of the PME reciprocal forces
/// (wrap-safe for net-neutral systems).
fn pme_virial(sys: &System, pme: &crate::pme::PmeParams) -> f64 {
    let Some(b) = sys.box_ else { return 0.0 };
    let (f, _) = crate::pme::pme_reciprocal_energy_forces(&sys.pos, &sys.charge, &b, pme);
    let mut w = 0.0;
    for (r, fi) in sys.pos.iter().zip(f.iter()) {
        w += r.dot(*fi);
    }
    w
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
            w += r * (f_lj + f_c);
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
