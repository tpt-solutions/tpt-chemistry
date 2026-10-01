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
            let lj = pair_lj(&sys.lj[a], &sys.lj[b]);
            let sr6 = (lj.sigma * lj.sigma / r2).powi(3);
            let e_lj = 4.0 * lj.epsilon * (sr6 * sr6 - sr6);
            let f_lj = 24.0 * lj.epsilon / r * (2.0 * sr6 * sr6 - sr6);
            let qq = sys.charge[a] * sys.charge[b];
            let e_c = COULOMB_PREFACTOR_KJ_ANG * qq / r;
            let f_c = COULOMB_PREFACTOR_KJ_ANG * qq / r2;
            energy += e_lj + e_c;
            // Radial force magnitude (positive = repulsive) along d̂.
            let fmag = f_lj + f_c;
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
        let lj = pair_lj(&sys.lj[a], &sys.lj[b]);
        let sr6 = (lj.sigma * lj.sigma / r2).powi(3);
        let e_lj = 4.0 * lj.epsilon * (sr6 * sr6 - sr6);
        let f_lj = 24.0 * lj.epsilon / r * (2.0 * sr6 * sr6 - sr6);
        let qq = sys.charge[a] * sys.charge[b];
        let e_c = COULOMB_PREFACTOR_KJ_ANG * qq / r;
        let f_c = COULOMB_PREFACTOR_KJ_ANG * qq / r2;
        energy += e_lj + e_c;
        let fmag = f_lj + f_c;
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
