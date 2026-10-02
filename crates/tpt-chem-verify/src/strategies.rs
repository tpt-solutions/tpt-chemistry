//! `proptest` strategies for generating *physically valid* inputs:
//! non-overlapping geometries, sensible force-field parameters, and
//! reasonable temperatures (spec.txt §4, `tpt-chem-verify`).

use proptest::prelude::*;

use tpt_chem_core::forcefield::LennardJones;
use tpt_chem_core::rng::Rng;
use tpt_chem_core::vec3::Vec3;
use tpt_chem_md::box3::Box3;
use tpt_chem_md::system::System;

/// Lennard-Jones parameters in physically plausible ranges:
/// σ ∈ [2.5, 4.0] Å, ε ∈ [0.05, 2.0] kJ·mol⁻¹.
pub fn arb_lj_params() -> impl Strategy<Value = LennardJones> {
    (2.5f64..4.0, 0.05f64..2.0).prop_map(|(sigma, epsilon)| LennardJones { sigma, epsilon })
}

/// A mass in amu within [1, 60].
pub fn arb_mass() -> impl Strategy<Value = f64> {
    1.0f64..60.0
}

/// A temperature in K within [1, 600].
pub fn arb_temperature() -> impl Strategy<Value = f64> {
    1.0f64..600.0
}

/// A deterministic seed for the in-house RNG.
pub fn arb_seed() -> impl Strategy<Value = u64> {
    any::<u64>()
}

/// A jittered simple-cubic lattice of `n` atoms with the given `spacing`
/// (Å) and per-component `jitter` bound (Å).
///
/// Jittered lattices are the standard trick for generating *valid*
/// molecular geometries: with `spacing ≥ 4.5` and `|jitter| ≤ 0.5` the
/// minimum interatomic distance stays above ~2.8 Å, so Lennard-Jones
/// clusters start without catastrophic overlaps no matter what the random
/// draw is.
pub fn lattice_geometry(n: usize, spacing: f64, jitter: f64, seed: u64) -> Vec<Vec3> {
    let mut rng = Rng::new(seed);
    let per_axis = (n as f64).cbrt().ceil() as usize;
    let mut out = Vec::with_capacity(n);
    let mut placed = 0;
    for i in 0..per_axis {
        for j in 0..per_axis {
            for k in 0..per_axis {
                if placed >= n {
                    return out;
                }
                out.push(Vec3::new(
                    2.0 + i as f64 * spacing + rng.range(-jitter, jitter),
                    2.0 + j as f64 * spacing + rng.range(-jitter, jitter),
                    2.0 + k as f64 * spacing + rng.range(-jitter, jitter),
                ));
                placed += 1;
            }
        }
    }
    out
}

/// A complete, ready-to-integrate MD system: `n ∈ [2, 12]` atoms of one
/// Lennard-Jones species on a jittered lattice in a periodic box, with
/// Maxwell–Boltzmann velocities at a random temperature.
///
/// The returned tuple is `(system, seed)`; the seed is echoed so failures
/// shrink to reproducible cases.
pub fn arb_md_system(max_atoms: usize) -> impl Strategy<Value = (System, u64)> {
    (
        2usize..=max_atoms,
        arb_lj_params(),
        arb_mass(),
        arb_temperature(),
        arb_seed(),
    )
        .prop_map(|(n, lj, mass, temperature, seed)| {
            let mut sys = System::new();
            for p in lattice_geometry(n, 4.5, 0.5, seed) {
                sys.add_atom(lj, 0.0, mass, p);
            }
            let per_axis = (n as f64).cbrt().ceil() as usize;
            sys.set_box(Box3::cubic(per_axis as f64 * 4.5 + 4.0));
            let mut rng = Rng::new(seed ^ 0xDEAD_BEEF);
            sys.init_velocities(temperature, &mut rng);
            (sys, seed)
        })
}

#[cfg(test)]
mod tests {
    use super::*;

    proptest! {
        #[test]
        fn lj_params_are_in_range((sigma, epsilon) in (2.5f64..4.0, 0.05f64..2.0)) {
            let lj = LennardJones { sigma, epsilon };
            prop_assert!(lj.sigma >= 2.5 && lj.sigma <= 4.0);
            prop_assert!(lj.epsilon >= 0.05 && lj.epsilon <= 2.0);
        }

        #[test]
        fn lattice_geometry_has_no_overlaps(n in 2usize..=12, seed in any::<u64>()) {
            let geo = lattice_geometry(n, 4.5, 0.5, seed);
            prop_assert_eq!(geo.len(), n);
            for i in 0..n {
                for j in (i + 1)..n {
                    prop_assert!(geo[i].dist(geo[j]) > 2.5,
                        "atoms {i},{j} too close: {}", geo[i].dist(geo[j]));
                }
            }
        }

        #[test]
        fn md_system_is_valid((sys, _seed) in arb_md_system(12)) {
            prop_assert!(sys.len() >= 2 && sys.len() <= 12);
            prop_assert!(sys.temperature() < 3000.0, "T = {}", sys.temperature());
            prop_assert!(sys.box_.is_some());
            // Total energy must be finite.
            prop_assert!(sys.kinetic_energy().is_finite());
        }
    }
}
