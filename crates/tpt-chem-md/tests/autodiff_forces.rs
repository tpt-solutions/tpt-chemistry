//! Forces from automatic differentiation of a generic energy must match the
//! analytic all-pairs Lennard-Jones forces of `ForceModel`.

use tpt_chem_core::autodiff::{gradient, Real};
use tpt_chem_core::forcefield::LennardJones;
use tpt_chem_core::vec3::Vec3;
use tpt_chem_md::forces::ForceModel;
use tpt_chem_md::system::System;

const SIGMA: f64 = 3.4;
const EPSILON: f64 = 0.997;

/// Lennard-Jones energy of flat coordinates `[x0, y0, z0, x1, …]`, generic
/// over the number type so it can be differentiated.
fn lj_energy<T: Real>(coords: &[T]) -> T {
    let n = coords.len() / 3;
    let mut total = T::constant(0.0);
    for i in 0..n {
        for j in (i + 1)..n {
            let dx = coords[3 * i] - coords[3 * j];
            let dy = coords[3 * i + 1] - coords[3 * j + 1];
            let dz = coords[3 * i + 2] - coords[3 * j + 2];
            let r2 = dx * dx + dy * dy + dz * dz;
            let r = r2.sqrt();
            let s6 = T::constant(SIGMA).powi(6) / r.powi(6);
            total = total + T::constant(4.0 * EPSILON) * (s6 * s6 - s6);
        }
    }
    total
}

#[test]
fn autodiff_forces_match_the_analytic_kernel() {
    let positions = [
        Vec3::new(0.0, 0.0, 0.0),
        Vec3::new(3.9, 0.2, -0.1),
        Vec3::new(0.3, 4.1, 0.5),
        Vec3::new(-0.2, 0.4, 3.8),
    ];
    let flat: Vec<f64> = positions.iter().flat_map(|p| [p.x, p.y, p.z]).collect();

    let (energy, grad) = gradient(lj_energy, &flat);

    let lj = LennardJones {
        sigma: SIGMA,
        epsilon: EPSILON,
    };
    let mut sys = System::new();
    for p in positions {
        sys.add_atom(lj, 0.0, 39.95, p);
    }
    let (forces, reference_energy) = ForceModel::AllPairs { cutoff: None }.evaluate(&sys);

    assert!(
        (energy - reference_energy).abs() < 1e-10,
        "{energy} vs {reference_energy}"
    );
    for (i, f) in forces.iter().enumerate() {
        for (axis, component) in [f.x, f.y, f.z].into_iter().enumerate() {
            let autodiff_force = -grad[3 * i + axis];
            assert!(
                (autodiff_force - component).abs() < 1e-9,
                "atom {i} axis {axis}: {autodiff_force} vs {component}"
            );
        }
    }
}
