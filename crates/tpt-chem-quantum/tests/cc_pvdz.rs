//! cc-pVDZ with spherical d functions: basis size, rotational invariance, and
//! the energy of H₂ and water.

use tpt_chem_core::molecule::Molecule;
use tpt_chem_core::vec3::Vec3;
use tpt_chem_quantum::basis::BasisSet;
use tpt_chem_quantum::hf::rhf_with_basis;

fn water(rotation: Option<[[f64; 3]; 3]>) -> Molecule {
    // Experimental geometry: r(OH) = 0.9572 Å, angle 104.52°.
    let r = 0.9572;
    let half = (104.52f64 / 2.0).to_radians();
    let raw = [
        Vec3::ZERO,
        Vec3::new(r * half.sin(), 0.0, r * half.cos()),
        Vec3::new(-r * half.sin(), 0.0, r * half.cos()),
    ];
    let z = [8u8, 1, 1];
    let mut mol = Molecule::new("water");
    for (&zi, p) in z.iter().zip(raw) {
        let q = match rotation {
            Some(m) => Vec3::new(
                m[0][0] * p.x + m[0][1] * p.y + m[0][2] * p.z,
                m[1][0] * p.x + m[1][1] * p.y + m[1][2] * p.z,
                m[2][0] * p.x + m[2][1] * p.y + m[2][2] * p.z,
            ),
            None => p,
        };
        mol.add_atom_raw(zi, q);
    }
    mol
}

/// Rotation by `angle` about the unit axis `(1, 2, 3)/√14` (Rodrigues).
fn rotation(angle: f64) -> [[f64; 3]; 3] {
    let n = 14f64.sqrt();
    let (x, y, z) = (1.0 / n, 2.0 / n, 3.0 / n);
    let (c, s) = (angle.cos(), angle.sin());
    let t = 1.0 - c;
    [
        [t * x * x + c, t * x * y - s * z, t * x * z + s * y],
        [t * x * y + s * z, t * y * y + c, t * y * z - s * x],
        [t * x * z - s * y, t * y * z + s * x, t * z * z + c],
    ]
}

#[test]
fn water_has_24_functions() {
    // O: 3s 2p 1d (3 + 6 + 5 = 14); each H: 2s 1p (2 + 3 = 5).
    let res = rhf_with_basis(&water(None), BasisSet::CcPvdz).unwrap();
    assert_eq!(res.n_basis, 24);
}

#[test]
fn water_energy_is_rotation_invariant() {
    // A wrong set of d harmonics would make the energy depend on orientation.
    // One rotation keeps this test quick in debug builds (about 25 s each).
    let reference = rhf_with_basis(&water(None), BasisSet::CcPvdz)
        .unwrap()
        .energy;
    let rotated = water(Some(rotation(0.7)));
    let e = rhf_with_basis(&rotated, BasisSet::CcPvdz).unwrap().energy;
    assert!((e - reference).abs() < 1e-8, "{e} vs {reference}");
}

#[test]
fn water_energy_regression_value() {
    // Recorded from this implementation, not from a reference program. The
    // literature RHF/cc-pVDZ value for water near this geometry is about
    // −76.0268 Eh; this test only locks the current result.
    let e = rhf_with_basis(&water(None), BasisSet::CcPvdz)
        .unwrap()
        .energy;
    assert!((e - (-76.0267986975)).abs() < 1e-8, "{e:.10}");
}

#[test]
fn gradients_refuse_spherical_bases() {
    use tpt_chem_quantum::gradient::rhf_gradient_with_basis;
    use tpt_chem_quantum::hf::HfError;
    let err = rhf_gradient_with_basis(&water(None), BasisSet::CcPvdz).unwrap_err();
    assert_eq!(err, HfError::SphericalUnsupported);
}
