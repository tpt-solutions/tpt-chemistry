//! Supercell and slab builders.

use tpt_chem_core::vec3::Vec3;

use crate::bravais::Lattice;
use crate::xrd::CrystalStructure;

/// Replicate a cell `n = [na, nb, nc]` times along its lattice vectors.
///
/// Fractional coordinates are rescaled to the new cell; site order is
/// image-major (all sites of image `(0,0,0)`, then the next image, ...).
///
/// # Panics
///
/// Panics if any repeat count is zero.
pub fn supercell(
    lat: &Lattice,
    structure: &CrystalStructure,
    n: [usize; 3],
) -> (Lattice, CrystalStructure) {
    assert!(
        n.iter().all(|&k| k > 0),
        "supercell: repeat counts must be >= 1"
    );
    let big = Lattice {
        a: lat.a * n[0] as f64,
        b: lat.b * n[1] as f64,
        c: lat.c * n[2] as f64,
    };
    let mut z = Vec::new();
    let mut pos = Vec::new();
    for i in 0..n[0] {
        for j in 0..n[1] {
            for k in 0..n[2] {
                for (zi, p) in structure
                    .atomic_numbers
                    .iter()
                    .zip(&structure.frac_positions)
                {
                    z.push(*zi);
                    pos.push(Vec3::new(
                        (p.x + i as f64) / n[0] as f64,
                        (p.y + j as f64) / n[1] as f64,
                        (p.z + k as f64) / n[2] as f64,
                    ));
                }
            }
        }
    }
    let mut out = CrystalStructure::new(z, pos);
    out.b_factor = structure.b_factor;
    (big, out)
}

/// Build a slab: `layers` unit cells stacked along `c`, plus `vacuum` Å of
/// empty space appended along `c`. No re-orientation is performed, so `c`
/// should already be the surface normal.
///
/// # Panics
///
/// Panics if `layers` is zero or `vacuum` is negative.
pub fn slab(
    lat: &Lattice,
    structure: &CrystalStructure,
    layers: usize,
    vacuum: f64,
) -> (Lattice, CrystalStructure) {
    assert!(vacuum >= 0.0, "slab: vacuum must be non-negative");
    let (stacked, s) = supercell(lat, structure, [1, 1, layers]);
    let len = stacked.c.norm();
    let new_len = len + vacuum;
    let c = stacked.c * (new_len / len);
    let scale = len / new_len;
    let pos = s
        .frac_positions
        .iter()
        .map(|p| Vec3::new(p.x, p.y, p.z * scale))
        .collect();
    let mut out = CrystalStructure::new(s.atomic_numbers, pos);
    out.b_factor = s.b_factor;
    (
        Lattice {
            a: stacked.a,
            b: stacked.b,
            c,
        },
        out,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nacl_like() -> (Lattice, CrystalStructure) {
        (
            Lattice::cubic(4.0),
            CrystalStructure::new(
                vec![11, 17],
                vec![Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.5, 0.5, 0.5)],
            ),
        )
    }

    #[test]
    fn supercell_counts_and_volume() {
        let (l, s) = nacl_like();
        let (bl, bs) = supercell(&l, &s, [2, 3, 4]);
        assert_eq!(bs.atomic_numbers.len(), 2 * 24);
        assert!((bl.volume() - l.volume() * 24.0).abs() < 1e-9);
        assert!(bs.frac_positions.iter().all(|p| (0.0..1.0).contains(&p.x)
            && (0.0..1.0).contains(&p.y)
            && (0.0..1.0).contains(&p.z)));
    }

    #[test]
    fn supercell_preserves_cartesian_images() {
        let (l, s) = nacl_like();
        let (bl, bs) = supercell(&l, &s, [2, 1, 1]);
        let cart: Vec<Vec3> = bs
            .frac_positions
            .iter()
            .map(|p| bl.to_cartesian(*p))
            .collect();
        // Second image of site 0 sits one lattice vector along a.
        assert!(cart[2].dist(l.a) < 1e-12);
    }

    #[test]
    fn slab_adds_vacuum_and_keeps_atoms_fixed() {
        let (l, s) = nacl_like();
        let (sl, ss) = slab(&l, &s, 2, 10.0);
        assert!((sl.c.norm() - 18.0).abs() < 1e-12);
        // Atom at layer 1, z = 0.5 cell → 6 Å Cartesian.
        let z_cart = sl.to_cartesian(ss.frac_positions[3]).z;
        assert!((z_cart - 6.0).abs() < 1e-12);
    }
}
