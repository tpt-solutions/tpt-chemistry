//! Property tests for boundary and corner cases: cell indexing, the Verlet
//! rebuild criterion, SCF density positivity, and integral evaluation at
//! coincident centres. They run under `cargo test --workspace`, so every PR
//! exercises them.
//!
//! This is randomized evidence, not an exhaustive proof over all `f64`:
//! inputs are drawn from ordinary ranges plus the corner values in
//! [`edge_f64`](crate::strategies::edge_f64).

use proptest::prelude::*;

use tpt_chem_core::vec3::Vec3;
use tpt_chem_md::box3::Box3;
use tpt_chem_md::cell::CellList;
use tpt_chem_md::neighbors::VerletList;
use tpt_chem_quantum::gaussian::Shell;
use tpt_chem_quantum::hf::density_from_orbitals;
use tpt_chem_quantum::integrals::{eri_shell, kinetic_shell, nuclear_shell};

use crate::strategies::{edge_f64, next_down};

/// A coordinate in `[0, l]` that favours the cell boundaries (0, l and the
/// largest `f64` below l), alongside interior points.
fn coord_in_box(l: f64) -> impl Strategy<Value = f64> {
    prop_oneof![
        8 => (0.0f64..1.0).prop_map(move |f| f * l),
        1 => Just(0.0),
        1 => Just(l),
        1 => Just(next_down(l)),
    ]
}

/// A position inside a cubic box of side `l`, boundaries included.
fn position_in_box(l: f64) -> impl Strategy<Value = Vec3> {
    (coord_in_box(l), coord_in_box(l), coord_in_box(l)).prop_map(|(x, y, z)| Vec3::new(x, y, z))
}

/// A unit-length-ish direction; degenerate near-zero vectors are rejected.
fn direction() -> impl Strategy<Value = Vec3> {
    (-1.0f64..1.0, -1.0f64..1.0, -1.0f64..1.0)
        .prop_filter("non-degenerate direction", |(x, y, z)| {
            x * x + y * y + z * z > 1e-3
        })
        .prop_map(|(x, y, z)| Vec3::new(x, y, z).normalize())
}

/// A displacement magnitude clear of the half-skin threshold (skin = 1 Å),
/// so the expected rebuild decision is unambiguous under rounding.
fn clear_of_half_skin() -> impl Strategy<Value = f64> {
    prop_oneof![0.0f64..0.49, 0.51f64..3.0]
}

/// Coefficients: ordinary values plus the edge corners, clamped to ±10 so
/// the density stays finite.
fn coefficient() -> impl Strategy<Value = f64> {
    prop_oneof![
        -10.0f64..10.0,
        edge_f64().prop_map(|c| c.clamp(-10.0, 10.0)),
    ]
}

/// A nuclear-centre-sized coordinate triple, clamped to ±5 Bohr, drawn from
/// the edge corners as well as the bulk.
fn centre() -> impl Strategy<Value = Vec3> {
    (edge_f64(), edge_f64(), edge_f64())
        .prop_map(|(x, y, z)| Vec3::new(x.clamp(-5.0, 5.0), y.clamp(-5.0, 5.0), z.clamp(-5.0, 5.0)))
}

fn s_shell(alpha: f64, center: Vec3) -> Shell {
    let mut sh = Shell::new(0, center, vec![(alpha, 1.0)]);
    sh.normalize();
    sh
}

proptest! {
    /// `cell_index_in_bounds`: for any grid of 1–8 cells per
    /// axis and any position in the box (boundaries included), the cell
    /// indices are in range and the flattened index is in bounds.
    #[test]
    fn cell_index_in_bounds(
        (l, n, p) in (0.5f64..100.0, 1usize..=8, 1usize..=8, 1usize..=8)
            .prop_flat_map(|(l, nx, ny, nz)| {
                (Just(l), Just([nx, ny, nz]), position_in_box(l))
            })
    ) {
        let list = CellList::build(Box3::cubic(l), n, &[]);
        let (i, j, k) = list.cell_of(&p);
        prop_assert!(i < n[0], "i = {i} out of {}", n[0]);
        prop_assert!(j < n[1], "j = {j} out of {}", n[1]);
        prop_assert!(k < n[2], "k = {k} out of {}", n[2]);
        prop_assert!(list.cell_index(i, j, k) < list.total_cells());
    }

    /// Verlet rebuild criterion, open boundaries: a move
    /// of at most skin/2 does not force a rebuild; a move beyond it does.
    #[test]
    fn verlet_open_rebuilds_exactly_past_half_skin(
        dir in direction(),
        mag in clear_of_half_skin(),
    ) {
        let reference = [Vec3::ZERO, Vec3::new(1.5, 0.0, 0.0)];
        let vl = VerletList::build(&reference, 4.0, 1.0, None);
        let moved = [dir.scale_to(mag), Vec3::new(1.5, 0.0, 0.0)];
        prop_assert_eq!(vl.needs_update(&moved, None), mag > 0.5);
    }

    /// The same criterion under periodic boundaries, where displacement is
    /// measured by minimum image. Positions are wrapped into the box, so a
    /// small move across the boundary must not trigger a rebuild.
    #[test]
    fn verlet_periodic_rebuilds_exactly_past_half_skin(
        dir in direction(),
        mag in clear_of_half_skin(),
    ) {
        let l = 20.0;
        let reference = [Vec3::new(0.1, 0.1, 0.1), Vec3::new(3.0, 1.0, 1.0)];
        let box_ = Some(Box3::cubic(l));
        let vl = VerletList::build(&reference, 4.0, 1.0, box_);
        let raw = reference[0] + dir.scale_to(mag);
        let wrapped = Vec3::new(
            raw.x.rem_euclid(l),
            raw.y.rem_euclid(l),
            raw.z.rem_euclid(l),
        );
        let moved = [wrapped, reference[1]];
        prop_assert_eq!(vl.needs_update(&moved, box_), mag > 0.5);
    }

    /// SCF density positivity, for 1–5 basis functions
    /// and any occupation count. The density is symmetric and its quadratic
    /// form `xᵀDx` is non-negative for random `x`, so it is positive
    /// semi-definite.
    #[test]
    fn density_from_orbitals_is_psd(
        (nb, n_occ, coeffs, x) in (1usize..=5).prop_flat_map(|nb| {
            (
                Just(nb),
                1usize..=nb,
                prop::collection::vec(coefficient(), nb * nb),
                prop::collection::vec(-1.0f64..1.0, nb),
            )
        })
    ) {
        let d = density_from_orbitals(&coeffs, n_occ, nb);
        let scale = d.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let tol = 1e-9 * (1.0 + scale);

        for i in 0..nb {
            prop_assert!(d[i * nb + i] >= -tol, "D[{i}][{i}] = {}", d[i * nb + i]);
            for j in 0..nb {
                prop_assert!((d[i * nb + j] - d[j * nb + i]).abs() <= tol, "asymmetric");
            }
        }

        let mut quad = 0.0;
        for i in 0..nb {
            for j in 0..nb {
                quad += x[i] * d[i * nb + j] * x[j];
            }
        }
        prop_assert!(quad >= -tol * nb as f64 * nb as f64, "xᵀDx = {quad}");
    }

    /// Coincident-centre s-shell integrals: the self-ERI and nuclear
    /// attraction of a normalized s shell are finite, the Coulomb
    /// self-energy is positive, and attraction is negative, including with
    /// the nucleus exactly on the electron centre (r → 0).
    #[test]
    fn s_shell_integrals_finite_at_coincident_centre(
        alpha in 0.1f64..10.0,
        c in centre(),
    ) {
        let s = s_shell(alpha, c);
        let v = eri_shell(&s, &s, &s, &s)[0][0][0][0];
        prop_assert!(v.is_finite(), "self-ERI = {v}");
        prop_assert!(v > 0.0, "self-ERI = {v}");

        let nuc = nuclear_shell(&s, &s, &[(1.0, c)])[0][0];
        prop_assert!(nuc.is_finite(), "attraction = {nuc}");
        prop_assert!(nuc < 0.0, "attraction = {nuc}");
    }

    /// Near-zero separation between shell centres (down to the smallest
    /// positive `f64`): integrals stay finite and keep their signs.
    #[test]
    fn s_shell_integrals_finite_at_near_zero_separation(
        alpha in 0.1f64..10.0,
        c in centre(),
        delta in prop_oneof![Just(0.0), Just(f64::MIN_POSITIVE), 1e-12f64..1e-6],
    ) {
        let a = s_shell(alpha, c);
        let b = s_shell(alpha, Vec3::new(c.x + delta, c.y, c.z));
        let v = eri_shell(&a, &a, &b, &b)[0][0][0][0];
        prop_assert!(v.is_finite() && v > 0.0, "(aa|bb) = {v}");

        let nuc = nuclear_shell(&a, &b, &[(1.0, c)])[0][0];
        prop_assert!(nuc.is_finite() && nuc < 0.0, "attraction = {nuc}");
    }

    /// Coincident-centre p-shell integrals: every (pp|pp)
    /// element at coincident centres is finite, and the kinetic diagonal of
    /// a normalized p shell is the analytic 5α/2.
    #[test]
    fn p_shell_coincident_integrals_finite_and_kinetic_analytic(
        alpha in 0.1f64..10.0,
        c in centre(),
    ) {
        let mut p = Shell::new(1, c, vec![(alpha, 1.0)]);
        p.normalize();

        let v = eri_shell(&p, &p, &p, &p);
        for r1 in &v {
            for r2 in r1 {
                for r3 in r2 {
                    for &val in r3 {
                        prop_assert!(val.is_finite(), "(pp|pp) element = {val}");
                    }
                }
            }
        }

        let t = kinetic_shell(&p, &p);
        for (k, row) in t.iter().enumerate() {
            prop_assert!(row[k].is_finite(), "T[{k}][{k}] = {}", row[k]);
            prop_assert!(
                (row[k] - 2.5 * alpha).abs() < 1e-9,
                "T[{k}][{k}] = {}, expected {}",
                row[k],
                2.5 * alpha
            );
        }
    }
}
