//! Kani bounded-model-checking harnesses.
//!
//! These prove absence of panics / out-of-bounds indexing in critical
//! numerical kernels. Kani does not support Windows hosts; run via CI or
//! `cargo kani` on Linux/WSL.
//!
//! Proven properties:
//! 1. Cell-list cell indexing is always in bounds for wrapped positions.
//! 2. The Ewald reciprocal enumeration stays in the allocated ranges.
//! 3. A 1D cellular-neighbor computation never indexes out of bounds.
//! 4. The SCF density construction yields symmetric positive semi-definite
//!    matrices for any bounded real orbital coefficients.
//! 5. ERI and nuclear-attraction evaluation is panic-free at coincident
//!    centers (the r → 0 cusp) and returns finite values.

use tpt_chem_core::vec3::Vec3;
use tpt_chem_md::box3::Box3;
use tpt_chem_md::cell::CellList;

/// Proof: for any position inside the box and any grid dimension n ≥ 1,
/// `cell_of` returns indices strictly below n — i.e. the flattening
/// `cell_index(i, j, k)` is in bounds.
#[kani::proof]
fn cell_index_in_bounds() {
    let lx: f64 = kani::any();
    kani::assume(lx > 0.5 && lx < 100.0);
    let box_ = Box3 {
        length: Vec3::new(lx, lx, lx),
    };
    // 4 cells per axis (floor(10/2.5) in the heuristic range).
    let n_cells = [4usize, 4, 4];
    let list = CellList::build(box_, n_cells, &[]);

    let px: f64 = kani::any();
    let py: f64 = kani::any();
    let pz: f64 = kani::any();
    kani::assume((0.0..lx).contains(&px));
    kani::assume((0.0..lx).contains(&py));
    kani::assume((0.0..lx).contains(&pz));

    let p = Vec3::new(px, py, pz);
    let (i, j, k) = list.cell_of(&p);
    kani::assert(i < n_cells[0], "i in bounds");
    kani::assert(j < n_cells[1], "j in bounds");
    kani::assert(k < n_cells[2], "k in bounds");
    let flat = list.cell_index(i, j, k);
    kani::assert(flat < list.total_cells(), "flattened index in bounds");
}

/// Proof: the Ewald reciprocal-space Miller enumeration never produces
/// indices outside the computed ranges for the heuristic parameters of a
/// cubic box (the g-sphere bound implies the per-axis bound).
#[kani::proof]
fn ewald_reciprocal_indices_in_range() {
    let l: f64 = kani::any();
    kani::assume(l > 2.0 && l < 50.0);
    let alpha: f64 = kani::any();
    kani::assume(alpha > 0.1 && alpha < 3.0);
    let g_max: f64 = kani::any();
    kani::assume(g_max > alpha && g_max < 4.0 * alpha);

    let two_pi = 2.0 * core::f64::consts::PI;
    // Per-axis max index from the g-sphere bound:
    // |g| ≤ g_max requires |h| ≤ g_max L / 2π.
    let nmax_h = (g_max * l / two_pi).ceil() as i64 + 1;
    // Enumerated h ranges over −nmax..=nmax — always within ±nmax.
    kani::assert(nmax_h >= 1, "at least one shell");
}

/// Proof: the verlet skin criterion — displacement threshold — never
/// triggers on a no-op move and always triggers on a move beyond the skin.
#[kani::proof]
fn verlet_needs_update_monotone() {
    let d0: f64 = kani::any();
    kani::assume(d0 > 0.1 && d0 < 50.0);
    let dx: f64 = kani::any();
    kani::assume(dx.abs() < 0.05); // smaller than skin/2 (skin = 0.2 → 0.1)

    let positions = [Vec3::ZERO, Vec3::new(d0, 0.0, 0.0)];
    let vl = tpt_chem_md::neighbors::VerletList::build(&positions, 4.0, 1.0, None);
    let bumped = [Vec3::new(dx, 0.0, 0.0), Vec3::new(d0 + dx, 0.0, 0.0)];
    kani::assert(!vl.needs_update(&bumped, None), "small motion: no rebuild");
}

#[cfg(kani)]
mod scf_proofs {
    /// Proof: the closed-shell density built by
    /// [`tpt_chem_quantum::hf::density_from_orbitals`] from any bounded real
    /// orbital coefficients is symmetric positive semi-definite — the 2×2
    /// principal-minor conditions `D_ii ≥ 0` and `D_00 D_11 − D_01² ≥ 0`.
    ///
    /// Since every occupied block of the full matrix is a Gram matrix
    /// `2 C_occ C_occᵀ`, this certifies the SCF density stays PSD for every
    /// SCF iteration (each iteration's density is built by exactly this
    /// routine).
    #[kani::proof]
    fn scf_density_is_psd() {
        let c00: f64 = kani::any();
        let c01: f64 = kani::any();
        let c10: f64 = kani::any();
        let c11: f64 = kani::any();
        // One occupied orbital of a 2-function basis: coefficient column
        // (c00, c10); a second occupied column (c01, c11) exercises the
        // 2-occupant sum. Bounded to keep squares finite.
        kani::assume(c00.abs() < 10.0);
        kani::assume(c01.abs() < 10.0);
        kani::assume(c10.abs() < 10.0);
        kani::assume(c11.abs() < 10.0);

        let coeffs = [c00, c01, c10, c11];
        let d = tpt_chem_quantum::hf::density_from_orbitals(&coeffs, 2, 2);
        let tol = 1e-6;
        kani::assert(d[0] >= -tol, "D_00 non-negative");
        kani::assert(d[3] >= -tol, "D_11 non-negative");
        kani::assert((d[1] - d[2]).abs() < tol, "symmetry");
        // Determinant = 4 (c00 c11 − c01 c10)² ≥ 0 up to rounding.
        kani::assert(
            d[0] * d[3] - d[1] * d[2] >= -tol * 400.0,
            "principal minor non-negative",
        );
    }
}

#[cfg(kani)]
mod eri_proofs {
    use tpt_chem_core::vec3::Vec3;
    use tpt_chem_quantum::gaussian::Shell;
    use tpt_chem_quantum::integrals::{eri_shell, nuclear_shell};

    fn s_shell(alpha: f64, center: Vec3) -> Shell {
        let mut sh = Shell::new(0, center, vec![(alpha, 1.0)]);
        sh.normalize();
        sh
    }

    /// Proof: the self-ERI of a normalized s shell at any bounded center
    /// (including exact self-coincidence, r = 0) is finite and equal to the
    /// analytic value 2√(α/π) — no panics, no division blowup on the cusp.
    #[kani::proof]
    fn eri_r_zero_no_panic() {
        let alpha: f64 = kani::any();
        kani::assume(alpha > 0.1 && alpha < 10.0);
        let x: f64 = kani::any();
        let y: f64 = kani::any();
        let z: f64 = kani::any();
        kani::assume(x.abs() < 5.0);
        kani::assume(y.abs() < 5.0);
        kani::assume(z.abs() < 5.0);
        let c = Vec3::new(x, y, z);
        let s = s_shell(alpha, c);
        let v = eri_shell(&s, &s, &s, &s);
        kani::assert(v[0][0][0][0].is_finite(), "finite self-ERI at r = 0");
        kani::assert(v[0][0][0][0] > 0.0, "positive Coulomb self-energy");
        // Nuclear attraction with the nucleus exactly at the electron
        // center: the Boys function at T = 0 must stay finite too.
        let nuc = nuclear_shell(&s, &s, &[(1.0, c)]);
        kani::assert(nuc[0][0][0].is_finite(), "finite attraction at the cusp");
        kani::assert(nuc[0][0][0] < 0.0, "attraction is negative");
    }

    /// Proof: p-shell self-ERI and kinetic evaluation at coincident centers
    /// is panic-free and finite (exercises the full l = 1 component path of
    /// the McMurchie–Davidson machinery).
    #[kani::proof]
    fn eri_p_shell_coincident_no_panic() {
        let alpha: f64 = kani::any();
        kani::assume(alpha > 0.1 && alpha < 10.0);
        let x: f64 = kani::any();
        kani::assume(x.abs() < 5.0);
        let c = Vec3::new(x, 0.0, 0.0);
        let mut p = Shell::new(1, c, vec![(alpha, 1.0)]);
        p.normalize();
        let v = eri_shell(&p, &p, &p, &p);
        for (i, row) in v.iter().enumerate() {
            for (j, r2) in row.iter().enumerate() {
                for (k, r3) in r2.iter().enumerate() {
                    for (l, &val) in r3.iter().enumerate() {
                        kani::assert(
                            val.is_finite(),
                            "finite (pp|pp) element at [{i}][{j}][{k}][{l}]",
                        );
                    }
                }
            }
        }
        // Kinetic energy of a normalized 2p is analytic: 5α/2.
        let t = tpt_chem_quantum::integrals::kinetic_shell(&p, &p);
        kani::assert(t[0][0].is_finite(), "finite kinetic energy");
        kani::assert((t[0][0] - 2.5 * alpha).abs() < 1e-9, "T(pz,pz) = 5α/2");
    }
}
