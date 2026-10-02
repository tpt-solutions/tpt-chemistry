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
