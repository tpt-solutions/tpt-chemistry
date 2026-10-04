//! Ewald vs Particle Mesh Ewald on an NaCl lattice: the two long-range
//! solvers must agree to the spline discretization level.
//!
//! Run: `cargo run -p tpt-chemistry --features md --example nacl_pme`

use tpt_chemistry::core::vec3::Vec3;
use tpt_chemistry::md::box3::Box3;
use tpt_chemistry::md::ewald::{ewald_energy_forces, EwaldParams};
use tpt_chemistry::md::pme::{pme_reciprocal_energy_forces, PmeParams};
use tpt_chemistry::md::system::System;

fn main() {
    // 3×3×3 conventional NaCl cell: 27 ions, lattice constant 5.64 Å
    // scaled to 3 cells → box 16.92 Å... keep a 2×2×2 supercell (8 ions,
    // box 11.28 Å) for a quick but well-converged check.
    let a = 5.64f64;
    let n = 2;
    let l = a * n as f64;
    let mut sys = System::new();
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let frac = [
                    (0.0, 0.0, 0.0),
                    (0.5, 0.5, 0.0),
                    (0.5, 0.0, 0.5),
                    (0.0, 0.5, 0.5),
                ];
                for (fi, fj, fk) in frac {
                    let is_na = (i + j + k + (fi * 2.0) as usize) % 2 == 0;
                    sys.add_atom(
                        tpt_chemistry::core::forcefield::LennardJones {
                            sigma: 2.0,
                            epsilon: 0.05,
                        },
                        if is_na { 1.0 } else { -1.0 },
                        if is_na { 22.99 } else { 35.45 },
                        Vec3::new(
                            (i as f64 + fi) * a,
                            (j as f64 + fj) * a,
                            (k as f64 + fk) * a,
                        ),
                    );
                }
            }
        }
    }
    let box_ = Box3::cubic(l);
    sys.set_box(box_);
    sys.wrap_in_box();
    println!("NaCl supercell: {} ions, box {l} Å", sys.len());

    let alpha = 0.35;
    let params = EwaldParams {
        alpha,
        r_cutoff: l / 2.0 - 0.01,
        g_max: 5.0,
    };
    let (f_ewald, e_ewald) = ewald_energy_forces(&sys.pos, &sys.charge, &box_, &params);
    println!("direct Ewald   : E = {e_ewald:.6} kJ/mol");

    for &mesh in &[16usize, 32, 64] {
        let pme = PmeParams {
            alpha,
            dims: [mesh, mesh, mesh],
            g_max: 0.0,
        };
        let (f_pme, _) = pme_reciprocal_energy_forces(&sys.pos, &sys.charge, &box_, &pme);
        // Isolate the reciprocal parts for an apples-to-apples comparison.
        let (_, e_recip_ewald) = {
            use tpt_chemistry::md::ewald::{ewald_real_space, ewald_self_energy};
            let (f_total, e_total) = ewald_energy_forces(&sys.pos, &sys.charge, &box_, &params);
            let (_, e_real) = ewald_real_space(&sys.pos, &sys.charge, &box_, &params);
            let e_self = ewald_self_energy(&sys.charge, alpha);
            let _ = f_total;
            (f_total, e_total - e_real - e_self)
        };
        let (_, e_pme) = pme_reciprocal_energy_forces(&sys.pos, &sys.charge, &box_, &pme);
        let max_df = f_pme
            .iter()
            .zip(f_ewald.iter())
            .map(|(a, b)| (*a - *b).norm())
            .fold(0.0f64, f64::max);
        println!(
            "PME {mesh:>2}³        : E_recip = {e_pme:10.4} vs {e_recip_ewald:10.4}  \
             (Δ = {:+.2e}, max ΔF = {:.4} kJ/mol/Å)",
            e_pme - e_recip_ewald,
            max_df
        );
    }
}
