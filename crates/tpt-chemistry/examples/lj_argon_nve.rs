//! Lennard-Jones argon NVE dynamics: energy drift over a trajectory.
//!
//! Run: `cargo run -p tpt-chemistry --features md --example lj_argon_nve`

use tpt_chemistry::core::forcefield::LennardJones;
use tpt_chemistry::core::rng::Rng;
use tpt_chemistry::core::vec3::Vec3;
use tpt_chemistry::md::box3::Box3;
use tpt_chemistry::md::ensemble::{run_nve, run_nvt, Thermostat};
use tpt_chemistry::md::forces::ForceModel;
use tpt_chemistry::md::integrator::VelocityVerlet;
use tpt_chemistry::md::system::System;

const ARGON: LennardJones = LennardJones {
    sigma: 3.4,
    epsilon: 0.997,
};

fn main() {
    // 4×4×4 argon lattice in a periodic box, ~94 K.
    let (n_side, l) = (4, 20.0f64);
    let mut sys = System::new();
    let mut rng = Rng::new(42);
    let spacing = l / n_side as f64;
    for i in 0..n_side {
        for j in 0..n_side {
            for k in 0..n_side {
                sys.add_atom(
                    ARGON,
                    0.0,
                    39.95,
                    Vec3::new(
                        (i as f64 + 0.5) * spacing + rng.range(-0.2, 0.2),
                        (j as f64 + 0.5) * spacing + rng.range(-0.2, 0.2),
                        (k as f64 + 0.5) * spacing + rng.range(-0.2, 0.2),
                    ),
                );
            }
        }
    }
    sys.set_box(Box3::cubic(l));
    sys.wrap_in_box();
    let t_realized = sys.init_velocities(94.4, &mut rng);
    println!("argon NVE: {} atoms, T0 = {t_realized:.1} K", sys.len());

    let model = ForceModel::AllPairs { cutoff: Some(9.0) };
    let vv = VelocityVerlet::new(1.0);

    // Relax with a thermostat, then switch it off.
    let nvt = run_nvt(
        &mut sys,
        &model,
        vv,
        2_000,
        94.4,
        Thermostat::Berendsen { tau: 100.0 },
        2_000,
    );
    println!(
        "equilibration (NVT, 2 ps): mean T = {:.1} K",
        nvt.mean_temperature
    );

    let report = run_nve(&mut sys, &model, vv, 5_000, 500);
    println!("production (NVE, 5 ps):");
    println!(
        "  {:>8}  {:>10}  {:>10}  {:>10}",
        "t (fs)", "E_pot", "E_kin", "E_total"
    );
    for s in &report.snapshots {
        println!(
            "  {t:8}  {p:10.3}  {k:10.3}  {tot:10.3}",
            t = s.time_fs,
            p = s.potential,
            k = s.kinetic,
            tot = s.potential + s.kinetic
        );
    }
    let e_first = {
        let s = &report.snapshots[0];
        s.potential + s.kinetic
    };
    let e_last = report
        .snapshots
        .last()
        .map(|s| s.potential + s.kinetic)
        .unwrap();
    println!(
        "\nrelative energy drift over 5 ps: {:.2e}",
        (e_last - e_first).abs() / e_first.abs()
    );
}
