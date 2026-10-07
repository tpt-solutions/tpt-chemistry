//! Rigid TIP3P water: 64 molecules, SHAKE/RATTLE constraints, PME
//! electrostatics. Velocity-rescale equilibration to 300 K, then NVE.
//!
//! Run: `cargo run --release -p tpt-chemistry --features core,md --example water_pme`

use tpt_chemistry::core::rng::Rng;
use tpt_chemistry::md::forces::ForceModel;
use tpt_chemistry::md::pme::PmeParams;
use tpt_chemistry::md::water::{constrained_temperature, water_box, TIP3P};

fn main() {
    let mut wb = water_box(&TIP3P, 4);
    let l = wb.system.box_.expect("periodic box").length.x;
    println!("{} molecules, box {:.3} A", wb.n_molecules, l);
    let model = ForceModel::LjPlusPme {
        lj_cutoff: 6.0,
        pme: PmeParams {
            alpha: 0.4,
            dims: [16, 16, 16],
            g_max: 0.0,
        },
    };
    wb.system.init_velocities(300.0, &mut Rng::new(2024));
    wb.constraints
        .rattle_velocities(&mut wb.system)
        .expect("RATTLE");
    let dt = 1.0;
    println!("step  T(K)   PE(kJ/mol)  E_total(kJ/mol)");
    for step in 0..=1500 {
        let pe = wb
            .constraints
            .step(&mut wb.system, &model, dt)
            .expect("constraints converge");
        // Equilibrate: rescale to 300 K every 10 steps for the first 1000.
        if step < 1000 && step % 10 == 0 {
            let t = constrained_temperature(&wb);
            let s = (300.0 / t).sqrt();
            for v in &mut wb.system.vel {
                *v = *v * s;
            }
        }
        if step % 250 == 0 {
            println!(
                "{step:5} {:6.1} {pe:12.2} {:12.2}",
                constrained_temperature(&wb),
                pe + wb.system.kinetic_energy()
            );
        }
    }
}
