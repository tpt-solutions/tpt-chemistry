//! Force-kernel throughput: all-pairs Lennard-Jones + Coulomb on a periodic argon lattice.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use tpt_chem_core::forcefield::LennardJones;
use tpt_chem_core::vec3::Vec3;
use tpt_chem_md::box3::Box3;
use tpt_chem_md::forces::ForceModel;
use tpt_chem_md::system::System;

/// Simple cubic argon lattice with `n` atoms per axis and `spacing` Å between sites.
fn argon_lattice(n: usize, spacing: f64) -> System {
    let lj = LennardJones {
        sigma: 3.4,
        epsilon: 0.997,
    };
    let mut sys = System::new();
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                let pos = Vec3::new(i as f64, j as f64, k as f64) * spacing;
                sys.add_atom(lj, 0.0, 39.95, pos);
            }
        }
    }
    sys.set_box(Box3::cubic(n as f64 * spacing));
    sys
}

fn all_pairs_argon(c: &mut Criterion) {
    let sys = argon_lattice(6, 3.8);
    let model = ForceModel::AllPairs { cutoff: Some(8.5) };
    c.bench_function("md/lj_all_pairs_216_atoms", |b| {
        b.iter(|| black_box(model.evaluate(black_box(&sys))))
    });
}

criterion_group!(benches, all_pairs_argon);
criterion_main!(benches);
