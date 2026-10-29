//! End-to-end RHF/STO-3G energy for water (integrals, SCF and Fock build).

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use tpt_chem_core::molecule::Molecule;
use tpt_chem_core::vec3::Vec3;
use tpt_chem_quantum::hf::rhf;

fn water() -> Molecule {
    let mut mol = Molecule::new("water");
    mol.add_atom_raw(8, Vec3::new(0.0, 0.0, 0.117));
    mol.add_atom_raw(1, Vec3::new(0.0, 0.757, -0.469));
    mol.add_atom_raw(1, Vec3::new(0.0, -0.757, -0.469));
    mol
}

fn rhf_water(c: &mut Criterion) {
    let mol = water();
    let mut group = c.benchmark_group("quantum");
    group.sample_size(10);
    group.bench_function("rhf_sto3g_water", |b| {
        b.iter(|| black_box(rhf(black_box(&mol)).unwrap().energy))
    });
    group.finish();
}

criterion_group!(benches, rhf_water);
criterion_main!(benches);
