//! Energies against OpenMM, an independent MD code. The fixtures in
//! `tests/data/` are written by `tools/openmm_reference.py`.

use serde_json::Value;

use tpt_chem_core::forcefield::LennardJones;
use tpt_chem_core::vec3::Vec3;
use tpt_chem_md::box3::Box3;
use tpt_chem_md::ewald::{ewald_energy_forces, EwaldParams};
use tpt_chem_md::forces::ForceModel;
use tpt_chem_md::system::System;

fn load(name: &str) -> Value {
    let path = format!("{}/tests/data/{name}.json", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    serde_json::from_str(&text).unwrap()
}

fn vec3s(v: &Value) -> Vec<Vec3> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|p| {
            let a = p.as_array().unwrap();
            Vec3::new(
                a[0].as_f64().unwrap(),
                a[1].as_f64().unwrap(),
                a[2].as_f64().unwrap(),
            )
        })
        .collect()
}

fn floats(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}

#[test]
fn periodic_coulomb_matches_openmm_pme() {
    let fx = load("periodic_ions");
    let box_len = fx["box_angstrom"].as_f64().unwrap();
    let positions = vec3s(&fx["positions_angstrom"]);
    let charges = floats(&fx["charges"]);
    let box_ = Box3::cubic(box_len);
    let params = EwaldParams::heuristic(&box_);
    let (_, ours) = ewald_energy_forces(&positions, &charges, &box_, &params);
    let reference = fx["openmm_energy_kj_mol"].as_f64().unwrap();
    // OpenMM's PME is set to 1e-8 error tolerance; ours is converged to the same order.
    assert!(
        (ours - reference).abs() < 1e-3,
        "Ewald {ours:.6} vs OpenMM PME {reference:.6} kJ/mol"
    );
}

#[test]
fn isolated_lj_and_coulomb_match_openmm_no_cutoff() {
    let fx = load("isolated_mixed");
    let positions = vec3s(&fx["positions_angstrom"]);
    let charges = floats(&fx["charges"]);
    let sigma = floats(&fx["sigma_angstrom"]);
    let epsilon = floats(&fx["epsilon_kj_mol"]);
    let mut sys = System::new();
    for i in 0..positions.len() {
        let lj = LennardJones {
            sigma: sigma[i],
            epsilon: epsilon[i],
        };
        sys.add_atom(lj, charges[i], 39.95, positions[i]);
    }
    let (_, ours) = ForceModel::AllPairs { cutoff: None }.evaluate(&sys);
    let reference = fx["openmm_energy_kj_mol"].as_f64().unwrap();
    assert!(
        (ours - reference).abs() < 1e-6 * reference.abs().max(1.0),
        "all-pairs {ours:.8} vs OpenMM {reference:.8} kJ/mol"
    );
}
