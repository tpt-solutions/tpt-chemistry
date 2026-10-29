//! The AMBER-style model against OpenMM's amber14 evaluation of the same
//! N-methylacetamide system. The fixture `tests/data/nma_amber.json` is written
//! by `tools/amber_reference.py`.

use serde_json::Value;

use tpt_chem_core::vec3::Vec3;
use tpt_chem_md::amber::{AmberModel, AngleTerm, BondTerm, TorsionTerm};

fn load() -> Value {
    let path = format!("{}/tests/data/nma_amber.json", env!("CARGO_MANIFEST_DIR"));
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    serde_json::from_str(&text).unwrap()
}

fn model(fx: &Value) -> AmberModel {
    let atoms = fx["atoms"].as_array().unwrap();
    let f = |v: &Value, key: &str| v[key].as_f64().unwrap();
    let u = |v: &Value, key: &str| v[key].as_u64().unwrap() as usize;
    AmberModel {
        charges: atoms.iter().map(|a| f(a, "charge")).collect(),
        sigma: atoms.iter().map(|a| f(a, "sigma")).collect(),
        epsilon: atoms.iter().map(|a| f(a, "epsilon")).collect(),
        bonds: fx["bonds"]
            .as_array()
            .unwrap()
            .iter()
            .map(|b| BondTerm {
                i: u(b, "i"),
                j: u(b, "j"),
                r0: f(b, "r0"),
                k: f(b, "k"),
            })
            .collect(),
        angles: fx["angles"]
            .as_array()
            .unwrap()
            .iter()
            .map(|a| AngleTerm {
                i: u(a, "i"),
                j: u(a, "j"),
                k: u(a, "k"),
                theta0: f(a, "theta0"),
                k_theta: f(a, "kk"),
            })
            .collect(),
        rb_torsions: Vec::new(),
        torsions: fx["torsions"]
            .as_array()
            .unwrap()
            .iter()
            .map(|t| TorsionTerm {
                atoms: [u(t, "i"), u(t, "j"), u(t, "k"), u(t, "l")],
                periodicity: t["periodicity"].as_u64().unwrap() as u32,
                phase: f(t, "phase"),
                k: f(t, "kt"),
            })
            .collect(),
        scee: f(fx, "scee"),
        scnb: f(fx, "scnb"),
    }
}

fn positions(fx: &Value) -> Vec<Vec3> {
    fx["positions_angstrom"]
        .as_array()
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

#[test]
fn energy_matches_openmm_amber14() {
    let fx = load();
    let (energy, _) = model(&fx).energy_forces(&positions(&fx));
    let reference = fx["openmm_energy_kj_mol"].as_f64().unwrap();
    assert!(
        (energy - reference).abs() < 1e-3,
        "energy {energy:.6} vs OpenMM {reference:.6} kJ/mol"
    );
}

#[test]
fn forces_match_openmm_amber14() {
    let fx = load();
    let (_, forces) = model(&fx).energy_forces(&positions(&fx));
    let reference = fx["openmm_forces_kj_mol_angstrom"].as_array().unwrap();
    for (i, (f, r)) in forces.iter().zip(reference).enumerate() {
        let r = r.as_array().unwrap();
        let expected = [
            r[0].as_f64().unwrap(),
            r[1].as_f64().unwrap(),
            r[2].as_f64().unwrap(),
        ];
        for (axis, (got, want)) in [f.x, f.y, f.z].into_iter().zip(expected).enumerate() {
            assert!(
                (got - want).abs() < 1e-3,
                "atom {i} axis {axis}: {got:.6} vs OpenMM {want:.6} kJ/mol/Å"
            );
        }
    }
}
