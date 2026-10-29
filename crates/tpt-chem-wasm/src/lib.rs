//! WebAssembly bindings for tpt-chemistry.
//!
//! Runs an RHF/STO-3G energy and SMILES helpers in a browser or Node with no
//! server. Build with `wasm-pack build crates/tpt-chem-wasm --target web`.

use tpt_chem_io::{parse_smiles, parse_xyz, write_smiles};
use tpt_chem_quantum::hf::rhf;
use wasm_bindgen::prelude::*;

/// RHF/STO-3G total energy (Hartree) for a geometry in XYZ format (coordinates in Å).
///
/// # Errors
/// Returns a JavaScript error when the XYZ text is malformed or the SCF fails.
#[wasm_bindgen]
pub fn rhf_energy_xyz(xyz: &str) -> Result<f64, JsError> {
    let mol = parse_xyz(xyz).map_err(|e| JsError::new(&e.to_string()))?;
    let result = rhf(&mol).map_err(|e| JsError::new(&e.to_string()))?;
    Ok(result.energy)
}

/// Molecular formula of a SMILES string, such as `C2H4O2` for acetic acid.
///
/// # Errors
/// Returns a JavaScript error when the SMILES string is malformed.
#[wasm_bindgen]
pub fn formula_from_smiles(smiles: &str) -> Result<String, JsError> {
    let mol = parse_smiles(smiles).map_err(|e| JsError::new(&e.to_string()))?;
    Ok(mol.formula())
}

/// Canonical SMILES of a SMILES string, so equivalent spellings compare equal.
///
/// # Errors
/// Returns a JavaScript error when the SMILES string is malformed or the
/// molecule carries a net charge, which the writer does not support.
#[wasm_bindgen]
pub fn canonical_smiles(smiles: &str) -> Result<String, JsError> {
    let mol = parse_smiles(smiles).map_err(|e| JsError::new(&e.to_string()))?;
    write_smiles(&mol).map_err(|e| JsError::new(&e.to_string()))
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER: &str = "3\nwater\nO 0 0 0.117\nH 0 0.757 -0.469\nH 0 -0.757 -0.469\n";

    #[test]
    fn water_energy_matches_the_native_solver() {
        let mol = parse_xyz(WATER).unwrap();
        let native = rhf(&mol).unwrap().energy;
        let wasm = rhf_energy_xyz(WATER).unwrap();
        assert!((wasm - native).abs() < 1e-12);
    }

    #[test]
    fn smiles_helpers_on_the_native_side() {
        let mol = parse_smiles("OC(C)=O").unwrap();
        assert_eq!(mol.formula(), "C2H4O2");
        assert_eq!(write_smiles(&mol).unwrap(), "CC(=O)O");
    }
}
