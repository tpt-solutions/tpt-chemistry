//! Python bindings for tpt-chemistry.
//!
//! ```python
//! import tpt_chem
//! tpt_chem.rhf_energy_xyz(open("water.xyz").read())   # Hartree
//! tpt_chem.canonical_smiles("OC(C)=O")                # "CC(=O)O"
//! ```

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use tpt_chem_core::units::BOHR_ANGSTROM;
use tpt_chem_io::{parse_smiles, parse_xyz, write_smiles};
use tpt_chem_quantum::gradient::rhf_gradient;
use tpt_chem_quantum::hf::rhf;

/// RHF/STO-3G total energy (Hartree) for a geometry in XYZ format (Å).
#[pyfunction]
fn rhf_energy_xyz(xyz: &str) -> PyResult<f64> {
    let mol = parse_xyz(xyz).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let result = rhf(&mol).map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(result.energy)
}

/// RHF/STO-3G energy (Hartree) and forces (Hartree/Å) for a geometry in XYZ
/// format. Returns `(energy, [(fx, fy, fz), …])`, one force per atom.
#[pyfunction]
fn rhf_energy_forces_xyz(xyz: &str) -> PyResult<(f64, Vec<(f64, f64, f64)>)> {
    let mol = parse_xyz(xyz).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let (result, gradient) = rhf_gradient(&mol).map_err(|e| PyValueError::new_err(e.to_string()))?;
    // Gradients are Hartree/Bohr; forces are −∇E in Hartree/Å.
    let forces = gradient
        .iter()
        .map(|g| {
            (
                -g.x / BOHR_ANGSTROM,
                -g.y / BOHR_ANGSTROM,
                -g.z / BOHR_ANGSTROM,
            )
        })
        .collect();
    Ok((result.energy, forces))
}

/// Molecular formula of a SMILES string, such as "C2H4O2".
#[pyfunction]
fn formula_from_smiles(smiles: &str) -> PyResult<String> {
    let mol = parse_smiles(smiles).map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(mol.formula())
}

/// Canonical SMILES of a SMILES string, so equivalent spellings compare equal.
#[pyfunction]
fn canonical_smiles(smiles: &str) -> PyResult<String> {
    let mol = parse_smiles(smiles).map_err(|e| PyValueError::new_err(e.to_string()))?;
    write_smiles(&mol).map_err(|e| PyValueError::new_err(e.to_string()))
}

#[pymodule]
fn tpt_chem(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(rhf_energy_xyz, m)?)?;
    m.add_function(wrap_pyfunction!(rhf_energy_forces_xyz, m)?)?;
    m.add_function(wrap_pyfunction!(formula_from_smiles, m)?)?;
    m.add_function(wrap_pyfunction!(canonical_smiles, m)?)?;
    Ok(())
}
