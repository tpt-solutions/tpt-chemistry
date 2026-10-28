//! End-to-end helpers that chain the `tpt-chem-io` readers into the
//! `tpt-chem-crystal` simulators (requires the `io` and `crystal` features).

use tpt_chem_crystal::bravais::Lattice;
use tpt_chem_crystal::xrd::{CrystalStructure, Reflection, XrdSimulator};
use tpt_chem_io::cif::{parse_cif, CifStructure};
use tpt_chem_io::error::IoError;

/// Build a [`Lattice`] from a parsed CIF cell (`a` along x, `b` in the xy-plane).
pub fn lattice_from_cif(cif: &CifStructure) -> Lattice {
    let [a, b, c] = cif.cell.to_cartesian_rows();
    Lattice { a, b, c }
}

/// Build the fractional-coordinate structure from a parsed CIF.
///
/// Sites are taken exactly as the file lists them: no symmetry expansion is
/// applied, so the CIF must contain every atom in the cell (for example a
/// P1 file) for the intensities to be physically meaningful.
pub fn structure_from_cif(cif: &CifStructure) -> CrystalStructure {
    CrystalStructure::new(
        cif.atoms.iter().map(|at| at.z).collect(),
        cif.atoms
            .iter()
            .map(|at| tpt_chem_core::vec3::Vec3::new(at.frac[0], at.frac[1], at.frac[2]))
            .collect(),
    )
}

/// Parse CIF text and simulate its powder pattern in one call.
///
/// # Errors
/// [`IoError`] if the CIF text cannot be parsed.
///
/// # Example
/// ```
/// # #[cfg(all(feature = "io", feature = "crystal", feature = "core"))]
/// # {
/// use tpt_chemistry::crystal::xrd::XrdSimulator;
/// use tpt_chemistry::pipeline::xrd_from_cif;
///
/// let cif = "data_nacl\n_cell_length_a 5.64\n_cell_length_b 5.64\n\
/// _cell_length_c 5.64\n_cell_angle_alpha 90\n_cell_angle_beta 90\n\
/// _cell_angle_gamma 90\nloop_\n_atom_site_label\n_atom_site_type_symbol\n\
/// _atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n\
/// Na1 Na 0 0 0\nCl1 Cl 0.5 0.5 0.5\n";
/// let peaks = xrd_from_cif(cif, &XrdSimulator::default()).unwrap();
/// assert!(!peaks.is_empty());
/// # }
/// ```
pub fn xrd_from_cif(text: &str, sim: &XrdSimulator) -> Result<Vec<Reflection>, IoError> {
    let cif = parse_cif(text)?;
    Ok(sim.simulate(&lattice_from_cif(&cif), &structure_from_cif(&cif)))
}

#[cfg(test)]
mod tests {
    use super::*;

    const QUARTZ_CIF: &str = r#"
data_quartz
_cell_length_a 4.9134
_cell_length_b 4.9134
_cell_length_c 5.4052
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 120.0
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 Si 0.4697 0.0000 0.0000
O1  O  0.4135 0.2669 0.1191
"#;

    #[test]
    fn quartz_lattice_matches_cell_volume() {
        let cif = parse_cif(QUARTZ_CIF).unwrap();
        let lat = lattice_from_cif(&cif);
        let vol = lat.a.dot(lat.b.cross(lat.c)).abs();
        assert!((vol - cif.cell.volume()).abs() < 1e-9);
    }

    #[test]
    fn quartz_strongest_lines_match_experiment() {
        let peaks = xrd_from_cif(QUARTZ_CIF, &XrdSimulator::default()).unwrap();
        let mut by_intensity = peaks.clone();
        by_intensity.sort_by(|x, y| y.intensity.total_cmp(&x.intensity));
        let top: Vec<f64> = by_intensity.iter().take(2).map(|p| p.two_theta).collect();
        assert!(top.iter().any(|t| (t - 20.9).abs() < 0.3), "{top:?}");
        assert!(top.iter().any(|t| (t - 26.7).abs() < 0.3), "{top:?}");
    }

    #[test]
    fn malformed_cif_is_an_error() {
        assert!(xrd_from_cif("not a cif", &XrdSimulator::default()).is_err());
    }
}
