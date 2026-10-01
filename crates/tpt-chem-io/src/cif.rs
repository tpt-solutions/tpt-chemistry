//! CIF (Crystallographic Information File) format.
//!
//! Parses cell parameters, the Hermann–Mauguin space-group symbol, and the
//! `loop_ _atom_site…` block (fractional coordinates). CIF structures live
//! in fractional space, so this module carries its own small
//! [`CellParams`] type for fractional→Cartesian conversion; `tpt-chem-crystal`
//! provides the full lattice machinery for everything beyond that.


use tpt_chem_core::element;
use tpt_chem_core::vec3::Vec3;

use crate::error::{bad_number, parse_err, IoError, Result};

const FORMAT: &str = "cif";

/// Unit-cell parameters: edge lengths in Å, angles in degrees.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CellParams {
    /// a edge (Å).
    pub a: f64,
    /// b edge (Å).
    pub b: f64,
    /// c edge (Å).
    pub c: f64,
    /// α angle (degrees), between b and c.
    pub alpha: f64,
    /// β angle (degrees), between a and c.
    pub beta: f64,
    /// γ angle (degrees), between a and b.
    pub gamma: f64,
}

impl CellParams {
    /// Cubic convenience constructor.
    pub fn cubic(a: f64) -> Self {
        CellParams {
            a,
            b: a,
            c: a,
            alpha: 90.0,
            beta: 90.0,
            gamma: 90.0,
        }
    }

    /// Unit-cell volume in Å³:
    /// `abc·√(1 − cos²α − cos²β − cos²γ + 2cosα cosβ cosγ)`.
    pub fn volume(&self) -> f64 {
        let (ca, cb, cg) = (
            self.alpha.to_radians().cos(),
            self.beta.to_radians().cos(),
            self.gamma.to_radians().cos(),
        );
        let s = 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg;
        self.a * self.b * self.c * s.max(0.0).sqrt()
    }

    /// Cartesian lattice vectors (rows), standard orientation: `a` along x,
    /// `b` in the xy-plane.
    pub fn to_cartesian_rows(&self) -> [Vec3; 3] {
        let (al, be, ga) = (
            self.alpha.to_radians(),
            self.beta.to_radians(),
            self.gamma.to_radians(),
        );
        let a_vec = Vec3::new(self.a, 0.0, 0.0);
        let b_vec = Vec3::new(self.b * ga.cos(), self.b * ga.sin(), 0.0);
        let c_x = self.c * be.cos();
        let c_y = self.c * (al.cos() - be.cos() * ga.cos()) / ga.sin();
        let c_z2 = self.c * self.c - c_x * c_x - c_y * c_y;
        let c_vec = Vec3::new(c_x, c_y, c_z2.max(0.0).sqrt());
        [a_vec, b_vec, c_vec]
    }

    /// Convert fractional coordinates to Cartesian (Å).
    pub fn frac_to_cart(&self, f: [f64; 3]) -> Vec3 {
        let [va, vb, vc] = self.to_cartesian_rows();
        va * f[0] + vb * f[1] + vc * f[2]
    }
}

/// One site from the `_atom_site` loop.
#[derive(Clone, Debug, PartialEq)]
pub struct CifAtom {
    /// Site label, e.g. `"Si1"`.
    pub label: String,
    /// Atomic number (from `_atom_site_type_symbol`).
    pub z: u8,
    /// Fractional coordinates.
    pub frac: [f64; 3],
    /// Site occupancy (defaults to 1.0 when absent).
    pub occupancy: f64,
}

/// A parsed CIF block: cell, symmetry symbol, and atom sites.
#[derive(Clone, Debug, PartialEq)]
pub struct CifStructure {
    /// Data-block name (text after `data_`).
    pub name: String,
    /// Unit-cell parameters.
    pub cell: CellParams,
    /// Hermann–Mauguin symbol, e.g. `"F d -3 m"` (empty if absent).
    pub space_group: String,
    /// Atom sites in fractional coordinates.
    pub atoms: Vec<CifAtom>,
}

/// Strip a trailing standard uncertainty: `5.431(2)` → `5.431`.
fn strip_esd(token: &str) -> &str {
    match token.find('(') {
        Some(i) => &token[..i],
        None => token,
    }
}

fn unquote(token: &str) -> String {
    let t = token.trim();
    if (t.starts_with('\'') && t.ends_with('\'') && t.len() >= 2)
        || (t.starts_with('"') && t.ends_with('"') && t.len() >= 2)
    {
        t[1..t.len() - 1].to_string()
    } else {
        t.to_string()
    }
}

/// Parse the first `data_` block of a CIF file.
///
/// # Errors
/// [`IoError`] if the block, cell parameters, or the `_atom_site` loop are
/// missing or malformed.
pub fn parse_cif(text: &str) -> Result<CifStructure> {
    let mut name = String::new();
    let mut cell = [0.0f64; 6]; // a b c alpha beta gamma
    let mut space_group = String::new();
    let mut in_loop = false;
    let mut columns: Vec<String> = Vec::new();
    let mut atoms: Vec<CifAtom> = Vec::new();

    for (i, raw) in text.lines().enumerate() {
        let line_no = i + 1;
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some(block) = line.strip_prefix("data_") {
            if !name.is_empty() {
                break; // only the first block
            }
            name = block.trim().to_string();
            continue;
        }
        if line == "loop_" {
            in_loop = true;
            columns.clear();
            continue;
        }
        if in_loop && line.starts_with('_') {
            columns.push(line.to_string());
            continue;
        }
        if in_loop && !columns.is_empty() {
            // A data row of the loop (or a non-loop tag ends it).
            if line.starts_with('_') {
                in_loop = false;
            } else {
                parse_atom_row(line, line_no, &columns, &mut atoms)?;
                continue;
            }
        }
        // Regular tag value pairs.
        let mut tokens = line.splitn(2, char::is_whitespace);
        let tag = tokens.next().unwrap_or("");
        let value = tokens.next().map(str::trim).unwrap_or("");
        match tag {
            "_cell_length_a" => cell[0] = parse_num(value, line_no)?,
            "_cell_length_b" => cell[1] = parse_num(value, line_no)?,
            "_cell_length_c" => cell[2] = parse_num(value, line_no)?,
            "_cell_angle_alpha" => cell[3] = parse_num(value, line_no)?,
            "_cell_angle_beta" => cell[4] = parse_num(value, line_no)?,
            "_cell_angle_gamma" => cell[5] = parse_num(value, line_no)?,
            "_symmetry_space_group_name_H-M" | "_space_group_name_H-M_alt" => {
                space_group = unquote(value);
            }
            "_atom_site_fract_x" | "_atom_site_type_symbol" | "_atom_site_label" => {
                // Tag with value on the same line inside a loop — ignored.
            }
            _ => {}
        }
    }

    if name.is_empty() {
        return Err(IoError::MissingSection {
            format: FORMAT,
            section: "data_".into(),
        });
    }
    if cell.iter().any(|&v| v == 0.0) {
        return Err(IoError::MissingSection {
            format: FORMAT,
            section: "_cell_*".into(),
        });
    }

    Ok(CifStructure {
        name,
        cell: CellParams {
            a: cell[0],
            b: cell[1],
            c: cell[2],
            alpha: cell[3],
            beta: cell[4],
            gamma: cell[5],
        },
        space_group,
        atoms,
    })
}

fn parse_num(value: &str, line_no: usize) -> Result<f64> {
    let unquoted = unquote(value);
    let cleaned = strip_esd(unquoted.trim());
    cleaned
        .parse::<f64>()
        .map_err(|_| bad_number(FORMAT, line_no, value))
}

fn parse_atom_row(
    line: &str,
    line_no: usize,
    columns: &[String],
    atoms: &mut Vec<CifAtom>,
) -> Result<()> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    if tokens.len() < columns.len() {
        return Err(parse_err(
            FORMAT,
            line_no,
            format!("atom_site row has {} fields, expected {}", tokens.len(), columns.len()),
        ));
    }
    let col = |name: &str| columns.iter().position(|c| c == name);
    let (Some(i_label), Some(i_sym), Some(ix), Some(iy), Some(iz)) = (
        col("_atom_site_label"),
        col("_atom_site_type_symbol"),
        col("_atom_site_fract_x"),
        col("_atom_site_fract_y"),
        col("_atom_site_fract_z"),
    ) else {
        return Ok(()); // loop without the columns we need — skip rows
    };
    let sym = unquote(tokens[i_sym]);
    let z = element::from_symbol(&sym).ok_or_else(|| IoError::UnknownElement {
        format: FORMAT,
        symbol: sym,
        line: Some(line_no),
    })?;
    let frac = [
        parse_num(strip_esd(tokens[ix]), line_no)?,
        parse_num(strip_esd(tokens[iy]), line_no)?,
        parse_num(strip_esd(tokens[iz]), line_no)?,
    ];
    let occupancy = match col("_atom_site_occupancy") {
        Some(io) => parse_num(strip_esd(tokens[io]), line_no)?,
        None => 1.0,
    };
    atoms.push(CifAtom {
        label: unquote(tokens[i_label]),
        z,
        frac,
        occupancy,
    });
    Ok(())
}

/// Serialize a [`CifStructure`] back to CIF text.
pub fn write_cif(s: &CifStructure) -> String {
    let mut out = String::new();
    out.push_str(&format!("data_{}\n", s.name));
    out.push_str(&format!("_cell_length_a {}\n", fmt_num(s.cell.a)));
    out.push_str(&format!("_cell_length_b {}\n", fmt_num(s.cell.b)));
    out.push_str(&format!("_cell_length_c {}\n", fmt_num(s.cell.c)));
    out.push_str(&format!("_cell_angle_alpha {}\n", fmt_num(s.cell.alpha)));
    out.push_str(&format!("_cell_angle_beta {}\n", fmt_num(s.cell.beta)));
    out.push_str(&format!("_cell_angle_gamma {}\n", fmt_num(s.cell.gamma)));
    if !s.space_group.is_empty() {
        out.push_str(&format!(
            "_symmetry_space_group_name_H-M '{}'\n",
            s.space_group
        ));
    }
    out.push_str("loop_\n");
    out.push_str("_atom_site_label\n");
    out.push_str("_atom_site_type_symbol\n");
    out.push_str("_atom_site_fract_x\n");
    out.push_str("_atom_site_fract_y\n");
    out.push_str("_atom_site_fract_z\n");
    out.push_str("_atom_site_occupancy\n");
    for a in &s.atoms {
        let sym = element::symbol(a.z).unwrap_or("X");
        out.push_str(&format!(
            "{} {} {} {} {} {}\n",
            a.label,
            sym,
            fmt_num(a.frac[0]),
            fmt_num(a.frac[1]),
            fmt_num(a.frac[2]),
            fmt_num(a.occupancy)
        ));
    }
    out
}

fn fmt_num(v: f64) -> String {
    format!("{v:.6}")
}

#[cfg(test)]
mod tests {
    use super::*;

    const NACL: &str = "\
# generated by tpt-chem-io test
data_sodium_chloride
_cell_length_a 5.6402(3)
_cell_length_b 5.6402
_cell_length_c 5.6402
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
_symmetry_space_group_name_H-M 'F m -3 m'
loop_
 _atom_site_label
 _atom_site_type_symbol
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
 Na1 Na 0.000000 0.000000 0.000000 1.0
 Cl1 Cl 0.500000 0.500000 0.500000 1.0
";

    #[test]
    fn parse_nacl() {
        let s = parse_cif(NACL).unwrap();
        assert_eq!(s.name, "sodium_chloride");
        assert_eq!(s.space_group, "F m -3 m");
        assert_eq!(s.atoms.len(), 2);
        // ESD stripped.
        assert!((s.cell.a - 5.6402).abs() < 1e-9);
        assert_eq!(s.atoms[0].z, 11);
        assert_eq!(s.atoms[1].z, 17);
    }

    #[test]
    fn fractional_to_cartesian_nacl() {
        let s = parse_cif(NACL).unwrap();
        let na = s.cell.frac_to_cart(s.atoms[0].frac);
        let cl = s.cell.frac_to_cart(s.atoms[1].frac);
        // Nearest Na–Cl distance is a√3/2 for the (½,½,½) site.
        assert!((na.dist(cl) - 5.6402 * 3.0f64.sqrt() / 2.0).abs() < 1e-9);
        assert!((s.cell.volume() - 5.6402f64.powi(3)).abs() < 1e-6);
    }

    #[test]
    fn triclinic_volume_and_vectors() {
        let cell = CellParams {
            a: 5.0,
            b: 6.0,
            c: 7.0,
            alpha: 80.0,
            beta: 95.0,
            gamma: 105.0,
        };
        let rows = cell.to_cartesian_rows();
        // |a| = a, a·b = ab cos γ, a·c = ac cos β.
        assert!((rows[0].norm() - 5.0).abs() < 1e-12);
        assert!((rows[0].dot(rows[1]) / (5.0 * 6.0) - 105f64.to_radians().cos()).abs() < 1e-12);
        assert!((rows[0].dot(rows[2]) / (5.0 * 7.0) - 95f64.to_radians().cos()).abs() < 1e-12);
        // The triple product reproduces the analytic volume.
        let vol = rows[0].dot(rows[1].cross(rows[2]));
        assert!((vol - cell.volume()).abs() < 1e-9);
    }

    #[test]
    fn roundtrip() {
        let s = parse_cif(NACL).unwrap();
        let text = write_cif(&s);
        let again = parse_cif(&text).unwrap();
        assert_eq!(again.name, "sodium_chloride");
        assert_eq!(again.space_group, "F m -3 m");
        assert_eq!(again.atoms.len(), 2);
        assert_eq!(again.atoms[1].frac, [0.5, 0.5, 0.5]);
    }

    #[test]
    fn missing_cell_is_error() {
        assert!(matches!(
            parse_cif("data_x\nloop_\n_atom_site_label\n"),
            Err(IoError::MissingSection { .. })
        ));
    }
}
