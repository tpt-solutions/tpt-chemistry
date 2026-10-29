//! Gaussian input and output: the subset needed to hand a geometry and run
//! settings to Gaussian, and to read the results back.
//!
//! **Input:** route section (`#` lines), title, charge and multiplicity, and
//! Cartesian coordinates in Å. `%` Link 0 lines before the route are skipped.
//!
//! **Output:** the last `SCF Done` energy (Hartree), the last
//! `Standard orientation` geometry (Å), and whether `Normal termination`
//! appears. Every other section is ignored.
//!
//! The layouts follow Gaussian's documented text output. The test fixtures
//! are hand-written in that layout, not copied from a Gaussian run.
//!
//! # Examples
//!
//! ```
//! use tpt_chem_io::gaussian::{parse_gaussian_input, write_gaussian_input, GaussianInput};
//! use tpt_chem_core::molecule::Molecule;
//! use tpt_chem_core::vec3::Vec3;
//!
//! let mut molecule = Molecule::new("water");
//! molecule.add_atom_raw(8, Vec3::new(0.0, 0.0, 0.117));
//! molecule.add_atom_raw(1, Vec3::new(0.0, 0.757, -0.469));
//! molecule.add_atom_raw(1, Vec3::new(0.0, -0.757, -0.469));
//! let input = GaussianInput {
//!     route: "#p HF/STO-3G SP".into(),
//!     title: "water".into(),
//!     charge: 0,
//!     multiplicity: 1,
//!     molecule,
//! };
//! let text = write_gaussian_input(&input).unwrap();
//! let again = parse_gaussian_input(&text).unwrap();
//! assert_eq!(again.molecule.formula(), "H2O");
//! ```

use tpt_chem_core::element;
use tpt_chem_core::molecule::Molecule;
use tpt_chem_core::vec3::Vec3;

use crate::error::{bad_number, parse_err, IoError, Result};

const FORMAT: &str = "Gaussian";

/// A Gaussian job input: route, title, charge, multiplicity and geometry.
#[derive(Clone, Debug)]
pub struct GaussianInput {
    /// Route section, starting with `#`. Multi-line routes are kept with `\n`.
    pub route: String,
    /// Title section (free text, may span lines).
    pub title: String,
    /// Total molecular charge written on the charge/multiplicity line.
    pub charge: i32,
    /// Spin multiplicity (2S + 1).
    pub multiplicity: u32,
    /// Geometry in Å. The molecule's own formal charge is not used.
    pub molecule: Molecule,
}

/// Results read from a Gaussian output file.
#[derive(Clone, Debug, Default)]
pub struct GaussianOutput {
    /// Energy from the last `SCF Done` line (Hartree).
    pub energy: Option<f64>,
    /// Geometry from the last `Standard orientation` block (Å).
    pub geometry: Option<Molecule>,
    /// Whether a `Normal termination` line was seen.
    pub normal_termination: bool,
}

/// Render a job input as Gaussian text.
///
/// # Errors
/// [`IoError::Unsupported`] when an atomic number has no element symbol.
pub fn write_gaussian_input(input: &GaussianInput) -> Result<String> {
    let mut out = String::new();
    out.push_str(&input.route);
    out.push_str("\n\n");
    out.push_str(&input.title);
    out.push_str("\n\n");
    out.push_str(&format!("{} {}\n", input.charge, input.multiplicity));
    for atom in input.molecule.atoms() {
        let symbol = element::symbol(atom.z).ok_or_else(|| IoError::Unsupported {
            format: FORMAT,
            message: format!("atomic number {} has no element symbol", atom.z),
        })?;
        let p = atom.pos;
        out.push_str(&format!(
            "{:<2} {:>14.8} {:>14.8} {:>14.8}\n",
            symbol, p.x, p.y, p.z
        ));
    }
    out.push('\n');
    Ok(out)
}

/// Read a job input written in Gaussian's layout.
///
/// # Errors
/// [`IoError`] when the route, charge/multiplicity line, or a coordinate
/// line is missing or malformed, or an element symbol is unknown.
pub fn parse_gaussian_input(text: &str) -> Result<GaussianInput> {
    let lines: Vec<&str> = text.lines().collect();
    let mut i = 0;

    // Link 0 lines and blanks before the route are skipped.
    while i < lines.len() && !lines[i].trim_start().starts_with('#') {
        let line = lines[i].trim();
        if !line.is_empty() && !line.starts_with('%') {
            return Err(parse_err(
                FORMAT,
                i + 1,
                "expected a route line starting with '#'",
            ));
        }
        i += 1;
    }
    let route = take_block(&lines, &mut i).join("\n");
    if route.is_empty() {
        return Err(parse_err(FORMAT, i + 1, "missing route section"));
    }

    let title = take_block(&lines, &mut i).join("\n");

    let charge_line = i + 1;
    let charge_text = lines
        .get(i)
        .copied()
        .ok_or_else(|| parse_err(FORMAT, charge_line, "missing charge/multiplicity line"))?;
    let mut fields = charge_text.split_whitespace();
    let charge_token = fields.next().unwrap_or("");
    let mult_token = fields.next().unwrap_or("");
    let charge = charge_token
        .parse::<i32>()
        .map_err(|_| bad_number(FORMAT, charge_line, charge_token))?;
    let multiplicity = mult_token
        .parse::<u32>()
        .map_err(|_| bad_number(FORMAT, charge_line, mult_token))?;
    i += 1;

    let mut molecule = Molecule::new(title.lines().next().unwrap_or(""));
    for (offset, line) in take_block(&lines, &mut i).iter().enumerate() {
        let line_no = charge_line + 1 + offset;
        let mut tokens = line.split_whitespace();
        let symbol_token = tokens
            .next()
            .ok_or_else(|| parse_err(FORMAT, line_no, "empty coordinate line"))?;
        let symbol: String = symbol_token
            .chars()
            .take_while(|c| c.is_ascii_alphabetic())
            .collect();
        let z = element::from_symbol(&symbol).ok_or_else(|| IoError::UnknownElement {
            format: FORMAT,
            symbol: symbol.clone(),
            line: Some(line_no),
        })?;
        let mut coord = [0.0f64; 3];
        for c in &mut coord {
            let token = tokens
                .next()
                .ok_or_else(|| parse_err(FORMAT, line_no, "expected x y z after the symbol"))?;
            *c = token
                .parse::<f64>()
                .map_err(|_| bad_number(FORMAT, line_no, token))?;
        }
        molecule.add_atom_raw(z, Vec3::new(coord[0], coord[1], coord[2]));
    }
    molecule.set_formal_charge(charge);

    Ok(GaussianInput {
        route,
        title,
        charge,
        multiplicity,
        molecule,
    })
}

/// Collect lines up to the next blank line, advancing `i` past the blank.
fn take_block<'a>(lines: &[&'a str], i: &mut usize) -> Vec<&'a str> {
    let mut block = Vec::new();
    while *i < lines.len() && !lines[*i].trim().is_empty() {
        block.push(lines[*i].trim());
        *i += 1;
    }
    // Skip the blank line that ends the block.
    *i += 1;
    block
}

/// Read energies, geometry and termination from Gaussian output text.
///
/// # Errors
/// [`IoError`] when an `SCF Done` energy or an orientation row has a
/// malformed number.
pub fn parse_gaussian_output(text: &str) -> Result<GaussianOutput> {
    let lines: Vec<&str> = text.lines().collect();
    let mut out = GaussianOutput::default();
    for (i, line) in lines.iter().enumerate() {
        if let Some(pos) = line.find("SCF Done:") {
            // "SCF Done:  E(RHF) =  -76.0267  A.U. after ..."
            let after = &line[pos..];
            if let Some(eq) = after.find('=') {
                let token = after[eq + 1..].split_whitespace().next().unwrap_or("");
                let value = token
                    .parse::<f64>()
                    .map_err(|_| bad_number(FORMAT, i + 1, token))?;
                out.energy = Some(value);
            }
        }
        if line.contains("Normal termination of Gaussian") {
            out.normal_termination = true;
        }
        if line.contains("Standard orientation:") {
            out.geometry = Some(read_orientation(&lines, i + 1)?);
        }
    }
    Ok(out)
}

/// Read the table that follows a `Standard orientation:` header.
///
/// Data rows have six fields: centre number, atomic number, atomic type,
/// and x, y, z. The table ends at the dashed rule after the rows.
fn read_orientation(lines: &[&str], start: usize) -> Result<Molecule> {
    let mut molecule = Molecule::new("");
    let mut rows = 0usize;
    for (offset, line) in lines[start..].iter().enumerate() {
        let trimmed = line.trim();
        if trimmed.starts_with("---") {
            if rows > 0 {
                break;
            }
            continue;
        }
        let fields: Vec<&str> = trimmed.split_whitespace().collect();
        let is_row = fields.len() == 6 && fields[0].parse::<u32>().is_ok();
        if !is_row {
            continue;
        }
        let line_no = start + offset + 1;
        let z = fields[1]
            .parse::<u8>()
            .map_err(|_| bad_number(FORMAT, line_no, fields[1]))?;
        let mut coord = [0.0f64; 3];
        for (c, token) in coord.iter_mut().zip(&fields[3..6]) {
            *c = token
                .parse::<f64>()
                .map_err(|_| bad_number(FORMAT, line_no, token))?;
        }
        molecule.add_atom_raw(z, Vec3::new(coord[0], coord[1], coord[2]));
        rows += 1;
    }
    Ok(molecule)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn water() -> Molecule {
        let mut mol = Molecule::new("water");
        mol.add_atom_raw(8, Vec3::new(0.0, 0.0, 0.117));
        mol.add_atom_raw(1, Vec3::new(0.0, 0.757, -0.469));
        mol.add_atom_raw(1, Vec3::new(0.0, -0.757, -0.469));
        mol
    }

    #[test]
    fn writes_gaussian_input_layout() {
        let input = GaussianInput {
            route: "#p HF/STO-3G Opt".into(),
            title: "water".into(),
            charge: 0,
            multiplicity: 1,
            molecule: water(),
        };
        let text = write_gaussian_input(&input).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "#p HF/STO-3G Opt");
        assert_eq!(lines[1], "");
        assert_eq!(lines[2], "water");
        assert_eq!(lines[3], "");
        assert_eq!(lines[4], "0 1");
        assert_eq!(lines[5], "O      0.00000000     0.00000000     0.11700000");
        assert_eq!(lines[6], "H      0.00000000     0.75700000    -0.46900000");
    }

    #[test]
    fn input_round_trips() {
        let input = GaussianInput {
            route: "#p B3LYP/6-31G(d) SP\nscf=tight".into(),
            title: "water cation".into(),
            charge: 1,
            multiplicity: 2,
            molecule: water(),
        };
        let text = write_gaussian_input(&input).unwrap();
        let again = parse_gaussian_input(&text).unwrap();
        assert_eq!(again.route, input.route);
        assert_eq!(again.title, input.title);
        assert_eq!(again.charge, 1);
        assert_eq!(again.multiplicity, 2);
        assert_eq!(again.molecule.formula(), "H2O");
        assert_eq!(again.molecule.formal_charge(), 1);
        for i in 0..3 {
            let a = input.molecule.get_atom(i).unwrap().pos;
            let b = again.molecule.get_atom(i).unwrap().pos;
            assert!(a.dist(b) < 1e-7, "atom {i}: {a:?} vs {b:?}");
        }
    }

    #[test]
    fn input_skips_link0_and_reads_bare_symbols() {
        let text = "%nprocshared=4\n%mem=1GB\n#p HF/STO-3G SP\n\nH2\n\n0 1\nH1 0 0 0\nH 0 0 0.74 extra\n\n";
        let input = parse_gaussian_input(text).unwrap();
        assert_eq!(input.route, "#p HF/STO-3G SP");
        assert_eq!(input.molecule.formula(), "H2");
    }

    #[test]
    fn malformed_input_is_rejected() {
        assert!(parse_gaussian_input("").is_err());
        assert!(parse_gaussian_input("#p HF/STO-3G SP\n\nt\n\nx 1\n").is_err());
        assert!(parse_gaussian_input("#p HF/STO-3G SP\n\nt\n\n0 1\nXx 0 0 0\n\n").is_err());
        assert!(parse_gaussian_input("#p HF/STO-3G SP\n\nt\n\n0 1\nO 0 0\n\n").is_err());
    }

    const OUTPUT: &str = "\
 SCF Done:  E(RHF) =  -75.9800000000     A.U. after    5 cycles
                          Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          8           0        0.000000    0.000000    0.100000
      2          1           0        0.000000    0.750000   -0.450000
      3          1           0        0.000000   -0.750000   -0.450000
 ---------------------------------------------------------------------
 SCF Done:  E(RHF) =  -75.9851234567     A.U. after    7 cycles
                          Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          8           0        0.000000    0.000000    0.117000
      2          1           0        0.000000    0.757000   -0.469000
      3          1           0        0.000000   -0.757000   -0.469000
 ---------------------------------------------------------------------
 Normal termination of Gaussian 16 at Sat Oct 10 12:00:00 2026.
";

    #[test]
    fn output_takes_last_energy_and_geometry() {
        let out = parse_gaussian_output(OUTPUT).unwrap();
        assert!((out.energy.unwrap() + 75.9851234567).abs() < 1e-12);
        assert!(out.normal_termination);
        let geom = out.geometry.unwrap();
        assert_eq!(geom.formula(), "H2O");
        let o = geom.get_atom(0).unwrap().pos;
        assert!((o.z - 0.117).abs() < 1e-9);
        let h = geom.get_atom(1).unwrap().pos;
        assert!((h.y - 0.757).abs() < 1e-9);
    }

    #[test]
    fn output_without_termination_reports_it() {
        let truncated: String = OUTPUT
            .lines()
            .filter(|l| !l.contains("Normal termination"))
            .collect::<Vec<_>>()
            .join("\n");
        let out = parse_gaussian_output(&truncated).unwrap();
        assert!(!out.normal_termination);
        assert!(out.energy.is_some());
    }

    #[test]
    fn empty_output_yields_nothing() {
        let out = parse_gaussian_output("").unwrap();
        assert!(out.energy.is_none());
        assert!(out.geometry.is_none());
        assert!(!out.normal_termination);
    }
}
