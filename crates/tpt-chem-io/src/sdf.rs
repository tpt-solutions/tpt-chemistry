//! MDL SDF (V2000 connection table) format.
//!
//! An SDF file is one or more MOL records, each closed by a `$$$$` line.
//! A record is a three-line header (name, program/timestamp, comment), a
//! counts line (`aaabbb…  V2000`), `aaa` atom lines, `bbb` bond lines, an
//! optional property block (only `M  CHG` is read), and `M  END`.
//! Coordinates are Å. Bond orders 1, 2 and 3 map to single, double and
//! triple; 4 maps to aromatic. Other orders fall back to single, as the
//! MOL2 reader does.
//!
//! The writer does not emit an `M  CHG` line, so the total formal charge is
//! read but not written back out.

use tpt_chem_core::element;
use tpt_chem_core::molecule::{BondOrder, Molecule};
use tpt_chem_core::vec3::Vec3;

use crate::error::{bad_number, parse_err, parse_f64, IoError, Result};

const FORMAT: &str = "sdf";

/// Parse the first MOL record in `text`.
///
/// # Errors
/// [`IoError`] on malformed records, unknown elements or out-of-range bonds.
pub fn parse_sdf(text: &str) -> Result<Molecule> {
    parse_sdf_all(text)?
        .into_iter()
        .next()
        .ok_or(IoError::UnexpectedEof { format: FORMAT })
}

/// Parse every record in `text`.
///
/// # Errors
/// [`IoError`] on malformed records, unknown elements or out-of-range bonds.
pub fn parse_sdf_all(text: &str) -> Result<Vec<Molecule>> {
    let lines: Vec<&str> = text.lines().collect();
    let mut mols = Vec::new();
    let mut start = 0;
    for (i, line) in lines.iter().enumerate() {
        if line.trim_end() == "$$$$" {
            mols.push(parse_record(&lines[start..i], start + 1)?);
            start = i + 1;
        }
    }
    // A final record may omit the `$$$$` terminator.
    let rest = &lines[start..];
    if rest.iter().any(|l| !l.trim().is_empty()) {
        mols.push(parse_record(rest, start + 1)?);
    }
    Ok(mols)
}

/// Serialize one molecule as a single SDF record, terminated by `$$$$`.
pub fn write_sdf(mol: &Molecule) -> String {
    let mut out = String::new();
    out.push_str(&format!("{}\n\n\n", mol.name()));
    out.push_str(&format!(
        "{:>3}{:>3}  0  0  0  0  0  0  0  0999 V2000\n",
        mol.len(),
        mol.n_bonds()
    ));
    for atom in mol.atoms() {
        let sym = element::symbol(atom.z).unwrap_or("X");
        out.push_str(&format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0  0  0  0  0  0  0  0  0  0  0\n",
            atom.pos.x, atom.pos.y, atom.pos.z, sym
        ));
    }
    for bond in mol.bonds() {
        let order = match bond.order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 4,
        };
        out.push_str(&format!(
            "{:>3}{:>3}{:>3}  0\n",
            bond.a.index() + 1,
            bond.b.index() + 1,
            order
        ));
    }
    out.push_str("M  END\n$$$$\n");
    out
}

/// Parse one record. `first_line` is the 1-based line number of `lines[0]`.
fn parse_record(lines: &[&str], first_line: usize) -> Result<Molecule> {
    if lines.len() < 4 {
        return Err(IoError::UnexpectedEof { format: FORMAT });
    }
    let name = lines[0].trim();
    let counts_no = first_line + 3;
    let counts: Vec<&str> = lines[3].split_whitespace().collect();
    if counts.len() < 2 {
        return Err(parse_err(
            FORMAT,
            counts_no,
            "counts line needs atom and bond counts",
        ));
    }
    let n_atoms: usize = counts[0]
        .parse()
        .map_err(|_| bad_number(FORMAT, counts_no, counts[0]))?;
    let n_bonds: usize = counts[1]
        .parse()
        .map_err(|_| bad_number(FORMAT, counts_no, counts[1]))?;

    let mut mol = Molecule::new(name);
    let atoms_start = 4;
    let bonds_start = atoms_start + n_atoms;
    if lines.len() < bonds_start + n_bonds {
        return Err(IoError::UnexpectedEof { format: FORMAT });
    }

    for i in 0..n_atoms {
        let line_no = first_line + atoms_start + i;
        let t: Vec<&str> = lines[atoms_start + i].split_whitespace().collect();
        if t.len() < 4 {
            return Err(parse_err(FORMAT, line_no, "atom line needs x y z symbol"));
        }
        let x = parse_f64(FORMAT, line_no, t[0])?;
        let y = parse_f64(FORMAT, line_no, t[1])?;
        let z = parse_f64(FORMAT, line_no, t[2])?;
        let elem = element::from_symbol(t[3]).ok_or_else(|| IoError::UnknownElement {
            format: FORMAT,
            symbol: t[3].to_string(),
            line: Some(line_no),
        })?;
        mol.add_atom_raw(elem, Vec3::new(x, y, z));
    }

    for i in 0..n_bonds {
        let line_no = first_line + bonds_start + i;
        let t: Vec<&str> = lines[bonds_start + i].split_whitespace().collect();
        if t.len() < 3 {
            return Err(parse_err(
                FORMAT,
                line_no,
                "bond line needs two atoms and an order",
            ));
        }
        let a: usize = t[0]
            .parse()
            .map_err(|_| bad_number(FORMAT, line_no, t[0]))?;
        let b: usize = t[1]
            .parse()
            .map_err(|_| bad_number(FORMAT, line_no, t[1]))?;
        let order = match t[2] {
            "2" => BondOrder::Double,
            "3" => BondOrder::Triple,
            "4" => BondOrder::Aromatic,
            _ => BondOrder::Single,
        };
        if a == 0 || b == 0 || a > n_atoms || b > n_atoms {
            return Err(parse_err(
                FORMAT,
                line_no,
                format!("bond refers to atom {a}/{b} outside 1..={n_atoms}"),
            ));
        }
        mol.add_bond((a - 1).into(), (b - 1).into(), order)
            .map_err(|e| parse_err(FORMAT, line_no, e.to_string()))?;
    }

    // Property block: read `M  CHG` entries into the total formal charge.
    let mut charge = 0;
    for (i, line) in lines[bonds_start + n_bonds..].iter().enumerate() {
        let line_no = first_line + bonds_start + n_bonds + i;
        if line.trim_end() == "M  END" {
            break;
        }
        if let Some(rest) = line.strip_prefix("M  CHG") {
            let t: Vec<&str> = rest.split_whitespace().collect();
            let count: usize = t
                .first()
                .ok_or_else(|| parse_err(FORMAT, line_no, "M  CHG without a count"))?
                .parse()
                .map_err(|_| bad_number(FORMAT, line_no, t[0]))?;
            for pair in t[1..].chunks(2).take(count) {
                if let [_, value] = pair {
                    charge += value
                        .parse::<i32>()
                        .map_err(|_| bad_number(FORMAT, line_no, value))?;
                }
            }
        }
    }
    mol.set_formal_charge(charge);
    Ok(mol)
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER: &str = "water
  test

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.1173 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.7572   -0.4692 H   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.7572   -0.4692 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
M  END
$$$$
";

    #[test]
    fn parses_water() {
        let mol = parse_sdf(WATER).unwrap();
        assert_eq!(mol.name(), "water");
        assert_eq!(mol.formula(), "H2O");
        assert_eq!(mol.n_bonds(), 2);
        assert!((mol.atom(0usize.into()).pos.z - 0.1173).abs() < 1e-12);
    }

    #[test]
    fn round_trip_keeps_atoms_bonds_and_orders() {
        let mut mol = Molecule::new("ethylene");
        let c1 = mol.add_atom_raw(6, Vec3::new(0.0, 0.0, 0.0));
        let c2 = mol.add_atom_raw(6, Vec3::new(1.33, 0.0, 0.0));
        mol.add_bond(c1, c2, BondOrder::Double).unwrap();
        let text = write_sdf(&mol);
        let back = parse_sdf(&text).unwrap();
        assert_eq!(back.name(), "ethylene");
        assert_eq!(back.len(), 2);
        assert_eq!(back.bonds().next().unwrap().order, BondOrder::Double);
        assert!((back.atom(1usize.into()).pos.x - 1.33).abs() < 1e-4);
    }

    #[test]
    fn multi_record_and_charge() {
        let charged = WATER.replace("M  END", "M  CHG  1   1   1\nM  END");
        let text = format!("{charged}{WATER}");
        let mols = parse_sdf_all(&text).unwrap();
        assert_eq!(mols.len(), 2);
        assert_eq!(mols[0].formal_charge(), 1);
        assert_eq!(mols[1].formal_charge(), 0);
    }

    #[test]
    fn final_record_without_terminator_is_read() {
        let text = WATER.trim_end().trim_end_matches("$$$$");
        assert_eq!(parse_sdf_all(text).unwrap().len(), 1);
    }

    #[test]
    fn rejects_bad_records() {
        assert!(matches!(
            parse_sdf(
                "name\n\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0 0 0 Zz\nM  END\n$$$$\n"
            ),
            Err(IoError::UnknownElement { .. })
        ));
        assert!(
            parse_sdf("a\n\n\n  2  1  0\n  0 0 0 H\n  1 0 0 H\n  1  5  1  0\nM  END\n$$$$\n")
                .is_err()
        );
        assert!(matches!(parse_sdf(""), Err(IoError::UnexpectedEof { .. })));
    }
}
