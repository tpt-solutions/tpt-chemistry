//! Tripos MOL2 format: atoms with sybyl types and charges, ordered bonds.
//!
//! Sections handled: `@<TRIPOS>MOLECULE`, `@<TRIPOS>ATOM`,
//! `@<TRIPOS>BOND`. Other sections (`SUBSTRUCTURE`, …) are skipped on
//! read; only the three above are written.

use tpt_chem_core::element;
use tpt_chem_core::molecule::{BondOrder, Molecule};
use tpt_chem_core::vec3::Vec3;

use crate::error::{bad_number, parse_err, parse_f64, IoError, Result};

const FORMAT: &str = "mol2";

/// Element atomic number from a sybyl atom type such as `C.3`, `N.ar`, or
/// `Cl`. Returns `None` for pseudo-types like `Du` (dummy) or `LP`
/// (lone pair).
pub fn element_from_sybyl(atom_type: &str) -> Option<u8> {
    let letters: String = atom_type
        .chars()
        .take_while(|c| c.is_ascii_alphabetic())
        .collect();
    match letters.as_str() {
        "Du" | "LP" | "Any" | "Hal" | "Het" | "Hev" | "H" => {
            // `H` alone is genuine hydrogen; handled by the fallthrough
            // below, only pseudo-types are filtered here.
            if letters == "H" {
                Some(1)
            } else {
                None
            }
        }
        s => element::from_symbol(s),
    }
}

fn bond_order_from_sybyl(t: &str) -> Option<BondOrder> {
    match t.trim() {
        "1" => Some(BondOrder::Single),
        "2" => Some(BondOrder::Double),
        "3" => Some(BondOrder::Triple),
        "ar" | "am" => Some(BondOrder::Aromatic),
        _ => None,
    }
}

/// Parse a MOL2 file.
///
/// Charges are stored on the atoms; the molecule name comes from the
/// `MOLECULE` block.
///
/// # Errors
/// [`IoError`] on missing sections, malformed rows, or unknown elements.
pub fn parse_mol2(text: &str) -> Result<Molecule> {
    let mut name = "mol2";
    let mut saw_atoms = false;
    let mut mol: Option<Molecule> = None;
    let mut section = "";

    for (i, raw) in text.lines().enumerate() {
        let line_no = i + 1;
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(marker) = line.strip_prefix("@<TRIPOS>") {
            section = marker.trim();
            if section == "MOLECULE" {
                // The name is the next non-empty line; skip the rest of the
                // block until the following marker.
                section = "MOLECULE_AWAIT_NAME";
            }
            continue;
        }
        match section {
            "MOLECULE_AWAIT_NAME" => {
                name = line;
                section = "MOLECULE_SKIP";
            }
            "MOLECULE_SKIP" => { /* counts/type lines — ignored */ }
            "ATOM" => {
                let tokens: Vec<&str> = line.split_whitespace().collect();
                if tokens.len() < 6 {
                    return Err(parse_err(
                        FORMAT,
                        line_no,
                        format!("ATOM row needs ≥6 fields, got {}", tokens.len()),
                    ));
                }
                let mol = mol.get_or_insert_with(|| Molecule::new(name));
                let x = parse_f64(FORMAT, line_no, tokens[2])?;
                let y = parse_f64(FORMAT, line_no, tokens[3])?;
                let z = parse_f64(FORMAT, line_no, tokens[4])?;
                let z_elem =
                    element_from_sybyl(tokens[5]).ok_or_else(|| IoError::UnknownElement {
                        format: FORMAT,
                        symbol: tokens[5].to_string(),
                        line: Some(line_no),
                    })?;
                let id = mol.add_atom_labeled(z_elem, Vec3::new(x, y, z), tokens[1]);
                if tokens.len() > 8 {
                    if let Ok(q) = parse_f64(FORMAT, line_no, tokens[8]) {
                        mol.atom_mut(id).charge = q;
                    }
                }
                saw_atoms = true;
            }
            "BOND" => {
                let tokens: Vec<&str> = line.split_whitespace().collect();
                if tokens.len() < 4 {
                    return Err(parse_err(
                        FORMAT,
                        line_no,
                        format!("BOND row needs ≥4 fields, got {}", tokens.len()),
                    ));
                }
                let mol = mol.as_mut().ok_or(IoError::MissingSection {
                    format: FORMAT,
                    section: "@<TRIPOS>ATOM".into(),
                })?;
                let a: usize = tokens[1]
                    .parse()
                    .map_err(|_| bad_number(FORMAT, line_no, tokens[1]))?;
                let b: usize = tokens[2]
                    .parse()
                    .map_err(|_| bad_number(FORMAT, line_no, tokens[2]))?;
                let order = bond_order_from_sybyl(tokens[3]).unwrap_or(BondOrder::Single);
                if a == 0 || b == 0 || a > mol.len() || b > mol.len() {
                    return Err(parse_err(
                        FORMAT,
                        line_no,
                        format!("bond refers to atom {a}/{b} outside 1..={}", mol.len()),
                    ));
                }
                mol.add_bond((a - 1).into(), (b - 1).into(), order)
                    .map_err(|e| parse_err(FORMAT, line_no, e.to_string()))?;
            }
            _ => {} // unknown sections skipped
        }
    }

    match (mol, saw_atoms) {
        (Some(m), true) => Ok(m),
        _ => Err(IoError::MissingSection {
            format: FORMAT,
            section: "@<TRIPOS>ATOM".into(),
        }),
    }
}

/// Serialize a molecule as a MOL2 string.
///
/// Sybyl types are emitted as bare element symbols; bond orders map back
/// to `1/2/3/ar`.
pub fn write_mol2(mol: &Molecule) -> String {
    let mut out = String::new();
    out.push_str("@<TRIPOS>MOLECULE\n");
    out.push_str(mol.name());
    out.push('\n');
    out.push_str(&format!("{} {}\n", mol.len(), mol.n_bonds()));
    out.push_str("SMALL\nUSER_CHARGES\n\n");
    out.push_str("@<TRIPOS>ATOM\n");
    for (i, atom) in mol.atoms().enumerate() {
        let sym = element::symbol(atom.z).unwrap_or("X");
        let name = atom.label.as_deref().unwrap_or(sym);
        out.push_str(&format!(
            "{:>7} {:<8} {:>10.4} {:>10.4} {:>10.4} {:<6} 1 MOL {:<10.4}\n",
            i + 1,
            name,
            atom.pos.x,
            atom.pos.y,
            atom.pos.z,
            sym,
            atom.charge
        ));
    }
    out.push_str("@<TRIPOS>BOND\n");
    for (i, bond) in mol.bonds().enumerate() {
        let t = match bond.order {
            BondOrder::Single => "1",
            BondOrder::Double => "2",
            BondOrder::Triple => "3",
            BondOrder::Aromatic => "ar",
        };
        out.push_str(&format!(
            "{:>6} {:>5} {:>5} {}\n",
            i + 1,
            bond.a.index() + 1,
            bond.b.index() + 1,
            t
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const METHANOL: &str = "\
@<TRIPOS>MOLECULE
methanol
 6 5 1 0 0
SMALL
USER_CHARGES

@<TRIPOS>ATOM
      1 C1          1.0000    -0.1000     0.0000 C.3       1 MOL      -0.1800
      2 O1          2.4500     0.0500     0.0000 O.3       1 MOL      -0.4500
      3 H1          0.6000     0.4000     0.9000 H         1 MOL       0.0400
      4 H2          0.6000    -1.1000     0.0000 H         1 MOL       0.0400
      5 H3          0.6000     0.4000    -0.9000 H         1 MOL       0.0400
      6 H4          2.7500    -0.8000     0.5000 H         1 MOL       0.4300
@<TRIPOS>BOND
     1     1     2 1
     2     1     3 1
     3     1     4 1
     4     1     5 1
     5     2     6 1
";

    #[test]
    fn parse_methanol() {
        let mol = parse_mol2(METHANOL).unwrap();
        assert_eq!(mol.name(), "methanol");
        assert_eq!(mol.formula(), "CH4O");
        assert_eq!(mol.n_bonds(), 5);
        assert!((mol.atom(0usize.into()).charge - (-0.18)).abs() < 1e-12);
        assert!((mol.atom(1usize.into()).pos.x - 2.45).abs() < 1e-12);
    }

    #[test]
    fn roundtrip() {
        let mol = parse_mol2(METHANOL).unwrap();
        let text = write_mol2(&mol);
        let again = parse_mol2(&text).unwrap();
        assert_eq!(again.len(), 6);
        assert_eq!(again.n_bonds(), 5);
        assert_eq!(again.name(), "methanol");
        for i in 0..6usize {
            assert!((again.atom(i.into()).pos.x - mol.atom(i.into()).pos.x).abs() < 1e-3);
            assert!((again.atom(i.into()).charge - mol.atom(i.into()).charge).abs() < 1e-3);
        }
    }

    #[test]
    fn bond_orders_map() {
        let text = "@<TRIPOS>MOLECULE\nethene\n2 1\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n1 C 0 0 0 C.2 1 R 0.0\n2 C 1.34 0 0 C.2 1 R 0.0\n@<TRIPOS>BOND\n1 1 2 2\n";
        let mol = parse_mol2(text).unwrap();
        assert_eq!(mol.bonds().next().unwrap().order, BondOrder::Double);
    }

    #[test]
    fn errors() {
        assert!(matches!(
            parse_mol2("no sections here"),
            Err(IoError::MissingSection { .. })
        ));
        let bad = "@<TRIPOS>MOLECULE\nx\n1 0\nSMALL\nNO_CHARGES\n\n@<TRIPOS>ATOM\n1 C 0 0 zz C.3 1 R 0.0\n";
        assert!(matches!(parse_mol2(bad), Err(IoError::BadNumber { .. })));
    }

    #[test]
    fn sybyl_elements() {
        assert_eq!(element_from_sybyl("C.3"), Some(6));
        assert_eq!(element_from_sybyl("N.ar"), Some(7));
        assert_eq!(element_from_sybyl("Cl"), Some(17));
        assert_eq!(element_from_sybyl("H"), Some(1));
        assert_eq!(element_from_sybyl("Du"), None);
    }
}
