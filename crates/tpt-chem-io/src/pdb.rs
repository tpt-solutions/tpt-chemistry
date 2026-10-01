//! PDB format: the protein Data Bank coordinate file.
//!
//! Hand-rolled fixed-column parsing of `ATOM`/`HETATM` records plus
//! `CONECT` bonds. Residue/chain metadata is intentionally not modeled on
//! [`Molecule`] (it is a flat molecular graph); atom names are preserved in
//! each atom's label, and the element is taken from columns 77–78 when
//! present with a PDB atom-name heuristic as fallback.

use tpt_chem_core::element;
use tpt_chem_core::molecule::{AtomId, BondOrder, Molecule};
use tpt_chem_core::vec3::Vec3;

use crate::error::{parse_err, parse_f64, IoError, Result};

const FORMAT: &str = "pdb";

/// 1-based inclusive column slice, PDB style. Returns `""` when the line is
/// shorter than the requested columns.
fn cols(line: &str, from: usize, to: usize) -> &str {
    if line.len() >= from {
        let end = to.min(line.len());
        &line[from - 1..end]
    } else {
        ""
    }
}

/// Parse a PDB file from text into a [`Molecule`].
///
/// `ATOM`/`HETATM` records supply atoms; `CONECT` records supply bonds
/// (written as [`BondOrder::Single`] — PDB does not carry bond orders).
///
/// # Errors
/// [`IoError`] on malformed coordinate fields or unknown elements.
pub fn parse_pdb(text: &str) -> Result<Molecule> {
    let mut mol = Molecule::new("pdb");
    // serial number -> atom index, for CONECT resolution.
    let mut serials: Vec<(i64, AtomId)> = Vec::new();

    for (i, raw) in text.lines().enumerate() {
        let line_no = i + 1;
        let record = cols(raw, 1, 6);
        match record.trim_end() {
            "ATOM" | "HETATM" => {
                let name = cols(raw, 13, 16).trim().to_string();
                let x = parse_f64(FORMAT, line_no, cols(raw, 31, 38))?;
                let y = parse_f64(FORMAT, line_no, cols(raw, 39, 46))?;
                let z = parse_f64(FORMAT, line_no, cols(raw, 47, 54))?;
                let elem_field = cols(raw, 77, 78).trim().to_string();
                let zelem = if elem_field.is_empty() {
                    element_from_name(&name).ok_or_else(|| IoError::UnknownElement {
                        format: FORMAT,
                        symbol: name.clone(),
                        line: Some(line_no),
                    })?
                } else {
                    element::from_symbol(&elem_field).ok_or(IoError::UnknownElement {
                        format: FORMAT,
                        symbol: elem_field,
                        line: Some(line_no),
                    })?
                };
                let id = mol.add_atom_labeled(zelem, Vec3::new(x, y, z), &name);
                if let Ok(serial) = cols(raw, 7, 11).trim().parse::<i64>() {
                    serials.push((serial, id));
                }
            }
            "CONECT" => {
                let nums: Vec<&str> = raw[6..].split_whitespace().collect();
                if nums.is_empty() {
                    return Err(parse_err(FORMAT, line_no, "CONECT without atoms"));
                }
                let from = nums[0]
                    .parse::<i64>()
                    .map_err(|_| crate::error::bad_number(FORMAT, line_no, nums[0]))?;
                let from_id = find_serial(&serials, from);
                for tok in &nums[1..] {
                    let to = tok
                        .parse::<i64>()
                        .map_err(|_| crate::error::bad_number(FORMAT, line_no, tok))?;
                    if let (Some(a), Some(b)) = (from_id, find_serial(&serials, to)) {
                        // Duplicate CONECT lines are common; ignore dups.
                        let _ = mol.add_bond(a, b, BondOrder::Single);
                    }
                }
            }
            _ => {} // HEADER, TITLE, TER, END, … — ignored
        }
    }

    if mol.is_empty() {
        return Err(IoError::Parse {
            format: FORMAT,
            line: None,
            message: "no ATOM/HETATM records found".into(),
        });
    }
    Ok(mol)
}

fn find_serial(serials: &[(i64, AtomId)], serial: i64) -> Option<AtomId> {
    serials
        .iter()
        .find(|(s, _)| *s == serial)
        .map(|(_, id)| *id)
}

/// PDB atom-name element heuristic for files missing columns 77–78.
///
/// Works on the raw 4-column name field: the element starts at the first
/// alphanumeric character; a digit there (e.g. `1H5'`) is skipped; a
/// lowercase second letter signals a two-letter element (`Cl`, `Fe`).
/// Right-justified uppercase names like `" FE "` remain ambiguous — such
/// files are expected to carry the element in columns 77–78.
pub fn element_from_name(name: &str) -> Option<u8> {
    let bytes = name.as_bytes();
    let idx = bytes.iter().position(|b| b.is_ascii_alphanumeric())?;
    let first = bytes[idx] as char;
    if first.is_ascii_digit() {
        let second = bytes.get(idx + 1).map(|b| *b as char)?;
        return element::from_symbol(&second.to_string());
    }
    match bytes.get(idx + 1).map(|b| *b as char) {
        Some(c) if c.is_ascii_lowercase() => element::from_symbol(&name[idx..idx + 2]),
        _ => element::from_symbol(&first.to_string()),
    }
}

/// Serialize a molecule as a minimal PDB file (ATOM + CONECT records).
///
/// Chain `A`, residue `MOL`, residue number 1, occupancy 1.0.
pub fn write_pdb(mol: &Molecule) -> String {
    let mut out = String::new();
    out.push_str(&format!(
        "REMARK   {} atoms, {} bonds\n",
        mol.len(),
        mol.n_bonds()
    ));
    for (i, atom) in mol.atoms().enumerate() {
        let sym = element::symbol(atom.z).unwrap_or("X");
        let name = atom.label.as_deref().unwrap_or(sym);
        // PDB convention: names of ≤3 characters are offset one column.
        let name_col = if name.len() < 4 {
            format!(" {name:<3}")
        } else {
            format!("{name:<4}")
        };
        out.push_str(&format!(
            "ATOM  {:>5} {} MOL A   1    {:>8.3}{:>8.3}{:>8.3}  1.00  0.00          {:>2}\n",
            i + 1,
            name_col,
            atom.pos.x,
            atom.pos.y,
            atom.pos.z,
            sym
        ));
    }
    for bond in mol.bonds() {
        out.push_str(&format!(
            "CONECT{:>5}{:>5}\n",
            bond.a.index() + 1,
            bond.b.index() + 1
        ));
    }
    out.push_str("END\n");
    out
}

/// Serialize and include explicit bond orders as REMARKs (PDB cannot store
/// them natively); lossless for [`BondOrder`] only through this side channel.
pub fn write_pdb_with_orders(mol: &Molecule) -> Result<String> {
    let mut out = write_pdb(mol);
    for (i, bond) in mol.bonds().enumerate() {
        let code = match bond.order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 4,
        };
        out.push_str(&format!(
            "REMARK BOND {} {} {} {}\n",
            i + 1,
            bond.a.index() + 1,
            bond.b.index() + 1,
            code
        ));
    }
    Ok(out)
}

/// Convenience: parse ignoring [`MolError`]s from duplicate CONECT lines.
pub fn parse_pdb_lenient(text: &str) -> Result<Molecule> {
    parse_pdb(text)
}

#[cfg(test)]
mod tests {
    use super::*;

    const GLYCINE: &str = "\
HEADER    TEST
ATOM      1  N   GLY A   1      -0.532   1.369   0.000  1.00  0.00           N
ATOM      2  CA  GLY A   1       0.414   0.275   0.000  1.00  0.00           C
ATOM      3  C   GLY A   1       1.766   0.698   0.000  1.00  0.00           C
ATOM      4  O   GLY A   1       2.745  -0.047   0.000  1.00  0.00           O
ATOM      5  Cl  GLY A   1       3.000   1.000   1.000  1.00  0.00          Cl
CONECT    1     2
CONECT    2     3
CONECT    3     4
END
";

    #[test]
    fn parse_glycine() {
        let mol = parse_pdb(GLYCINE).unwrap();
        assert_eq!(mol.len(), 5);
        assert_eq!(mol.formula(), "C2ClNO"); // no hydrogens in this stub
        assert_eq!(mol.atom(0usize.into()).label.as_deref(), Some("N"));
        assert!((mol.atom(1usize.into()).pos.x - 0.414).abs() < 1e-12);
        assert_eq!(mol.n_bonds(), 3);
        assert_eq!(mol.atom(4usize.into()).z, 17); // Cl from columns 77-78
    }

    #[test]
    fn roundtrip() {
        let mol = parse_pdb(GLYCINE).unwrap();
        let text = write_pdb(&mol);
        let again = parse_pdb(&text).unwrap();
        assert_eq!(again.len(), mol.len());
        assert_eq!(again.n_bonds(), mol.n_bonds());
        for i in 0..mol.len() {
            let (a, b) = (mol.atom(i.into()), again.atom(i.into()));
            assert_eq!(a.z, b.z);
            for c in 0..3 {
                assert!(
                    (a.pos.get(c) - b.pos.get(c)).abs() < 5e-4,
                    "atom {i} coord {c}"
                );
            }
        }
    }

    #[test]
    fn element_heuristic() {
        assert_eq!(element_from_name(" CA "), Some(6)); // alpha carbon
        assert_eq!(element_from_name("Cl1 "), Some(17));
        assert_eq!(element_from_name("N   "), Some(7));
        assert_eq!(element_from_name("1H  "), Some(1));
        assert_eq!(element_from_name("Fe  "), Some(26)); // lowercase second letter
    }

    #[test]
    fn no_atoms_is_error() {
        assert!(parse_pdb("HEADER nothing\nEND\n").is_err());
    }
}
