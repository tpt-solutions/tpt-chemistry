//! GROMACS `.gro` structure format: fixed-column reader and writer.
//!
//! Line layout (all coordinates in **nanometers** on disk; this API takes
//! and returns Ångström and converts at the boundary):
//!
//! ```text
//! <title>
//! <n_atoms>
//! %5d%-5s%5s%5d%8.3f%8.3f%8.3f   (resid, resname, atomname, atomnum, x, y, z)
//! <box: 3 or 9 floats>            (v1x v2y v3z [v1y v1z v2x v2z v3x v3y])
//! ```
//!
//! Atom elements are recovered from the atom-name heuristic shared with the
//! PDB reader (leading alphabetic characters matched against the periodic
//! table).

use crate::error::{IoError, Result};
use tpt_chem_core::element;
use tpt_chem_core::molecule::Molecule;
use tpt_chem_core::vec3::Vec3;

/// nm → Å.
const NM_TO_ANG: f64 = 10.0;

/// Parse a `.gro` structure. The box line is returned alongside the
/// molecule as the orthorhombic diagonal `(v1x, v2y, v3z)` in Å, when
/// present.
pub fn parse_gro(text: &str) -> Result<(Molecule, Option<Vec3>)> {
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let title = lines
        .next()
        .ok_or(IoError::UnexpectedEof { format: "gro" })?
        .to_string();
    let n_atoms: usize = lines
        .next()
        .ok_or(IoError::UnexpectedEof { format: "gro" })?
        .trim()
        .parse()
        .map_err(|_| IoError::BadNumber {
            format: "gro",
            line: Some(2),
            token: "atom count".into(),
        })?;

    let mut mol = Molecule::new(&title);
    for (k, line) in lines.by_ref().take(n_atoms).enumerate() {
        let lineno = k + 3;
        // Fixed columns: resid(5) name(5) atomname(5) num(5) then 3×8.3f.
        if line.len() < 44 {
            return Err(IoError::Parse {
                format: "gro",
                line: Some(lineno),
                message: format!("atom line is {} chars, need 44", line.len()),
            });
        }
        let atom_name = line[10..15].trim();
        let z = element::from_symbol(&strip_digits(atom_name)).ok_or_else(|| {
            IoError::UnknownElement {
                format: "gro",
                symbol: atom_name.to_string(),
                line: Some(lineno),
            }
        })?;
        let parse = |s: &str| -> Result<f64> {
            s.trim().parse().map_err(|_| IoError::BadNumber {
                format: "gro",
                line: Some(lineno),
                token: s.trim().to_string(),
            })
        };
        let x = parse(&line[20..28])? * NM_TO_ANG;
        let y = parse(&line[28..36])? * NM_TO_ANG;
        let zc = parse(&line[36..44])? * NM_TO_ANG;
        mol.add_atom_raw(z, Vec3::new(x, y, zc));
    }
    if mol.len() != n_atoms {
        return Err(IoError::Parse {
            format: "gro",
            line: None,
            message: format!("declared {n_atoms} atoms, found {}", mol.len()),
        });
    }

    // Optional box line (one line of 3 or 9 floats).
    let box_line = lines.next();
    let box_ = box_line.and_then(|l| {
        let v: Vec<f64> = l
            .split_whitespace()
            .filter_map(|t| t.parse().ok())
            .collect();
        if v.len() >= 3 {
            Some(Vec3::new(
                v[0] * NM_TO_ANG,
                v[1] * NM_TO_ANG,
                v[2] * NM_TO_ANG,
            ))
        } else {
            None
        }
    });
    Ok((mol, box_))
}

/// Write a `.gro` structure: residue 1 named after the molecule (truncated
/// to 5 characters), atom names from element symbols, positions in Å
/// converted to nm. `box_diag` is the optional orthorhombic box (Å).
pub fn write_gro(mol: &Molecule, box_diag: Option<Vec3>) -> String {
    let mut out = String::new();
    out.push_str(mol.name());
    out.push('\n');
    out.push_str(&format!("{:5}\n", mol.len()));
    for (k, atom) in mol.atoms().enumerate() {
        let name = atom
            .label
            .clone()
            .unwrap_or_else(|| element::symbol(atom.z).unwrap_or("X").to_string());
        let p = atom.pos / NM_TO_ANG;
        out.push_str(&format!(
            "{:5}{:<5}{:>5}{:5}{:8.3}{:8.3}{:8.3}\n",
            1,
            truncate(mol.name(), 5),
            truncate(&name, 5),
            k + 1,
            p.x,
            p.y,
            p.z
        ));
    }
    if let Some(b) = box_diag {
        out.push_str(&format!(
            "   {:10.5}{:10.5}{:10.5}\n",
            b.x / NM_TO_ANG,
            b.y / NM_TO_ANG,
            b.z / NM_TO_ANG
        ));
    }
    out
}

fn truncate(s: &str, n: usize) -> String {
    s.chars().take(n).collect()
}

/// Strip trailing digits from an atom name: `"CA12"` → `"CA"`, `"H1"` → `"H"`.
fn strip_digits(name: &str) -> String {
    name.trim_end_matches(|c: char| c.is_ascii_digit())
        .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn water() -> Molecule {
        let mut m = Molecule::new("water");
        m.add_atom_raw(8, Vec3::new(0.0, 0.0, 0.117));
        m.add_atom_raw(1, Vec3::new(0.0, 0.757, -0.469));
        m.add_atom_raw(1, Vec3::new(0.0, -0.757, -0.469));
        m
    }

    #[test]
    fn round_trip() {
        let mol = water();
        let text = write_gro(&mol, Some(Vec3::new(30.0, 30.0, 30.0)));
        let (parsed, box_) = parse_gro(&text).unwrap();
        assert_eq!(parsed.formula(), "H2O");
        assert_eq!(parsed.name(), "water");
        for k in 0..3usize {
            let a = parsed.atom((k as u16).into()).pos;
            let b = mol.atom((k as u16).into()).pos;
            assert!((a - b).norm() < 5e-3, "atom {k}: {a:?} vs {b:?}");
        }
        let b = box_.unwrap();
        assert!((b.x - 30.0).abs() < 1e-6);
    }

    #[test]
    fn parses_reference_layout() {
        // A hand-written GROMACS-style file (nm coordinates, fixed columns).
        let text = "\
Argon dimer
    2
    1Ar      Ar    1   0.000   0.000   0.000
    1Ar      Ar    2   0.380   0.000   0.000
   2.00000   2.00000   2.00000
";
        let (mol, box_) = parse_gro(text).unwrap();
        assert_eq!(mol.len(), 2);
        assert_eq!(mol.atom(0usize.into()).z, 18);
        // 0.380 nm = 3.8 Å.
        assert!((mol.atom(1usize.into()).pos.x - 3.8).abs() < 1e-9);
        assert!((box_.unwrap().x - 20.0).abs() < 1e-9);
    }

    #[test]
    fn rejects_short_atom_line() {
        let text = "bad\n2\n  1Ar  Ar  1  0.0\n";
        assert!(parse_gro(text).is_err());
    }
}
