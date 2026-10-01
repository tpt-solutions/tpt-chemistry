//! XYZ format: the minimal, ubiquitous geometry format.
//!
//! Layout: an atom-count line, a free-text comment line, then one line per
//! atom (`symbol x y z [extras…]`, whitespace-separated, Å). Multiple
//! frames may be concatenated — that is the de-facto MD trajectory format,
//! and [`super::trajectory::XyzTrajectoryReader`] streams it.

use std::io::{BufRead, Write};

use tpt_chem_core::element;
use tpt_chem_core::molecule::Molecule;
use tpt_chem_core::vec3::Vec3;

use crate::error::{io_err, parse_err, parse_f64, IoError, Result};

const FORMAT: &str = "xyz";

/// Parse the first XYZ frame from `text`.
///
/// The comment line becomes the molecule name.
///
/// # Errors
/// [`IoError`] on malformed input or unknown element symbols.
pub fn parse_xyz(text: &str) -> Result<Molecule> {
    let mut reader = text.as_bytes();
    read_frame(&mut reader)?.ok_or(IoError::UnexpectedEof { format: FORMAT })
}

/// Parse every XYZ frame in `text` (a trajectory).
///
/// # Errors
/// [`IoError`] on malformed input or unknown element symbols.
pub fn parse_xyz_all(text: &str) -> Result<Vec<Molecule>> {
    let mut reader = text.as_bytes();
    let mut frames = Vec::new();
    while let Some(mol) = read_frame(&mut reader)? {
        frames.push(mol);
    }
    Ok(frames)
}

/// Serialize one frame as an XYZ string.
pub fn write_xyz(mol: &Molecule, comment: &str) -> String {
    let mut out = String::new();
    out.push_str(&format!("{}\n{}\n", mol.len(), comment));
    for atom in mol.atoms() {
        let sym = element::symbol(atom.z).unwrap_or("X");
        out.push_str(&format!(
            "{sym:<4} {:>12.6} {:>12.6} {:>12.6}\n",
            atom.pos.x, atom.pos.y, atom.pos.z
        ));
    }
    out
}

/// Write one frame to a writer (used by the trajectory writer).
pub(crate) fn write_frame_to<W: Write>(
    w: &mut W,
    mol: &Molecule,
    comment: &str,
) -> std::io::Result<()> {
    write!(w, "{}\n{}\n", mol.len(), comment)?;
    for atom in mol.atoms() {
        let sym = element::symbol(atom.z).unwrap_or("X");
        writeln!(
            w,
            "{sym:<4} {:>12.6} {:>12.6} {:>12.6}",
            atom.pos.x, atom.pos.y, atom.pos.z
        )?;
    }
    Ok(())
}

/// Read one XYZ frame; `Ok(None)` at clean end of stream.
///
/// Blank lines between frames are tolerated (they are not part of the
/// format but appear in real trajectory files).
pub(crate) fn read_frame<R: BufRead>(reader: &mut R) -> Result<Option<Molecule>> {
    let mut count_line = String::new();
    loop {
        count_line.clear();
        let n = reader.read_line(&mut count_line).map_err(io_err)?;
        if n == 0 {
            return Ok(None); // clean EOF
        }
        if count_line.trim().is_empty() {
            continue; // skip inter-frame blank line
        }
        break;
    }
    let natoms: usize = count_line.trim().parse().map_err(|_| IoError::BadNumber {
        format: FORMAT,
        line: None,
        token: count_line.trim().to_string(),
    })?;

    let mut comment = String::new();
    reader.read_line(&mut comment).map_err(io_err)?;

    let mut mol = Molecule::new(comment.trim());
    for i in 0..natoms {
        let mut line = String::new();
        let n = reader.read_line(&mut line).map_err(io_err)?;
        if n == 0 {
            return Err(IoError::Parse {
                format: FORMAT,
                line: Some(i + 3),
                message: "unexpected EOF inside frame".into(),
            });
        }
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.len() < 4 {
            return Err(parse_err(
                FORMAT,
                i + 3,
                format!("expected ≥4 fields, got {}", tokens.len()),
            ));
        }
        let z = element::from_symbol(tokens[0]).ok_or_else(|| IoError::UnknownElement {
            format: FORMAT,
            symbol: tokens[0].to_string(),
            line: Some(i + 3),
        })?;
        let x = parse_f64(FORMAT, i + 3, tokens[1])?;
        let y = parse_f64(FORMAT, i + 3, tokens[2])?;
        let zc = parse_f64(FORMAT, i + 3, tokens[3])?;
        mol.add_atom_raw(z, Vec3::new(x, y, zc));
    }
    Ok(Some(mol))
}

#[cfg(test)]
mod tests {
    use super::*;

    const WATER: &str =
        "3\nwater\nO  0.000  0.000  0.117\nH  0.000  0.757 -0.469\nH  0.000 -0.757 -0.469\n";

    #[test]
    fn parse_write_roundtrip() {
        let mol = parse_xyz(WATER).unwrap();
        assert_eq!(mol.name(), "water");
        assert_eq!(mol.formula(), "H2O");
        assert!((mol.atom(0usize.into()).pos.z - 0.117).abs() < 1e-12);

        let text = write_xyz(&mol, "rt");
        let again = parse_xyz(&text).unwrap();
        assert_eq!(again.len(), 3);
        for i in 0..3usize {
            assert_eq!(again.atom(i.into()).pos, mol.atom(i.into()).pos);
            assert_eq!(again.atom(i.into()).z, mol.atom(i.into()).z);
        }
    }

    #[test]
    fn multi_frame_trajectory() {
        let text = format!("{WATER}{WATER}");
        let frames = parse_xyz_all(&text).unwrap();
        assert_eq!(frames.len(), 2);
        assert_eq!(frames[1].formula(), "H2O");
    }

    #[test]
    fn extras_and_blank_lines_tolerated() {
        let text = "1\nwith extras\nH 0.0 0.0 0.0 extra1 extra2\n\n1\nsecond\nHe 1 2 3\n";
        let frames = parse_xyz_all(text).unwrap();
        assert_eq!(frames.len(), 2);
        assert_eq!(frames[0].atom(0usize.into()).z, 1);
        assert!((frames[1].atom(0usize.into()).pos.z - 3.0).abs() < 1e-12);
    }

    #[test]
    fn errors() {
        // Not enough atoms provided.
        assert!(parse_xyz("2\nshort\nH 0 0 0\n").is_err());
        // Unknown element.
        assert!(matches!(
            parse_xyz("1\nx\nZz 0 0 0\n"),
            Err(IoError::UnknownElement { .. })
        ));
        // Bad count.
        assert!(matches!(
            parse_xyz("two\nx\nH 0 0 0\n"),
            Err(IoError::BadNumber { .. })
        ));
        // Empty input.
        assert!(matches!(parse_xyz(""), Err(IoError::UnexpectedEof { .. })));
    }
}
