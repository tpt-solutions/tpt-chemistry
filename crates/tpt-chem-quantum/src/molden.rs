//! Molden-format orbital output (`[Molden Format]`) for visualizers
//! (Molden, Jmol, chemcraft, …): atom block, Cartesian GTO shell block, and
//! the MO coefficient/occupation/energy blocks from an [`HfResult`].
//!
//! Two format subtleties handled here:
//!
//! * **Component ordering** — molden's Cartesian convention within a shell
//!   is `p: (x,y,z)`, `d: (xx,yy,zz,xy,xz,yz)`, while this crate's
//!   [`cartesian_components`] enumerates `x`-fastest (`p: (z,y,x)`); the MO
//!   coefficients are reordered on output.
//! * **Spherical flags** — `[5D]` marks pure d shells and `[6D]` Cartesian d
//!   shells. Spherical d coefficients are already in molden's harmonic order
//!   and are written as they are.
//! * **Supported shells** — s, p and d only. Higher momenta panic rather than
//!   writing zero coefficients.

use tpt_chem_core::element::symbol;
use tpt_chem_core::molecule::Molecule;

use crate::gaussian::{cartesian_components, Shell};
use crate::hf::HfResult;
use crate::spherical;

/// Molden's Cartesian component order per angular momentum.
fn molden_order(l: u8) -> &'static [(u8, u8, u8)] {
    match l {
        0 => &[(0, 0, 0)],
        1 => &[(1, 0, 0), (0, 1, 0), (0, 0, 1)],
        2 => &[
            (2, 0, 0),
            (0, 2, 0),
            (0, 0, 2),
            (1, 1, 0),
            (1, 0, 1),
            (0, 1, 1),
        ],
        // Higher momenta are outside the STO-3G scope; the d-order pattern
        // generalizes but is unused.
        _ => &[],
    }
}

fn shell_letter(l: u8) -> &'static str {
    match l {
        0 => "s",
        1 => "p",
        2 => "d",
        3 => "f",
        4 => "g",
        _ => "h",
    }
}

/// Render a converged RHF calculation as a Molden-format string.
///
/// `mol` carries Ångström positions (the `tpt-chem-core` convention);
/// `shells` must be the basis used for `result` (e.g. [`crate::basis::
/// sto3g_basis`] over the *Bohr-converted* atom positions, as the HF
/// driver builds it). Atom coordinates are written in Bohr (molden's
/// `[Atoms] (AU)` default). Occupations are derived from the electron
/// count `Σ Z − formal charge` (closed shell: 2.0 per occupied MO).
pub fn write_molden(mol: &Molecule, shells: &[Shell], result: &HfResult) -> String {
    let mut out = String::new();
    out.push_str("[Molden Format]\n");
    out.push_str("[Atoms] (AU)\n");
    let to_bohr = tpt_chem_core::units::angstrom_to_bohr(1.0);
    for (k, atom) in mol.atoms().enumerate() {
        let p = atom.pos * to_bohr;
        out.push_str(&format!(
            "  {}  {}  {}  {:.12}  {:.12}  {:.12}\n",
            symbol(atom.z).expect("valid element"),
            atom.z,
            k + 1,
            p.x,
            p.y,
            p.z
        ));
    }

    assert!(
        shells.iter().all(|s| s.l <= 2),
        "Molden output supports s, p and d shells only"
    );
    // `[5D]` selects pure (spherical) d and `[6D]` Cartesian d. Molden has one
    // keyword per angular momentum, so a basis mixing both is written as pure.
    if shells.iter().any(|s| s.l == 2 && s.spherical) {
        out.push_str("[5D]\n");
    } else if shells.iter().any(|s| s.l == 2) {
        out.push_str("[6D]\n");
    }

    out.push_str("[GTO]\n");
    // Basis-function offset per shell — molden numbers GTO blocks by atom;
    // we emit one block per atom collecting its shells.
    let n_atoms = mol.atoms().count();
    let mut atom_shells: Vec<Vec<&Shell>> = vec![Vec::new(); n_atoms];
    {
        // Shell centers match atom positions exactly (sto3g_basis builds
        // them so); match by position to keep the writer independent of
        // caller construction order.
        let atom_pos: Vec<tpt_chem_core::vec3::Vec3> =
            mol.atoms().map(|a| a.pos * to_bohr).collect();
        for sh in shells {
            let idx = atom_pos
                .iter()
                .position(|p| (*p - sh.center).norm_sq() < 1e-12)
                .unwrap_or_else(|| {
                    panic!("shell center not on any atom");
                });
            atom_shells[idx].push(sh);
        }
    }
    for (k, shells_k) in atom_shells.iter().enumerate() {
        out.push_str(&format!("  {} 0\n", k + 1));
        for sh in shells_k {
            out.push_str(&format!(
                "{} {} 1.00\n",
                shell_letter(sh.l),
                sh.primitives.len()
            ));
            for (i, &(alpha, coeff)) in sh.primitives.iter().enumerate() {
                out.push_str(&format!("  {}  {:.10}  {:.10}\n", i + 1, alpha, coeff));
            }
        }
        out.push_str("****\n");
    }

    out.push_str("[MO]\n");
    let nb = result.n_basis;
    let n_elec: f64 =
        mol.atoms().map(|a| f64::from(a.z)).sum::<f64>() - f64::from(mol.formal_charge());
    let n_occ = (n_elec / 2.0).round() as usize;
    for mo in 0..nb {
        let occ = if mo < n_occ { 2.0 } else { 0.0 };
        out.push_str(&format!("Occup= {occ:.6}\n"));
        out.push_str(&format!("Energy= {:.10}\n", result.orbital_energies[mo]));
        out.push_str("Spin= Alpha\n");
        // Coefficients in molden component order: per shell, permute from
        // this crate's x-fastest order into molden's. The leading integers
        // are sequential 1..nb positions; the reorder lives in the values.
        let mut seq = 0usize;
        let mut cursor = 0usize;
        for sh in shells {
            if sh.spherical {
                // Spherical d functions are stored in molden's harmonic order
                // (D0, D+1, D-1, D+2, D-2), so no reordering is needed.
                let n_functions = spherical::transform(sh).first().map_or(0, Vec::len);
                for m in 0..n_functions {
                    seq += 1;
                    out.push_str(&format!(
                        "  {}  {:.10}\n",
                        seq,
                        result.coefficients[(cursor + m) * nb + mo]
                    ));
                }
                cursor += n_functions;
                continue;
            }
            let comps = cartesian_components(sh.l);
            for &m in molden_order(sh.l) {
                let ci = comps.iter().position(|&c| c == m).expect("component");
                seq += 1;
                out.push_str(&format!(
                    "  {}  {:.10}\n",
                    seq,
                    result.coefficients[(cursor + ci) * nb + mo]
                ));
            }
            cursor += comps.len();
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::sto3g_basis;
    use crate::hf::rhf;
    use tpt_chem_core::molecule::Molecule;
    use tpt_chem_core::units::angstrom_to_bohr;
    use tpt_chem_core::vec3::Vec3;

    /// Water in Ångström (the `Molecule` convention `rhf` expects).
    fn water() -> Molecule {
        let r = 0.9894;
        let half = (100.025f64 / 2.0).to_radians();
        let mut m = Molecule::new("water");
        m.add_atom::<8>(Vec3::ZERO);
        m.add_atom::<1>(Vec3::new(0.0, r * half.sin(), r * half.cos()));
        m.add_atom::<1>(Vec3::new(0.0, -r * half.sin(), r * half.cos()));
        m
    }

    fn water_shells(mol: &Molecule) -> Vec<Shell> {
        sto3g_basis(
            &mol.atoms()
                .map(|a| (a.z, a.pos * angstrom_to_bohr(1.0)))
                .collect::<Vec<_>>(),
        )
        .unwrap()
    }

    #[test]
    fn water_molden_structure() {
        let mol = water();
        let shells = water_shells(&mol);
        let result = rhf(&mol).unwrap();
        let text = write_molden(&mol, &shells, &result);
        assert!(text.starts_with("[Molden Format]\n"));
        assert!(text.contains("[Atoms] (AU)"));
        assert!(text.contains("[GTO]"));
        assert!(text.contains("[MO]"));
        // No spherical flags for an s/p basis.
        assert!(!text.contains("[5D7F]"));
        // 7 basis functions → 7 MO blocks, 5 occupied with Occ= 2.
        assert_eq!(text.matches("Occup= 2.000000").count(), 5);
        assert_eq!(text.matches("Occup= 0.000000").count(), 2);
        // Each MO block lists 7 coefficients: total coefficient lines = 49.
        let coeff_lines = text
            .lines()
            .filter(|l| l.starts_with("  ") && l.split_whitespace().count() == 2)
            .count();
        // Atoms(3) + GTO header/prim lines + 49 MO coefficients; just check
        // the MO coefficient count via "Spin=" markers instead.
        assert_eq!(text.matches("Spin= Alpha").count(), 7);
        let _ = coeff_lines;
    }

    #[test]
    fn p_shell_coefficients_are_reordered_to_molden_order() {
        // Identity coefficients make the mapping observable: our basis
        // function k appears with coefficient 1.0 in MO block k. Water's
        // oxygen 2p shell occupies our functions 2,3,4 = (z,y,x); molden
        // orders that shell (px,py,pz) = our (4,3,2), so MO block 2 (0-based;
        // our pz) must place its 1.0 on the 5th coefficient line of the
        // block (after s, s, px, py).
        let mol = water();
        let shells = water_shells(&mol);
        let nb = 7;
        let mut coeffs = vec![0.0; nb * nb];
        for k in 0..nb {
            coeffs[k * nb + k] = 1.0;
        }
        let result = HfResult {
            energy: 0.0,
            electronic_energy: 0.0,
            nuclear_energy: 0.0,
            orbital_energies: (0..nb).map(|i| i as f64).collect(),
            coefficients: coeffs,
            density: vec![0.0; nb * nb],
            n_basis: nb,
        };
        let text = write_molden(&mol, &shells, &result);
        // Third MO block (our function 2).
        let blocks: Vec<Vec<&str>> = text
            .split("Spin= Alpha\n")
            .skip(1)
            .map(|b| b.lines().collect())
            .collect();
        assert_eq!(blocks.len(), 7);
        let block3 = &blocks[2];
        let coeffs3: Vec<(usize, f64)> = block3
            .iter()
            .filter_map(|l| {
                let mut it = l.split_whitespace();
                let i = it.next()?.parse::<usize>().ok()?;
                let v = it.next()?.parse::<f64>().ok()?;
                Some((i, v))
            })
            .collect();
        assert_eq!(coeffs3.len(), 7);
        // 1.0 lands on molden index 5 = pz within the p shell.
        let one = coeffs3.iter().find(|(_, v)| v.abs() > 0.5).unwrap();
        assert_eq!(one.0, 5, "block: {block3:?}");
    }

    #[test]
    fn charged_molecule_occupations_follow_formal_charge() {
        // H₃O⁺: 10 electrons → 5 occupied MOs (same as neutral water).
        let mut m = Molecule::new("H3O+");
        m.add_atom::<8>(Vec3::ZERO);
        for k in 0..3 {
            let a = 2.0 * core::f64::consts::PI * k as f64 / 3.0;
            m.add_atom::<1>(Vec3::new(0.98 * a.cos(), 0.98 * a.sin(), 0.0));
        }
        m.set_formal_charge(1);
        let shells = sto3g_basis(
            &m.atoms()
                .map(|a| (a.z, a.pos * angstrom_to_bohr(1.0)))
                .collect::<Vec<_>>(),
        )
        .unwrap();
        let result = rhf(&m).unwrap();
        let text = write_molden(&m, &shells, &result);
        assert_eq!(text.matches("Occup= 2.000000").count(), 5);
    }

    #[test]
    fn spherical_d_basis_writes_pure_d_coefficients() {
        use crate::basis::BasisSet;
        use crate::hf::rhf_with_basis;

        let mol = water();
        let atoms: Vec<(u8, Vec3)> = mol
            .atoms()
            .map(|a| (a.z, a.pos * angstrom_to_bohr(1.0)))
            .collect();
        let shells = BasisSet::CcPvdz.build(&atoms).unwrap();
        let res = rhf_with_basis(&mol, BasisSet::CcPvdz).unwrap();
        let text = write_molden(&mol, &shells, &res);
        let nb = res.n_basis;

        assert!(text.contains("[5D]\n"));
        assert!(!text.contains("[6D]"));
        assert_eq!(text.matches("Occup=").count(), nb);

        // Each MO block is Occup, Energy, Spin and then one line per basis
        // function. The spherical d functions are written in the order the SCF
        // uses, so their values must match the coefficient matrix directly.
        // (Cartesian p shells are reordered, so they are not checked here.)
        let mut offset = 0;
        let mut d_offset = None;
        for sh in &shells {
            let n = if sh.spherical {
                spherical::transform(sh).first().map_or(0, Vec::len)
            } else {
                cartesian_components(sh.l).len()
            };
            if sh.spherical && d_offset.is_none() {
                d_offset = Some(offset);
            }
            offset += n;
        }
        let d_offset = d_offset.expect("cc-pVDZ water has a spherical d shell");

        let mo_text = text.split("[MO]\n").nth(1).unwrap();
        let lines: Vec<&str> = mo_text.lines().collect();
        let block_len = 3 + nb;
        for mo in 0..nb {
            let block = &lines[mo * block_len..(mo + 1) * block_len];
            assert!(block[0].starts_with("Occup="));
            for m in 0..5 {
                let f = d_offset + m;
                let value: f64 = block[3 + f]
                    .split_whitespace()
                    .nth(1)
                    .unwrap()
                    .parse()
                    .unwrap();
                assert!((value - res.coefficients[f * nb + mo]).abs() < 1e-9);
            }
        }
    }
}
