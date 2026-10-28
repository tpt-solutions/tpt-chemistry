//! The sample structures in `data/` must parse and have the expected contents.

use tpt_chem_io::cif::parse_cif;
use tpt_chem_io::parse_xyz;
use tpt_chem_io::pdb::parse_pdb;

#[test]
fn water_sample() {
    let mol = parse_xyz(include_str!("../../../data/water.xyz")).unwrap();
    assert_eq!(mol.formula(), "H2O");
    assert_eq!(mol.len(), 3);
}

#[test]
fn benzene_sample() {
    let mol = parse_xyz(include_str!("../../../data/benzene.xyz")).unwrap();
    assert_eq!(mol.formula(), "C6H6");
    // Ring C–C bonds are 1.39 Å.
    let d = mol.distance(0usize.into(), 1usize.into());
    assert!((d - 1.39).abs() < 1e-3, "C–C = {d}");
}

#[test]
fn nacl_sample() {
    let mol = parse_xyz(include_str!("../../../data/nacl_diatomic.xyz")).unwrap();
    assert_eq!(mol.len(), 2);
    let d = mol.distance(0usize.into(), 1usize.into());
    assert!((d - 2.36).abs() < 1e-9);
}

/// Smallest pairwise distance in a molecule (Å).
fn min_pair_distance(mol: &tpt_chem_core::prelude::Molecule) -> f64 {
    let n = mol.len();
    let mut best = f64::INFINITY;
    for i in 0..n {
        for j in (i + 1)..n {
            best = best.min(mol.distance(i.into(), j.into()));
        }
    }
    best
}

#[test]
fn argon_fcc_sample() {
    let mol = parse_xyz(include_str!("../../../data/argon_fcc.xyz")).unwrap();
    assert_eq!(mol.len(), 256);
    assert_eq!(mol.formula(), "Ar256");
    // FCC nearest neighbours sit at a/√2 with a = 5.26 Å.
    let expected = 5.26 / 2f64.sqrt();
    assert!((min_pair_distance(&mol) - expected).abs() < 1e-6);
}

#[test]
fn silicon_diamond_sample() {
    let mol = parse_xyz(include_str!("../../../data/silicon_diamond.xyz")).unwrap();
    assert_eq!(mol.len(), 64);
    assert_eq!(mol.formula(), "Si64");
    // Diamond-cubic nearest neighbours sit at a√3/4 with a = 5.431 Å.
    let expected = 5.431 * 3f64.sqrt() / 4.0;
    assert!((min_pair_distance(&mol) - expected).abs() < 1e-6);
}

#[test]
fn water_pdb_matches_xyz() {
    let pdb = parse_pdb(include_str!("../../../data/water.pdb")).unwrap();
    let xyz = parse_xyz(include_str!("../../../data/water.xyz")).unwrap();
    assert_eq!(pdb.formula(), "H2O");
    assert_eq!(pdb.len(), 3);
    let d_pdb = pdb.distance(0usize.into(), 1usize.into());
    let d_xyz = xyz.distance(0usize.into(), 1usize.into());
    assert!(
        (d_pdb - d_xyz).abs() < 1e-3,
        "O–H pdb {d_pdb} vs xyz {d_xyz}"
    );
}

#[test]
fn benzene_pdb_matches_xyz() {
    let pdb = parse_pdb(include_str!("../../../data/benzene.pdb")).unwrap();
    assert_eq!(pdb.formula(), "C6H6");
    assert_eq!(pdb.len(), 12);
    // Six ring bonds plus six C–H bonds, all from CONECT records.
    assert_eq!(pdb.n_bonds(), 12);
    // CONECT records carry the ring bonds, so the C–C distance is 1.39 Å.
    let d = pdb.distance(0usize.into(), 1usize.into());
    assert!((d - 1.39).abs() < 1e-3, "C–C = {d}");
}

/// Smallest pairwise distance between the Cartesian sites of a CIF (Å).
fn min_cif_pair_distance(cif: &tpt_chem_io::cif::CifStructure) -> f64 {
    let pos: Vec<_> = cif
        .atoms
        .iter()
        .map(|at| cif.cell.frac_to_cart(at.frac))
        .collect();
    let mut best = f64::INFINITY;
    for i in 0..pos.len() {
        for j in (i + 1)..pos.len() {
            best = best.min(pos[i].dist(pos[j]));
        }
    }
    best
}

#[test]
fn argon_fcc_cif_sample() {
    let cif = parse_cif(include_str!("../../../data/argon_fcc.cif")).unwrap();
    assert_eq!(cif.atoms.len(), 256);
    assert!((cif.cell.a - 21.04).abs() < 1e-9);
    let expected = 5.26 / 2f64.sqrt();
    assert!((min_cif_pair_distance(&cif) - expected).abs() < 1e-5);
}

#[test]
fn silicon_diamond_cif_sample() {
    let cif = parse_cif(include_str!("../../../data/silicon_diamond.cif")).unwrap();
    assert_eq!(cif.atoms.len(), 64);
    assert!((cif.cell.a - 10.862).abs() < 1e-9);
    let expected = 5.431 * 3f64.sqrt() / 4.0;
    assert!((min_cif_pair_distance(&cif) - expected).abs() < 1e-5);
}

#[test]
fn nacl_pdb_matches_xyz() {
    let pdb = parse_pdb(include_str!("../../../data/nacl_diatomic.pdb")).unwrap();
    assert_eq!(pdb.formula(), "ClNa");
    let d = pdb.distance(0usize.into(), 1usize.into());
    assert!((d - 2.36).abs() < 1e-3, "Na–Cl = {d}");
}
