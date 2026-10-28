//! The sample structures in `data/` must parse and have the expected contents.

use tpt_chem_io::parse_xyz;

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
