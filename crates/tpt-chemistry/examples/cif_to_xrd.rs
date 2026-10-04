//! End-to-end pipeline: CIF text → parsed structure → powder XRD.
//!
//! Run: `cargo run -p tpt-chemistry --features io,crystal --example cif_to_xrd`

use tpt_chemistry::core::vec3::Vec3;
use tpt_chemistry::crystal::bravais::Lattice;
use tpt_chemistry::crystal::xrd::{CrystalStructure, XrdSimulator};
use tpt_chemistry::io::cif::parse_cif;

const QUARTZ_CIF: &str = r#"
data_quartz
_cell_length_a 4.9134
_cell_length_b 4.9134
_cell_length_c 5.4052
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 120.0
_symmetry_space_group_name_H-M 'P 32 2 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 Si 0.4697 0.0000 0.0000
O1  O  0.4135 0.2669 0.1191
"#;

fn main() {
    let cif = parse_cif(QUARTZ_CIF).expect("parses the embedded quartz CIF");
    println!(
        "parsed '{}': {} sites, cell a = {}, c = {}",
        cif.name,
        cif.atoms.len(),
        cif.cell.a,
        cif.cell.c
    );

    // Convert the hexagonal cell to lattice vectors and simulate.
    let (a, b, c) = (cif.cell.a, cif.cell.b, cif.cell.c);
    let ga = cif.cell.gamma.to_radians();
    let lat = Lattice {
        a: Vec3::new(a, 0.0, 0.0),
        b: Vec3::new(b * ga.cos(), b * ga.sin(), 0.0),
        c: Vec3::new(0.0, 0.0, c),
    };
    let structure = CrystalStructure::new(
        cif.atoms.iter().map(|at| at.z).collect(),
        cif.atoms
            .iter()
            .map(|at| Vec3::new(at.frac[0], at.frac[1], at.frac[2]))
            .collect(),
    );
    let peaks = XrdSimulator::default().simulate(&lat, &structure);

    println!("\nα-quartz powder pattern (Cu Kα):");
    println!("two_theta_deg,d_spacing_A,relative_intensity,miller");
    for p in peaks.iter().take(10) {
        println!(
            "{:.4},{:.4},{:.6},({}{}{})",
            p.two_theta, p.d, p.intensity, p.miller.h, p.miller.k, p.miller.l
        );
    }
    println!("\n(experiment: strongest lines at 20.9° and 26.7°)");
}
