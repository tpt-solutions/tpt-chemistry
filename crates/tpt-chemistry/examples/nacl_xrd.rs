//! NaCl (rock salt) powder XRD: extinction rules in action — mixed
//! (h+k+l) odd reflections absent, fcc-allowed families present as CSV.
//!
//! Run: `cargo run -p tpt-chemistry --features crystal --example nacl_xrd`

use tpt_chemistry::core::vec3::Vec3;
use tpt_chemistry::crystal::bravais::Lattice;
use tpt_chemistry::crystal::xrd::{CrystalStructure, XrdSimulator};

fn main() {
    let a = 5.64; // Å
    let frac = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5], // Na fcc
        [0.5, 0.0, 0.0],
        [0.0, 0.5, 0.0],
        [0.0, 0.0, 0.5],
        [0.5, 0.5, 0.5], // Cl fcc
    ];
    let mut zs = vec![11u8; 4];
    zs.extend(vec![17u8; 4]);
    let structure = CrystalStructure::new(
        zs,
        frac.iter().map(|f| Vec3::new(f[0], f[1], f[2])).collect(),
    );
    let lat = Lattice::cubic(a);
    let peaks = XrdSimulator::default().simulate(&lat, &structure);

    println!("two_theta_deg,d_spacing_A,relative_intensity,miller");
    for p in &peaks {
        println!(
            "{:.4},{:.4},{:.6},({}{}{})",
            p.two_theta, p.d, p.intensity, p.miller.h, p.miller.k, p.miller.l
        );
    }
    let families: Vec<String> = peaks
        .iter()
        .map(|p| format!("({}{}{})", p.miller.h, p.miller.k, p.miller.l))
        .collect();
    println!("\n{families:?}");
    println!("experiment (Cu Kα): (111) 27.4°, (200) 31.7°, (220) 45.5°, (222) 56.5°");
}
