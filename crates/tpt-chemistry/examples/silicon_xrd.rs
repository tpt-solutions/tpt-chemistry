//! Silicon diamond-structure powder XRD: the classic peak pattern as CSV.
//!
//! Run: `cargo run -p tpt-chemistry --features crystal --example silicon_xrd`

use tpt_chemistry::core::vec3::Vec3;
use tpt_chemistry::crystal::bravais::Lattice;
use tpt_chemistry::crystal::xrd::{CrystalStructure, XrdSimulator};

fn main() {
    // Diamond structure: fcc Si + the (1/4,1/4,1/4) basis.
    let a = 5.431; // Å
    let frac = [
        [0.0, 0.0, 0.0],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
        [0.25, 0.25, 0.25],
        [0.75, 0.75, 0.25],
        [0.75, 0.25, 0.75],
        [0.25, 0.75, 0.75],
    ];
    let structure = CrystalStructure::new(
        vec![14u8; frac.len()],
        frac.iter().map(|f| Vec3::new(f[0], f[1], f[2])).collect(),
    );
    let lat = Lattice::cubic(a);
    let sim = XrdSimulator::default(); // Cu Kα
    let peaks = sim.simulate(&lat, &structure);

    println!("two_theta_deg,d_spacing_A,relative_intensity,miller");
    for p in &peaks {
        println!(
            "{:.4},{:.4},{:.6},({}{}{})",
            p.two_theta, p.d, p.intensity, p.miller.h, p.miller.k, p.miller.l
        );
    }
    println!(
        "\n{} reflections; strongest at 2θ = {:.3}° \
         (experiment: (111) at 28.44°)",
        peaks.len(),
        peaks[0].two_theta
    );
}
