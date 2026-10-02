use tpt_chem_core::vec3::Vec3;
use tpt_chem_crystal::bravais::Lattice;
use tpt_chem_crystal::xrd::{CrystalStructure, XrdSimulator};

fn main() {
    let lat = Lattice::cubic(5.431);
    let mut si = CrystalStructure::new(
        vec![14; 8],
        vec![
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.0, 0.5, 0.5),
            Vec3::new(0.5, 0.0, 0.5),
            Vec3::new(0.5, 0.5, 0.0),
            Vec3::new(0.25, 0.25, 0.25),
            Vec3::new(0.25, 0.75, 0.75),
            Vec3::new(0.75, 0.25, 0.75),
            Vec3::new(0.75, 0.75, 0.25),
        ],
    );
    si.b_factor = 0.5;
    let sim = XrdSimulator { wavelength: 1.5406, min_intensity: 1e-5, max_index: 6 };
    for p in sim.simulate(&lat, &si) {
        println!("({} {} {})  d={:.4} 2th={:.3} I={:.4}",
            p.miller.h, p.miller.k, p.miller.l, p.d, p.two_theta, p.intensity);
    }
}
