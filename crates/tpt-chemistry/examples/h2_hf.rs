//! H₂ Hartree–Fock potential-energy scan (RHF/STO-3G).
//!
//! Run: `cargo run -p tpt-chemistry --features quantum --example h2_hf`

fn main() {
    use tpt_chemistry::core::molecule::Molecule;
    use tpt_chemistry::core::vec3::Vec3;
    use tpt_chemistry::quantum::hf::rhf_energy;

    println!("R(Angstrom)  E(RHF/STO-3G)   D_e vs atoms (Hartree)");
    println!("---------------------------------------------------");
    let mut best = (f64::INFINITY, 0.0);
    for &r in &[
        0.5, 0.6, 0.7, 0.7414, 0.8, 0.9, 1.0, 1.2, 1.4, 1.8, 2.4, 3.0, 5.0,
    ] {
        let mut m = Molecule::new("H2");
        m.add_atom::<1>(Vec3::new(0.0, 0.0, -r / 2.0));
        m.add_atom::<1>(Vec3::new(0.0, 0.0, r / 2.0));
        let e = rhf_energy(&m).expect("SCF converges on the H2 scan");
        // Dissociation limit: two isolated H atoms in STO-3G, E = -0.466582.
        let de = e - 2.0 * (-0.466_582_0);
        println!("{r:9.4}  {e:13.8}  {de:22.8}");
        if e < best.0 {
            best = (e, r);
        }
    }
    println!(
        "\nminimum: E = {:.8} Hartree at R = {:.4} Bohr \
         (literature STO-3G: -1.11750 at 1.346 Bohr)",
        best.0, best.1
    );
}
