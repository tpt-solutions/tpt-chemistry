//! Water Hartree–Fock: energies, orbital spectrum, and a Molden-format
//! orbital dump for external visualizers.
//!
//! Run: `cargo run -p tpt-chemistry --features quantum --example water_hf`

fn main() {
    use tpt_chemistry::core::molecule::Molecule;
    use tpt_chemistry::core::units::angstrom_to_bohr;
    use tpt_chemistry::core::vec3::Vec3;
    use tpt_chemistry::quantum::basis::sto3g_basis;
    use tpt_chemistry::quantum::hf::rhf;
    use tpt_chemistry::quantum::molden::write_molden;

    // Szabo–Ostlund near-equilibrium geometry.
    let r = 0.9894; // Å
    let half = (100.025f64 / 2.0).to_radians();
    let mut water = Molecule::new("water");
    water.add_atom::<8>(Vec3::ZERO);
    water.add_atom::<1>(Vec3::new(0.0, r * half.sin(), r * half.cos()));
    water.add_atom::<1>(Vec3::new(0.0, -r * half.sin(), r * half.cos()));

    let res = rhf(&water).expect("SCF converges");
    println!("RHF/STO-3G water (r = {r} A, HOH = 100.025 deg)");
    println!("  total energy     : {:.8} Hartree", res.energy);
    println!("  electronic       : {:.8} Hartree", res.electronic_energy);
    println!("  nuclear repulsion: {:.8} Hartree", res.nuclear_energy);
    println!("  basis functions  : {}", res.n_basis);
    println!("  orbital energies (Hartree):");
    for (i, eps) in res.orbital_energies.iter().enumerate() {
        let occ = if i < res.orbital_energies.len() - 2 {
            "2.0"
        } else {
            "0.0"
        };
        println!("    {i}: {eps:12.6}   occ = {occ}");
    }

    // Molden dump for e.g. Jmol / chemcraft.
    let shells = sto3g_basis(
        &water
            .atoms()
            .map(|a| (a.z, a.pos * angstrom_to_bohr(1.0)))
            .collect::<Vec<_>>(),
    )
    .expect("STO-3G covers H, O");
    let molden = write_molden(&water, &shells, &res);
    let out = "water_sto3g.molden";
    std::fs::write(out, &molden).expect("write molden file");
    println!("\nwrote orbital file: {out} ({} bytes)", molden.len());
}
