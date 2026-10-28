//! Water RHF/STO-3G electron density and HOMO as Gaussian `.cube` files.
//!
//! Writes `water_density.cube` and `water_homo.cube` to the current
//! directory, and prints ∫ρ dV as a check (it should match the 10 electrons).
//! The 0.1 Bohr spacing resolves the oxygen 1s core, whose STO-3G exponents
//! reach about 130 Bohr⁻²; coarser grids miss electron density.
//!
//! Run: `cargo run -p tpt-chemistry --features core,quantum --example water_cube`

fn main() {
    use std::path::Path;
    use tpt_chemistry::core::molecule::Molecule;
    use tpt_chemistry::core::units::angstrom_to_bohr;
    use tpt_chemistry::core::vec3::Vec3;
    use tpt_chemistry::manifest::RunManifest;
    use tpt_chemistry::quantum::basis::sto3g_basis;
    use tpt_chemistry::quantum::cube::{density_on_grid, orbital_on_grid, write_cube, CubeGrid};
    use tpt_chemistry::quantum::hf::rhf;

    // Same Szabo–Ostlund geometry as the water_hf example.
    let r = 0.9894; // Å
    let half = (100.025f64 / 2.0).to_radians();
    let mut water = Molecule::new("water");
    water.add_atom::<8>(Vec3::ZERO);
    water.add_atom::<1>(Vec3::new(0.0, r * half.sin(), r * half.cos()));
    water.add_atom::<1>(Vec3::new(0.0, -r * half.sin(), r * half.cos()));

    let res = rhf(&water).expect("SCF converges");
    let bohr_atoms: Vec<(u8, Vec3)> = water
        .atoms()
        .map(|a| (a.z, a.pos * angstrom_to_bohr(1.0)))
        .collect();
    let shells = sto3g_basis(&bohr_atoms).expect("STO-3G covers H, O");

    // 10 electrons → five doubly occupied orbitals; the HOMO is index 4.
    let homo = 4;
    let positions: Vec<Vec3> = bohr_atoms.iter().map(|(_, p)| *p).collect();
    let grid = CubeGrid::enclosing(&positions, 4.0, 0.1);

    let rho = density_on_grid(&shells, &res.density, res.n_basis, &grid);
    let total: f64 = rho.iter().sum::<f64>() * grid.voxel_volume();
    println!("grid: {:?} points, ∫ρ dV = {total:.4} (expect 10)", grid.n);

    let psi = orbital_on_grid(&shells, &res.coefficients, res.n_basis, homo, &grid);
    println!("HOMO energy: {:.6} Hartree", res.orbital_energies[homo]);

    for (name, values, comment) in [
        ("water_density.cube", &rho, "RHF/STO-3G electron density"),
        ("water_homo.cube", &psi, "RHF/STO-3G HOMO"),
    ] {
        let text = write_cube(
            &water,
            &grid,
            values,
            &format!("{comment}\nwater, Bohr grid"),
        );
        std::fs::write(name, &text).expect("write cube file");
        RunManifest::new("water_cube", env!("CARGO_PKG_VERSION"))
            .param("basis", "STO-3G")
            .param("grid_margin_bohr", 4.0)
            .param("grid_spacing_bohr", 0.1)
            .output_bytes(name, text.as_bytes())
            .write_alongside(Path::new(name))
            .expect("write manifest");
        println!("wrote {name} ({} bytes)", text.len());
    }
}
