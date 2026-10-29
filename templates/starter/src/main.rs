//! Starter program: read a geometry from an XYZ file and print its RHF/STO-3G
//! energy. Edit this file to start your own calculation.
//!
//! Run: `cargo run -- path/to/molecule.xyz`

use std::error::Error;

use tpt_chemistry::io::xyz::parse_xyz;
use tpt_chemistry::quantum::hf::rhf;

fn main() -> Result<(), Box<dyn Error>> {
    let path = std::env::args()
        .nth(1)
        .ok_or("usage: {{project-name}} <file.xyz>")?;
    let text = std::fs::read_to_string(&path)?;
    // XYZ coordinates are Å, which is the crate-wide convention.
    let mol = parse_xyz(&text)?;
    let res = rhf(&mol)?;
    println!("molecule : {} ({})", mol.name(), mol.formula());
    println!("atoms    : {}", mol.len());
    println!("basis fns: {}", res.n_basis);
    println!("energy   : {:.8} Hartree", res.energy);
    Ok(())
}
