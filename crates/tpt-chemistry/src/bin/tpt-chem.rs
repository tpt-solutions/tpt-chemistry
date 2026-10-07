//! `tpt-chem` command-line front end.
//!
//! ```text
//! tpt-chem hf  <file.xyz> [--charge Q] [--mult M] [--basis sto-3g|6-31g] [--grad] [--opt] [--out opt.xyz]
//! tpt-chem xrd <file.cif> [--max N]
//! ```
//!
//! Build: `cargo run -p tpt-chemistry --features core,io,quantum,crystal --bin tpt-chem -- hf water.xyz`

use std::process::ExitCode;

use tpt_chemistry::core::vec3::Vec3;
use tpt_chemistry::crystal::bravais::Lattice;
use tpt_chemistry::crystal::xrd::{CrystalStructure, XrdSimulator};
use tpt_chemistry::io::cif::parse_cif;
use tpt_chemistry::io::xyz::{parse_xyz, write_xyz};
use tpt_chemistry::quantum::basis::BasisSet;
use tpt_chemistry::quantum::gradient::{optimize_geometry, rhf_gradient_with_basis, OptSettings};
use tpt_chemistry::quantum::hf::rhf_with_basis;
use tpt_chemistry::quantum::uhf::uhf_with_basis;

const USAGE: &str = "usage:
  tpt-chem hf  <file.xyz> [--charge Q] [--mult M] [--basis sto-3g|6-31g] [--grad] [--opt] [--out opt.xyz]
  tpt-chem xrd <file.cif> [--max N]
Energies are RHF/UHF STO-3G (Hartree); gradients in Hartree/Bohr; XRD for Cu K-alpha.";

struct Args {
    rest: Vec<String>,
}

impl Args {
    fn flag(&mut self, name: &str) -> bool {
        match self.rest.iter().position(|a| a == name) {
            Some(i) => {
                self.rest.remove(i);
                true
            }
            None => false,
        }
    }

    fn value(&mut self, name: &str) -> Result<Option<String>, String> {
        match self.rest.iter().position(|a| a == name) {
            Some(i) => {
                if i + 1 >= self.rest.len() {
                    return Err(format!("{name} needs a value"));
                }
                let v = self.rest.remove(i + 1);
                self.rest.remove(i);
                Ok(Some(v))
            }
            None => Ok(None),
        }
    }

    fn positional(&mut self) -> Result<String, String> {
        if self.rest.is_empty() {
            Err("missing input file".into())
        } else {
            Ok(self.rest.remove(0))
        }
    }
}

fn parse_num<T: std::str::FromStr>(v: Option<String>, default: T, name: &str) -> Result<T, String> {
    match v {
        None => Ok(default),
        Some(s) => s
            .parse()
            .map_err(|_| format!("invalid value for {name}: {s}")),
    }
}

fn cmd_hf(mut a: Args) -> Result<(), String> {
    let charge: i32 = parse_num(a.value("--charge")?, 0, "--charge")?;
    let mult: u32 = parse_num(a.value("--mult")?, 1, "--mult")?;
    let basis = match a.value("--basis")?.as_deref() {
        None | Some("sto-3g") => BasisSet::Sto3g,
        Some("6-31g") => BasisSet::Pople631g,
        Some(other) => return Err(format!("unknown basis '{other}' (sto-3g, 6-31g)")),
    };
    let label = if basis == BasisSet::Sto3g {
        "STO-3G"
    } else {
        "6-31G"
    };
    let out = a.value("--out")?;
    let want_grad = a.flag("--grad");
    let want_opt = a.flag("--opt");
    let path = a.positional()?;
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{path}: {e}"))?;
    let mut mol = parse_xyz(&text).map_err(|e| format!("{path}: {e:?}"))?;
    mol.set_formal_charge(charge);
    let closed = mult == 1;

    if want_opt {
        if !closed {
            return Err("--opt needs a closed-shell (--mult 1) calculation".into());
        }
        let r = optimize_geometry(
            &mol,
            OptSettings {
                basis,
                ..Default::default()
            },
        )
        .map_err(|e| e.to_string())?;
        println!(
            "optimisation: {} steps, converged = {}, max |grad| = {:.2e}",
            r.steps, r.converged, r.max_gradient
        );
        mol = r.molecule;
        if let Some(o) = &out {
            std::fs::write(
                o,
                write_xyz(&mol, &format!("optimised by tpt-chem (RHF/{label})")),
            )
            .map_err(|e| format!("{o}: {e}"))?;
            println!("wrote {o}");
        } else {
            print!(
                "{}",
                write_xyz(&mol, &format!("optimised by tpt-chem (RHF/{label})"))
            );
        }
    }

    if closed {
        let res = rhf_with_basis(&mol, basis).map_err(|e| e.to_string())?;
        println!("method: RHF/{label}");
        println!("energy_hartree: {:.10}", res.energy);
        println!("nuclear_repulsion_hartree: {:.10}", res.nuclear_energy);
        println!("n_basis: {}", res.n_basis);
        if want_grad {
            let (_, g) = rhf_gradient_with_basis(&mol, basis).map_err(|e| e.to_string())?;
            println!("gradient_hartree_per_bohr:");
            for (atom, v) in mol.atoms().zip(&g) {
                println!("  Z={} {:.8} {:.8} {:.8}", atom.z, v.x, v.y, v.z);
            }
        }
    } else {
        if want_grad {
            return Err("--grad is only available for closed-shell RHF".into());
        }
        let res = uhf_with_basis(&mol, mult, basis).map_err(|e| e.to_string())?;
        println!("method: UHF/{label} (multiplicity {mult})");
        println!("energy_hartree: {:.10}", res.energy);
        println!("s_squared: {:.6}", res.s_squared);
        println!("n_alpha: {}  n_beta: {}", res.n_alpha, res.n_beta);
    }
    Ok(())
}

fn cmd_xrd(mut a: Args) -> Result<(), String> {
    let max: usize = parse_num(a.value("--max")?, 20, "--max")?;
    let path = a.positional()?;
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{path}: {e}"))?;
    let cif = parse_cif(&text).map_err(|e| format!("{path}: {e:?}"))?;
    let (la, lb, lc) = (cif.cell.a, cif.cell.b, cif.cell.c);
    let lat = Lattice::triclinic(la, lb, lc, cif.cell.alpha, cif.cell.beta, cif.cell.gamma);
    let structure = CrystalStructure::new(
        cif.atoms.iter().map(|s| s.z).collect(),
        cif.atoms
            .iter()
            .map(|s| Vec3::new(s.frac[0], s.frac[1], s.frac[2]))
            .collect(),
    );
    let peaks = XrdSimulator::default().simulate(&lat, &structure);
    println!("two_theta_deg,d_spacing_A,relative_intensity,hkl");
    for p in peaks.iter().take(max) {
        println!(
            "{:.4},{:.4},{:.6},{} {} {}",
            p.two_theta, p.d, p.intensity, p.miller.h, p.miller.k, p.miller.l
        );
    }
    Ok(())
}

fn main() -> ExitCode {
    let mut argv: Vec<String> = std::env::args().skip(1).collect();
    if argv.is_empty() || argv[0] == "-h" || argv[0] == "--help" {
        println!("{USAGE}");
        return ExitCode::SUCCESS;
    }
    let cmd = argv.remove(0);
    let args = Args { rest: argv };
    let result = match cmd.as_str() {
        "hf" => cmd_hf(args),
        "xrd" => cmd_xrd(args),
        other => Err(format!("unknown subcommand '{other}'")),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}\n{USAGE}");
            ExitCode::FAILURE
        }
    }
}
