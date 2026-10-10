//! `tpt-chem` command-line front end.
//!
//! ```text
//! tpt-chem hf  <file.xyz> [--charge Q] [--mult M] [--basis sto-3g|6-31g] [--grad] [--opt] [--out opt.xyz]
//! tpt-chem xrd <file.cif> [--max N]
//! tpt-chem md  water <n_per_axis> [--model tip3p|spc] [--steps N] [--temp K]
//! tpt-chem md  <file.xyz> [--steps N] [--dt FS] [--temp K] [--minimise] [--out traj.xyz]
//! tpt-chem kinetics <network.txt> [--ssa] [--seed N] [--points N]
//! ```
//!
//! Build: `cargo run -p tpt-chemistry --features core,io,quantum,crystal,md,kinetics --bin tpt-chem -- hf water.xyz`

use std::process::ExitCode;

use tpt_chemistry::core::molecule::Molecule;
use tpt_chemistry::core::rng::Rng;
use tpt_chemistry::core::vec3::Vec3;
use tpt_chemistry::crystal::bravais::Lattice;
use tpt_chemistry::crystal::xrd::{CrystalStructure, XrdSimulator};
use tpt_chemistry::io::cif::parse_cif;
use tpt_chemistry::io::xyz::{parse_xyz, write_xyz};
use tpt_chemistry::kinetics::{network::ReactionNetwork, solver};
use tpt_chemistry::md::forces::ForceModel;
use tpt_chemistry::md::integrator::VelocityVerlet;
use tpt_chemistry::md::minimiser::minimise;
use tpt_chemistry::md::pme::PmeParams;
use tpt_chemistry::md::topology::{build_system, TopologySettings};
use tpt_chemistry::md::water::{constrained_temperature, water_box, SPC, TIP3P};
use tpt_chemistry::quantum::basis::BasisSet;
use tpt_chemistry::quantum::gradient::{optimize_geometry, rhf_gradient_with_basis, OptSettings};
use tpt_chemistry::quantum::hf::rhf_with_basis;
use tpt_chemistry::quantum::uhf::uhf_with_basis;

const USAGE: &str = "usage:
  tpt-chem hf  <file.xyz> [--charge Q] [--mult M] [--basis sto-3g|6-31g] [--grad] [--opt] [--out opt.xyz]
  tpt-chem xrd <file.cif> [--max N]
  tpt-chem md  water <n_per_axis> [--model tip3p|spc] [--steps N] [--temp K]
  tpt-chem md  <file.xyz> [--steps N] [--dt FS] [--temp K] [--minimise] [--out traj.xyz]
  tpt-chem kinetics <network.txt> [--ssa] [--seed N] [--points N]
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

fn cmd_md(mut a: Args) -> Result<(), String> {
    let steps: usize = parse_num(a.value("--steps")?, 1000, "--steps")?;
    let temp: f64 = parse_num(a.value("--temp")?, 300.0, "--temp")?;
    let model_name = a.value("--model")?;
    let dt_opt = a.value("--dt")?;
    let out = a.value("--out")?;
    let do_min = a.flag("--minimise");
    let first = a.positional()?;
    if first == "water" {
        let n: usize = a
            .positional()?
            .parse()
            .map_err(|_| "n_per_axis must be a positive integer".to_string())?;
        if n == 0 {
            return Err("n_per_axis must be at least 1".into());
        }
        let dt: f64 = parse_num(dt_opt, 1.0, "--dt")?;
        let model = match model_name.as_deref() {
            None | Some("tip3p") => &TIP3P,
            Some("spc") => &SPC,
            Some(o) => return Err(format!("unknown water model '{o}' (tip3p, spc)")),
        };
        let mut wb = water_box(model, n);
        let l = wb.system.box_.map(|b| b.length.x).unwrap_or(0.0);
        let fm = ForceModel::LjPlusPme {
            lj_cutoff: (0.5 * l).min(9.0),
            pme: PmeParams {
                alpha: 0.4,
                dims: [(l.ceil() as usize).next_power_of_two().max(8); 3],
                g_max: 0.0,
            },
        };
        wb.system.init_velocities(temp, &mut Rng::new(1));
        wb.constraints
            .rattle_velocities(&mut wb.system)
            .map_err(|e| e.to_string())?;
        println!(
            "# {} x {} molecules, box {l:.3} A, dt {dt} fs",
            model.name, wb.n_molecules
        );
        println!("step,temperature_K,potential_kJmol,total_kJmol");
        let equil = steps / 2;
        for step in 0..=steps {
            let pe = wb
                .constraints
                .step(&mut wb.system, &fm, dt)
                .map_err(|e| e.to_string())?;
            if step < equil && step % 10 == 0 {
                let s = (temp / constrained_temperature(&wb)).sqrt();
                for v in &mut wb.system.vel {
                    *v = *v * s;
                }
            }
            if step % (steps / 10).max(1) == 0 {
                println!(
                    "{step},{:.2},{pe:.3},{:.3}",
                    constrained_temperature(&wb),
                    pe + wb.system.kinetic_energy()
                );
            }
        }
        println!("# first half velocity-rescaled to {temp} K, second half NVE");
        return Ok(());
    }

    let text = std::fs::read_to_string(&first).map_err(|e| format!("{first}: {e}"))?;
    let mol = parse_xyz(&text).map_err(|e| format!("{first}: {e:?}"))?;
    let dt: f64 = parse_num(dt_opt, 0.5, "--dt")?;
    let mut sys = build_system(&mol, &TopologySettings::default());
    let fm = ForceModel::AllPairs { cutoff: None };
    if do_min {
        let rep = minimise(&mut sys, &fm, 5000, 1e-2, 0.01);
        println!(
            "# minimised: E {:.4} -> {:.4} kJ/mol in {} steps (converged: {})",
            rep.initial_energy, rep.final_energy, rep.steps, rep.converged
        );
    }
    sys.init_velocities(temp, &mut Rng::new(1));
    sys.remove_com_motion();
    sys.update_forces();
    let vv = VelocityVerlet::new(dt);
    let mut traj = String::new();
    println!("step,temperature_K,potential_kJmol,total_kJmol");
    for step in 0..=steps {
        if step > 0 {
            vv.step(&mut sys, None);
        }
        if step % (steps / 10).max(1) == 0 {
            println!(
                "{step},{:.2},{:.4},{:.4}",
                sys.temperature(),
                sys.potential_energy,
                sys.total_energy()
            );
        }
        if out.is_some() && step % (steps / 100).max(1) == 0 {
            let mut frame = Molecule::new(mol.name());
            for (atom, p) in mol.atoms().zip(&sys.pos) {
                frame.add_atom_raw(atom.z, *p);
            }
            traj.push_str(&write_xyz(&frame, &format!("step {step}")));
        }
    }
    if let Some(o) = out {
        std::fs::write(&o, traj).map_err(|e| format!("{o}: {e}"))?;
        println!("# wrote trajectory {o}");
    }
    Ok(())
}

type Side = Vec<(String, f64)>;

/// One side of a reaction, e.g. `2 A + B`.
fn parse_side(s: &str) -> Result<Side, String> {
    let mut out = Vec::new();
    for term in s
        .split('+')
        .map(str::trim)
        .filter(|t| !t.is_empty() && *t != "0")
    {
        let (n, name) = match term.find(|c: char| c.is_alphabetic() || c == '_') {
            Some(0) => (1.0, term),
            Some(i) => (
                term[..i]
                    .trim()
                    .parse::<f64>()
                    .map_err(|_| format!("bad stoichiometry in '{term}'"))?,
                term[i..].trim(),
            ),
            None => return Err(format!("bad term '{term}'")),
        };
        out.push((name.to_string(), n));
    }
    Ok(out)
}

/// Parse the plain-text reaction-network format:
///
/// ```text
/// species A B C
/// init A=100 B=0
/// t_end 10
/// reaction A + B -> C ; 0.01
/// ```
fn parse_network(text: &str) -> Result<(ReactionNetwork, Vec<String>, Vec<f64>, f64), String> {
    let mut species: Vec<String> = Vec::new();
    let mut init: Vec<(String, f64)> = Vec::new();
    let mut t_end = 1.0;
    let mut rxns: Vec<(Side, Side, f64)> = Vec::new();
    for (ln, raw) in text.lines().enumerate() {
        let line = raw.split('#').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        let err = |m: &str| format!("line {}: {m}", ln + 1);
        let (kw, rest) = line.split_once(char::is_whitespace).unwrap_or((line, ""));
        match kw {
            "species" => species.extend(rest.split_whitespace().map(String::from)),
            "init" => {
                for kv in rest.split_whitespace() {
                    let (k, v) = kv
                        .split_once('=')
                        .ok_or_else(|| err("expected NAME=VALUE"))?;
                    init.push((k.into(), v.parse().map_err(|_| err("bad number"))?));
                }
            }
            "t_end" => t_end = rest.trim().parse().map_err(|_| err("bad t_end"))?,
            "reaction" => {
                let (lhs_rhs, k) = rest
                    .split_once(';')
                    .ok_or_else(|| err("missing '; rate'"))?;
                let (l, r) = lhs_rhs
                    .split_once("->")
                    .ok_or_else(|| err("missing '->'"))?;
                let k: f64 = k.trim().parse().map_err(|_| err("bad rate"))?;
                rxns.push((
                    parse_side(l).map_err(|e| err(&e))?,
                    parse_side(r).map_err(|e| err(&e))?,
                    k,
                ));
            }
            other => return Err(err(&format!("unknown keyword '{other}'"))),
        }
    }
    if species.is_empty() {
        return Err("no 'species' line".into());
    }
    let mut y0 = vec![0.0; species.len()];
    for (n, v) in &init {
        let i = species
            .iter()
            .position(|s| s == n)
            .ok_or_else(|| format!("init: unknown species '{n}'"))?;
        y0[i] = *v;
    }
    let names: Vec<&str> = species.iter().map(String::as_str).collect();
    let mut b = ReactionNetwork::builder().species(&names);
    for (re, pr, k) in &rxns {
        for (n, _) in re.iter().chain(pr.iter()) {
            if !species.iter().any(|s| s == n) {
                return Err(format!("reaction uses undeclared species '{n}'"));
            }
        }
        let r: Vec<(&str, f64)> = re.iter().map(|(n, m)| (n.as_str(), *m)).collect();
        let p: Vec<(&str, f64)> = pr.iter().map(|(n, m)| (n.as_str(), *m)).collect();
        b = b.reaction(&r, &p, *k);
    }
    Ok((b.build(), species, y0, t_end))
}

fn cmd_kinetics(mut a: Args) -> Result<(), String> {
    let ssa = a.flag("--ssa");
    let seed: u64 = parse_num(a.value("--seed")?, 1, "--seed")?;
    let points: usize = parse_num(a.value("--points")?, 100, "--points")?;
    let path = a.positional()?;
    let text = std::fs::read_to_string(&path).map_err(|e| format!("{path}: {e}"))?;
    let (net, species, y0, t_end) = parse_network(&text).map_err(|e| format!("{path}: {e}"))?;
    println!("time,{}", species.join(","));
    if ssa {
        let mut rng = Rng::new(seed);
        for (t, y) in solver::ssa_events(&net, &y0, t_end, &mut rng) {
            let row: Vec<String> = y.iter().map(|v| format!("{v}")).collect();
            println!("{t:.6},{}", row.join(","));
        }
    } else {
        let traj = solver::integrate(&net, &y0, 0.0, t_end, points.max(1));
        let n = traj.len().max(2) - 1;
        for (i, y) in traj.iter().enumerate() {
            let row: Vec<String> = y.iter().map(|v| format!("{v:.8}")).collect();
            println!("{:.6},{}", t_end * i as f64 / n as f64, row.join(","));
        }
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
        "md" => cmd_md(args),
        "kinetics" => cmd_kinetics(args),
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn network_format_parses() {
        let (net, sp, y0, t) = parse_network(
            "species A B C  # comment\ninit A=10 B=5\nt_end 2.5\nreaction A + 2 B -> C ; 0.1\nreaction C -> A + 2 B ; 0.01\n",
        )
        .unwrap();
        assert_eq!(sp, ["A", "B", "C"]);
        assert_eq!(y0, [10.0, 5.0, 0.0]);
        assert_eq!(t, 2.5);
        assert_eq!(net.n_reactions(), 2);
        assert_eq!(net.reactions[0].stoichiometry, [-1.0, -2.0, 1.0]);
    }

    #[test]
    fn network_errors_are_reported() {
        assert!(parse_network("species A\nreaction A -> Z ; 1\n").is_err());
        assert!(parse_network("species A\nbogus 1\n").is_err());
        assert!(parse_network("reaction A -> B ; 1\n").is_err());
    }
}
