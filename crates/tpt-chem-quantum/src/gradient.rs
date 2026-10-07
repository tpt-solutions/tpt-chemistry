//! Analytic nuclear gradients of the RHF energy.
//!
//! Basis-function derivatives use the Gaussian identity
//! `d/dA_x [x^l e^(-a r^2)] = 2a * g(l+1) - l * g(l-1)`, so every
//! derivative integral reduces to ordinary integrals over *raw*
//! (unnormalized, single-primitive) shells with angular momentum `l +- 1`,
//! evaluated by the existing McMurchie-Davidson routines.
//!
//! The gradient is
//! `dE/dX = sum D.dH + 1/2 sum [D_ij D_kl - 1/2 D_ik D_jl] d(ij|kl)
//!          - sum W.dS + dE_nuc/dX`
//! with `W = 2 sum_occ eps C C^T`. Operator (Hellmann-Feynman) terms of the
//! nuclear attraction follow from translational invariance:
//! `d<a|V_C|b>/dC = -(d/dA + d/dB) <a|V_C|b>`; the fourth ERI centre and
//! the overlap/kinetic partner centre use the same invariance.
//!
//! Units: Hartree per Bohr.

use tpt_chem_core::molecule::Molecule;
use tpt_chem_core::potential::Potential;
use tpt_chem_core::units::{angstrom_to_bohr, HARTREE_KJ_PER_MOL};
use tpt_chem_core::vec3::Vec3;

use crate::basis::BasisSet;
use crate::gaussian::{cartesian_components, n_cartesian, primitive_normalization, Shell};
use crate::hf::{rhf_with_basis, HfError, HfResult};
use crate::integrals::{eri_shell, kinetic_shell, nuclear_shell, overlap_shell};

/// Single-primitive shell whose integrals are those of the unnormalized
/// Cartesian primitive `x^i y^j z^k e^(-alpha r^2)`.
fn raw_shell(l: u8, center: Vec3, alpha: f64) -> Shell {
    let mut sh = Shell::new(l, center, vec![(alpha, 1.0)]);
    sh.normalizations = cartesian_components(l)
        .iter()
        .map(|&c| 1.0 / primitive_normalization(alpha, c))
        .collect();
    sh
}

fn comp_axis(c: (u8, u8, u8), ax: usize) -> u8 {
    [c.0, c.1, c.2][ax]
}

fn comp_shift(c: (u8, u8, u8), ax: usize, delta: i8) -> (u8, u8, u8) {
    let mut a = [c.0, c.1, c.2];
    a[ax] = (i16::from(a[ax]) + i16::from(delta)) as u8;
    (a[0], a[1], a[2])
}

/// Derivative (along `ax`, w.r.t. the shell centre) of a flat tensor built
/// by `eval`, where `shell` occupies tensor position `pos` and `dims` are
/// the tensor dimensions with the original shell's component count at
/// `pos`. `eval` receives a raw shell to substitute at that position.
fn deriv_tensor(
    shell: &Shell,
    ax: usize,
    pos: usize,
    dims: [usize; 4],
    eval: &dyn Fn(&Shell) -> Vec<f64>,
) -> Vec<f64> {
    let outer: usize = dims[..pos].iter().product();
    let inner: usize = dims[pos + 1..].iter().product();
    let dpos = dims[pos];
    let mut out = vec![0.0; dims.iter().product()];
    let comps = cartesian_components(shell.l);
    let plus_comps = cartesian_components(shell.l + 1);
    let minus_comps = if shell.l > 0 {
        cartesian_components(shell.l - 1)
    } else {
        Vec::new()
    };
    let (dp, dm) = (
        n_cartesian(shell.l + 1),
        if shell.l > 0 {
            n_cartesian(shell.l - 1)
        } else {
            0
        },
    );
    for &(alpha, c) in &shell.primitives {
        let tp = eval(&raw_shell(shell.l + 1, shell.center, alpha));
        let tm = if shell.l > 0 {
            Some(eval(&raw_shell(shell.l - 1, shell.center, alpha)))
        } else {
            None
        };
        for (ic, &comp) in comps.iter().enumerate() {
            let base = c * shell.normalizations[ic] * primitive_normalization(alpha, comp);
            let ip = plus_comps
                .iter()
                .position(|&x| x == comp_shift(comp, ax, 1))
                .expect("plus component");
            let wp = base * 2.0 * alpha;
            let l_ax = comp_axis(comp, ax);
            let minus = if l_ax > 0 {
                let im = minus_comps
                    .iter()
                    .position(|&x| x == comp_shift(comp, ax, -1))
                    .expect("minus component");
                Some((im, -base * f64::from(l_ax)))
            } else {
                None
            };
            for o in 0..outer {
                for n in 0..inner {
                    out[(o * dpos + ic) * inner + n] += wp * tp[(o * dp + ip) * inner + n];
                    if let (Some((im, wm)), Some(tm)) = (minus, tm.as_ref()) {
                        out[(o * dpos + ic) * inner + n] += wm * tm[(o * dm + im) * inner + n];
                    }
                }
            }
        }
    }
    out
}

fn flat2(m: Vec<Vec<f64>>) -> Vec<f64> {
    m.into_iter().flatten().collect()
}

fn flat4(t: Vec<Vec<Vec<Vec<f64>>>>) -> Vec<f64> {
    t.into_iter().flatten().flatten().flatten().collect()
}

/// Derivatives of a two-centre matrix w.r.t. the centres of `a` and `b`:
/// `[ax]` -> flat `[ia][ib]` tensors.
fn d2(
    a: &Shell,
    b: &Shell,
    op: &dyn Fn(&Shell, &Shell) -> Vec<Vec<f64>>,
) -> ([Vec<f64>; 3], [Vec<f64>; 3]) {
    let dims = [n_cartesian(a.l), n_cartesian(b.l), 1, 1];
    let da = [0, 1, 2].map(|ax| deriv_tensor(a, ax, 0, dims, &|raw| flat2(op(raw, b))));
    let db = [0, 1, 2].map(|ax| deriv_tensor(b, ax, 1, dims, &|raw| flat2(op(a, raw))));
    (da, db)
}

/// RHF energy gradient (Hartree/Bohr) per atom, plus the converged SCF.
///
/// # Errors
/// See [`crate::hf::rhf_energy`].
pub fn rhf_gradient(mol: &Molecule) -> Result<(HfResult, Vec<Vec3>), HfError> {
    rhf_gradient_with_basis(mol, BasisSet::Sto3g)
}

/// [`rhf_gradient`] with an explicit basis set.
///
/// # Errors
/// See [`rhf_gradient`].
pub fn rhf_gradient_with_basis(
    mol: &Molecule,
    basis: BasisSet,
) -> Result<(HfResult, Vec<Vec3>), HfError> {
    let res = rhf_with_basis(mol, basis)?;
    let nb = res.n_basis;
    let atoms: Vec<(u8, Vec3)> = mol
        .atoms()
        .map(|a| (a.z, a.pos * angstrom_to_bohr(1.0)))
        .collect();
    let na = atoms.len();
    let nuclei: Vec<(f64, Vec3)> = atoms.iter().map(|&(z, r)| (f64::from(z), r)).collect();

    // Shells with their owning atom and basis offset.
    let mut shells: Vec<Shell> = Vec::new();
    let mut owner: Vec<usize> = Vec::new();
    for (ia, &(z, r)) in atoms.iter().enumerate() {
        for sh in basis
            .shells(z, r)
            .map_err(|_| HfError::UnsupportedElement)?
        {
            shells.push(sh);
            owner.push(ia);
        }
    }
    let mut offsets = Vec::new();
    let mut n = 0;
    for sh in &shells {
        offsets.push(n);
        n += n_cartesian(sh.l);
    }
    debug_assert_eq!(n, nb);

    let d = &res.density;
    let n_occ = {
        let ne: f64 = nuclei.iter().map(|&(z, _)| z).sum::<f64>() - f64::from(mol.formal_charge());
        (ne / 2.0).round() as usize
    };
    // Energy-weighted density W = 2 sum_occ eps_i C_mi C_ni.
    let mut w = vec![0.0; nb * nb];
    for m in 0..nb {
        for k in 0..nb {
            w[m * nb + k] = 2.0
                * (0..n_occ)
                    .map(|i| {
                        res.orbital_energies[i]
                            * res.coefficients[m * nb + i]
                            * res.coefficients[k * nb + i]
                    })
                    .sum::<f64>();
        }
    }

    let mut grad = vec![Vec3::ZERO; na];
    let add = |grad: &mut Vec<Vec3>, atom: usize, ax: usize, v: f64| {
        let mut g = grad[atom];
        g.set(ax, g.get(ax) + v);
        grad[atom] = g;
    };

    // One-electron terms.
    let ns = shells.len();
    for a in 0..ns {
        for b in 0..ns {
            let (sa, sb) = (&shells[a], &shells[b]);
            let (oa, ob) = (offsets[a], offsets[b]);
            let (dna, dnb) = (n_cartesian(sa.l), n_cartesian(sb.l));
            let contract = |t: &[f64], dens: &[f64]| -> f64 {
                let mut s = 0.0;
                for i in 0..dna {
                    for j in 0..dnb {
                        s += dens[(oa + i) * nb + ob + j] * t[i * dnb + j];
                    }
                }
                s
            };
            // Overlap: -W . dS.
            let (sda, sdb) = d2(sa, sb, &overlap_shell);
            // Kinetic.
            let (tda, tdb) = d2(sa, sb, &kinetic_shell);
            for ax in 0..3 {
                add(
                    &mut grad,
                    owner[a],
                    ax,
                    -contract(&sda[ax], &w) + contract(&tda[ax], d),
                );
                add(
                    &mut grad,
                    owner[b],
                    ax,
                    -contract(&sdb[ax], &w) + contract(&tdb[ax], d),
                );
            }
            // Nuclear attraction, per nucleus.
            for (ic, nuc) in nuclei.iter().enumerate() {
                let one = [*nuc];
                let (vda, vdb) = d2(sa, sb, &|x, y| nuclear_shell(x, y, &one));
                for ax in 0..3 {
                    let ga = contract(&vda[ax], d);
                    let gb = contract(&vdb[ax], d);
                    add(&mut grad, owner[a], ax, ga);
                    add(&mut grad, owner[b], ax, gb);
                    add(&mut grad, ic, ax, -(ga + gb));
                }
            }
        }
    }

    // Two-electron terms: 1/2 sum Gamma_ijkl d(ij|kl), Gamma = D_ij D_kl - 1/2 D_ik D_jl.
    for a in 0..ns {
        for b in 0..ns {
            for c in 0..ns {
                for dd in 0..ns {
                    let q = [a, b, c, dd];
                    let sh = [&shells[a], &shells[b], &shells[c], &shells[dd]];
                    let dims = [
                        n_cartesian(sh[0].l),
                        n_cartesian(sh[1].l),
                        n_cartesian(sh[2].l),
                        n_cartesian(sh[3].l),
                    ];
                    let off = [offsets[a], offsets[b], offsets[c], offsets[dd]];
                    // Skip quartets whose density weight is negligible.
                    let mut gam = vec![0.0; dims.iter().product()];
                    let mut idx = 0;
                    let mut any = false;
                    for i in 0..dims[0] {
                        for j in 0..dims[1] {
                            for k in 0..dims[2] {
                                for l in 0..dims[3] {
                                    let (gi, gj, gk, gl) =
                                        (off[0] + i, off[1] + j, off[2] + k, off[3] + l);
                                    let g = 0.5
                                        * (d[gi * nb + gj] * d[gk * nb + gl]
                                            - 0.5 * d[gi * nb + gk] * d[gj * nb + gl]);
                                    any |= g.abs() > 1e-14;
                                    gam[idx] = g;
                                    idx += 1;
                                }
                            }
                        }
                    }
                    if !any {
                        continue;
                    }
                    for ax in 0..3 {
                        let mut sum_abc = 0.0;
                        for pos in 0..3 {
                            let eval = |raw: &Shell| -> Vec<f64> {
                                let mut s = [sh[0], sh[1], sh[2], sh[3]];
                                s[pos] = raw;
                                flat4(eri_shell(s[0], s[1], s[2], s[3]))
                            };
                            let t = deriv_tensor(sh[pos], ax, pos, dims, &eval);
                            let g: f64 = t.iter().zip(&gam).map(|(x, y)| x * y).sum();
                            add(&mut grad, owner[q[pos]], ax, g);
                            sum_abc += g;
                        }
                        add(&mut grad, owner[dd], ax, -sum_abc);
                    }
                }
            }
        }
    }

    // Nuclear repulsion.
    for i in 0..na {
        for j in 0..na {
            if i != j {
                let r = nuclei[i].1 - nuclei[j].1;
                let f = -nuclei[i].0 * nuclei[j].0 / r.norm().powi(3);
                for ax in 0..3 {
                    add(&mut grad, i, ax, f * r.get(ax));
                }
            }
        }
    }
    Ok((res, grad))
}

/// Settings for [`optimize_geometry`].
#[derive(Clone, Copy, Debug)]
pub struct OptSettings {
    /// Maximum number of accepted steps.
    pub max_steps: usize,
    /// Convergence threshold on the largest gradient component
    /// (Hartree/Bohr).
    pub grad_tol: f64,
    /// Largest allowed displacement of any coordinate per step (Bohr).
    pub max_step: f64,
    /// Basis set used for energies and gradients.
    pub basis: BasisSet,
}

impl Default for OptSettings {
    fn default() -> Self {
        OptSettings {
            max_steps: 100,
            grad_tol: 3e-4,
            max_step: 0.3,
            basis: BasisSet::Sto3g,
        }
    }
}

/// Result of a geometry optimisation.
#[derive(Clone, Debug)]
pub struct OptResult {
    /// The optimised molecule (positions in Angstrom).
    pub molecule: Molecule,
    /// Final RHF energy (Hartree).
    pub energy: f64,
    /// Accepted steps taken.
    pub steps: usize,
    /// Whether `grad_tol` was reached.
    pub converged: bool,
    /// Final largest gradient component (Hartree/Bohr).
    pub max_gradient: f64,
}

fn with_positions(mol: &Molecule, x_bohr: &[f64]) -> Molecule {
    let inv = 1.0 / angstrom_to_bohr(1.0);
    let mut m = Molecule::new(mol.name());
    m.set_formal_charge(mol.formal_charge());
    let ids: Vec<_> = mol
        .atoms()
        .enumerate()
        .map(|(i, a)| {
            let p = Vec3::new(x_bohr[3 * i], x_bohr[3 * i + 1], x_bohr[3 * i + 2]) * inv;
            m.add_atom_raw(a.z, p)
        })
        .collect();
    for b in mol.bonds() {
        // Bonds copied from a valid molecule cannot fail.
        let _ = m.add_bond(ids[b.a.index()], ids[b.b.index()], b.order);
    }
    m
}

/// Minimise the RHF energy over all nuclear coordinates with BFGS on the
/// analytic gradient (backtracking line search, step-length cap).
///
/// # Errors
/// Any [`HfError`] from the underlying SCF.
pub fn optimize_geometry(mol: &Molecule, settings: OptSettings) -> Result<OptResult, HfError> {
    let n = 3 * mol.len();
    let mut x: Vec<f64> = mol
        .atoms()
        .flat_map(|a| {
            let p = a.pos * angstrom_to_bohr(1.0);
            [p.x, p.y, p.z]
        })
        .collect();
    let flat = |g: &[Vec3]| -> Vec<f64> { g.iter().flat_map(|v| [v.x, v.y, v.z]).collect() };
    let ident = || {
        let mut h = vec![0.0; n * n];
        for i in 0..n {
            h[i * n + i] = 1.0;
        }
        h
    };
    let mut cur = with_positions(mol, &x);
    let (res, g0) = rhf_gradient_with_basis(&cur, settings.basis)?;
    let (mut e, mut g) = (res.energy, flat(&g0));
    let mut hinv = ident();
    let mut steps = 0;
    let max_g = |g: &[f64]| g.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    while max_g(&g) > settings.grad_tol && steps < settings.max_steps {
        // Quasi-Newton direction.
        let mut dir: Vec<f64> = (0..n)
            .map(|i| -(0..n).map(|j| hinv[i * n + j] * g[j]).sum::<f64>())
            .collect();
        if dir.iter().zip(&g).map(|(d, gg)| d * gg).sum::<f64>() >= 0.0 {
            hinv = ident();
            dir = g.iter().map(|v| -v).collect();
        }
        let mut scale = (settings.max_step / max_g(&dir).max(1e-300)).min(1.0);
        let mut accepted = None;
        for _ in 0..10 {
            let xn: Vec<f64> = x.iter().zip(&dir).map(|(a, d)| a + scale * d).collect();
            let trial = with_positions(mol, &xn);
            let (r, gt) = rhf_gradient_with_basis(&trial, settings.basis)?;
            if r.energy < e + 1e-12 {
                accepted = Some((xn, r.energy, flat(&gt), trial));
                break;
            }
            scale *= 0.5;
        }
        let Some((xn, en, gn, trial)) = accepted else {
            if hinv != ident() {
                hinv = ident();
                continue;
            }
            break;
        };
        // BFGS inverse-Hessian update.
        let s_: Vec<f64> = xn.iter().zip(&x).map(|(a, b)| a - b).collect();
        let y: Vec<f64> = gn.iter().zip(&g).map(|(a, b)| a - b).collect();
        let sy: f64 = s_.iter().zip(&y).map(|(a, b)| a * b).sum();
        if sy > 1e-12 {
            let hy: Vec<f64> = (0..n)
                .map(|i| (0..n).map(|j| hinv[i * n + j] * y[j]).sum())
                .collect();
            let yhy: f64 = y.iter().zip(&hy).map(|(a, b)| a * b).sum();
            for i in 0..n {
                for j in 0..n {
                    hinv[i * n + j] += (1.0 + yhy / sy) * s_[i] * s_[j] / sy
                        - (hy[i] * s_[j] + s_[i] * hy[j]) / sy;
                }
            }
        }
        x = xn;
        e = en;
        g = gn;
        cur = trial;
        steps += 1;
    }
    let max_gradient = max_g(&g);
    Ok(OptResult {
        molecule: cur,
        energy: e,
        steps,
        converged: max_gradient <= settings.grad_tol,
        max_gradient,
    })
}

/// RHF/STO-3G as a [`Potential`]: energies in kJ/mol, forces in
/// kJ/mol/Angstrom from the analytic gradient. The molecule supplies the
/// elements, charge and bonds; positions come from each call.
#[derive(Clone, Debug)]
pub struct RhfPotential {
    /// Template molecule (elements and formal charge are used).
    pub molecule: Molecule,
    /// Basis set.
    pub basis: BasisSet,
}

impl RhfPotential {
    /// Wrap a molecule (STO-3G).
    pub fn new(molecule: Molecule) -> Self {
        RhfPotential {
            molecule,
            basis: BasisSet::Sto3g,
        }
    }

    /// Use another basis set.
    pub fn with_basis(mut self, basis: BasisSet) -> Self {
        self.basis = basis;
        self
    }
}

impl Potential for RhfPotential {
    type Error = HfError;

    fn n_atoms(&self) -> usize {
        self.molecule.len()
    }

    fn energy_forces(&mut self, positions: &[Vec3]) -> Result<(f64, Vec<Vec3>), HfError> {
        let bohr_per_a = angstrom_to_bohr(1.0);
        let x: Vec<f64> = positions
            .iter()
            .flat_map(|p| [p.x * bohr_per_a, p.y * bohr_per_a, p.z * bohr_per_a])
            .collect();
        let (res, grad) = rhf_gradient_with_basis(&with_positions(&self.molecule, &x), self.basis)?;
        let to_force = -HARTREE_KJ_PER_MOL * bohr_per_a;
        Ok((
            res.energy * HARTREE_KJ_PER_MOL,
            grad.into_iter().map(|g| g * to_force).collect(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hf::rhf_energy;

    fn water(shift: f64) -> Molecule {
        let mut m = Molecule::new("H2O");
        m.add_atom::<8>(Vec3::new(0.0, 0.0, 0.1173));
        m.add_atom::<1>(Vec3::new(0.0, 0.7572 + shift, -0.4692));
        m.add_atom::<1>(Vec3::new(0.0, -0.7572, -0.4692));
        m
    }

    fn fd_check(mol: &Molecule, coords: &[(usize, usize)]) {
        let (_, g) = rhf_gradient(mol).unwrap();
        let h = 1e-3; // Angstrom
        for &(atom, ax) in coords {
            {
                let shifted = |s: f64| {
                    let mut m = Molecule::new("x");
                    for (i, a) in mol.atoms().enumerate() {
                        let mut p = a.pos;
                        if i == atom {
                            p.set(ax, p.get(ax) + s);
                        }
                        match a.z {
                            1 => m.add_atom::<1>(p),
                            8 => m.add_atom::<8>(p),
                            _ => unreachable!(),
                        };
                    }
                    rhf_energy(&m).unwrap()
                };
                // Hartree/Angstrom -> Hartree/Bohr.
                let fd = (shifted(h) - shifted(-h)) / (2.0 * h) / angstrom_to_bohr(1.0);
                assert!(
                    (g[atom].get(ax) - fd).abs() < 2e-5,
                    "atom {atom} axis {ax}: analytic {} vs fd {fd}",
                    g[atom].get(ax)
                );
            }
        }
    }

    #[test]
    fn h2_gradient_matches_finite_difference() {
        let mut m = Molecule::new("H2");
        m.add_atom::<1>(Vec3::new(0.0, 0.0, -0.4));
        m.add_atom::<1>(Vec3::new(0.0, 0.1, 0.4));
        fd_check(&m, &[(0, 1), (0, 2), (1, 2)]);
    }

    #[test]
    fn water_gradient_matches_finite_difference() {
        // A subset of coordinates keeps the debug-build runtime sane.
        fd_check(&water(0.05), &[(0, 1), (0, 2), (1, 1), (1, 2)]);
    }

    #[test]
    fn gradient_sums_to_zero() {
        let (_, g) = rhf_gradient(&water(0.05)).unwrap();
        let s = g.iter().fold(Vec3::ZERO, |a, b| a + *b);
        assert!(s.norm() < 1e-8, "net force {}", s.norm());
    }

    #[test]
    fn h2_optimises_to_sto3g_bond_length() {
        let mut m = Molecule::new("H2");
        m.add_atom::<1>(Vec3::new(0.0, 0.0, -0.5));
        m.add_atom::<1>(Vec3::new(0.0, 0.0, 0.5));
        let r = optimize_geometry(&m, OptSettings::default()).unwrap();
        assert!(r.converged);
        let at: Vec<Vec3> = r.molecule.atoms().map(|a| a.pos).collect();
        // HF/STO-3G H2 equilibrium: 0.712 A, E = -1.11751 Eh.
        assert!(((at[1] - at[0]).norm() - 0.712).abs() < 0.003);
        assert!((r.energy - (-1.11751)).abs() < 1e-4, "{}", r.energy);
    }

    #[test]
    #[ignore = "slow in debug builds; run with --release -- --ignored"]
    fn water_optimises_to_sto3g_geometry() {
        let r = optimize_geometry(&water(0.15), OptSettings::default()).unwrap();
        assert!(r.converged, "max grad {}", r.max_gradient);
        let at: Vec<Vec3> = r.molecule.atoms().map(|a| a.pos).collect();
        let (d1, d2) = ((at[1] - at[0]).norm(), (at[2] - at[0]).norm());
        // HF/STO-3G water: r(OH) = 0.9894 A, angle 100.0 deg.
        assert!(
            (d1 - 0.9894).abs() < 0.005 && (d2 - 0.9894).abs() < 0.005,
            "{d1} {d2}"
        );
        let cos = (at[1] - at[0]).dot(at[2] - at[0]) / (d1 * d2);
        assert!((cos.acos().to_degrees() - 100.0).abs() < 0.5);
        assert!((r.energy - (-74.9659)).abs() < 2e-3, "{}", r.energy);
    }

    #[test]
    fn rhf_potential_drives_generic_minimiser() {
        use tpt_chem_core::potential::{minimize_bfgs, MinimizeSettings};
        let mut m = Molecule::new("H2");
        m.add_atom::<1>(Vec3::new(0.0, 0.0, -0.5));
        m.add_atom::<1>(Vec3::new(0.0, 0.0, 0.5));
        let mut pot = RhfPotential::new(m);
        let mut pos = [Vec3::new(0.0, 0.0, -0.5), Vec3::new(0.0, 0.0, 0.5)];
        let rep = minimize_bfgs(
            &mut pot,
            &mut pos,
            MinimizeSettings {
                force_tol: 0.5,
                ..Default::default()
            },
        )
        .unwrap();
        assert!(rep.converged, "{rep:?}");
        assert!(((pos[1] - pos[0]).norm() - 0.712).abs() < 0.003);
        assert!((rep.energy / HARTREE_KJ_PER_MOL - (-1.11751)).abs() < 1e-4);
    }

    #[test]
    fn h2_631g_gradient_matches_finite_difference() {
        use crate::hf::rhf_with_basis;
        let mk = |dz: f64| {
            let mut m = Molecule::new("H2");
            m.add_atom::<1>(Vec3::new(0.0, 0.0, -0.4 + dz));
            m.add_atom::<1>(Vec3::new(0.0, 0.1, 0.4));
            m
        };
        let (_, g) = rhf_gradient_with_basis(&mk(0.0), BasisSet::Pople631g).unwrap();
        let h = 1e-3;
        let e = |dz| rhf_with_basis(&mk(dz), BasisSet::Pople631g).unwrap().energy;
        let fd = (e(h) - e(-h)) / (2.0 * h) / angstrom_to_bohr(1.0);
        assert!((g[0].z - fd).abs() < 2e-5, "{} vs {fd}", g[0].z);
    }
}
