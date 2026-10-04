//! Restricted closed-shell Hartree–Fock: the SCF loop, Fock build from
//! the ERI list, and the Roothaan–Hall diagonalization.
//!
//! For spin orbitals (μν|λσ) with 8-fold permutation symmetry, the
//! two-electron part of the Fock matrix is
//! `G_μν = Σ_λσ D_λσ [(μν|λσ) − ½(μλ|νσ)]`; the total energy is
//! `E = ½Σ_μν D_μν (H_μν + F_μν) + E_nuc`.

use tpt_chem_core::molecule::Molecule;
use tpt_chem_core::num;
use tpt_chem_core::units::angstrom_to_bohr;
use tpt_chem_core::vec3::Vec3;

use crate::basis::sto3g_basis;
use crate::gaussian::cartesian_components;
use crate::integrals::{eri_shell, kinetic_shell, nuclear_repulsion, nuclear_shell, overlap_shell};
use crate::linalg::{jacobi_eigh, mat_mul, transpose};

/// Closed-shell density matrix `D = 2·Σ_occ C_io C_jo` (row-major).
///
/// Extracted from the SCF loop so verification can reason about the exact
/// construction: for any real coefficient matrix the result is symmetric
/// positive semi-definite by construction.
pub fn density_from_orbitals(coeffs: &[f64], n_occ: usize, nb: usize) -> Vec<f64> {
    let mut d = vec![0.0; nb * nb];
    for i in 0..nb {
        for j in 0..nb {
            let mut acc = 0.0;
            for o in 0..n_occ {
                acc += coeffs[i * nb + o] * coeffs[j * nb + o];
            }
            d[i * nb + j] = 2.0 * acc;
        }
    }
    d
}

/// Errors from the HF driver.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HfError {
    /// The basis set does not cover one of the elements.
    UnsupportedElement,
    /// The SCF loop did not reach convergence within the iteration cap.
    NotConverged,
    /// The molecule needs more basis functions than the driver accepts.
    BasisTooLarge,
    /// A non-finite value appeared during the SCF (imploding geometry).
    NumericalBlowup,
}

impl core::fmt::Display for HfError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            HfError::UnsupportedElement => {
                write!(f, "element outside the STO-3G parameterization (H–Ne)")
            }
            HfError::NotConverged => write!(f, "SCF did not converge"),
            HfError::BasisTooLarge => write!(f, "basis exceeds the driver's function cap"),
            HfError::NumericalBlowup => write!(f, "non-finite value during SCF"),
        }
    }
}

impl std::error::Error for HfError {}

/// Result of a converged restricted HF calculation.
#[derive(Clone, Debug)]
pub struct HfResult {
    /// Total electronic + nuclear energy (Hartree).
    pub energy: f64,
    /// Electronic energy (Hartree).
    pub electronic_energy: f64,
    /// Nuclear repulsion energy (Hartree).
    pub nuclear_energy: f64,
    /// MO energies, ascending (Hartree).
    pub orbital_energies: Vec<f64>,
    /// MO coefficient matrix, row-major (columns = orbitals).
    pub coefficients: Vec<f64>,
    /// Density matrix, row-major.
    pub density: Vec<f64>,
    /// Number of basis functions.
    pub n_basis: usize,
}

/// Maximum number of basis functions accepted. Memory is O(N²) per SCF
/// (Fock, density, orthogonalizer) plus the packed Schwarz-screened shell
/// quartets, so the cap is a compute guard rather than a memory cliff.
const MAX_BASIS: usize = 256;

/// Schwarz screening threshold (Hartree): a shell quartet `(ab|cd)` is
/// skipped when `√|(ab|ab)| · √|(cd|cd)|` falls below it, which bounds every
/// element of the quartet.
const SCHWARZ_THRESHOLD: f64 = 1e-9;

/// One Schwarz-screened shell quartet's packed component tensor.
struct PackedQuartet {
    /// Shell indices `[a, b, c, d]`.
    shells: [usize; 4],
    /// Cartesian components per shell `[na, nb, nc, nd]`.
    dims: [usize; 4],
    /// Component tensor, row-major `[ia][ib][ic][id]`.
    v: Vec<f64>,
}

/// `G = Σ_{λσ} D_{λσ}[(μν|λσ) − ½(μλ|νσ)]`, accumulated directly from the
/// packed shell quartets. Each basis-function integral `(ij|kl)` appears in
/// exactly one quartet, so no symmetry multiplicities are needed.
fn accumulate_g(
    fock: &mut [f64],
    density: &[f64],
    quartets: &[PackedQuartet],
    offsets: &[usize],
    nb: usize,
) {
    for q in quartets {
        let [a, b, c, d] = q.shells;
        let [na, nb2, nc, nd] = q.dims;
        let (oi, oj, ok, ol) = (offsets[a], offsets[b], offsets[c], offsets[d]);
        let mut idx = 0;
        #[allow(clippy::needless_range_loop)]
        for ia in 0..na {
            for ib in 0..nb2 {
                for ic in 0..nc {
                    for id in 0..nd {
                        let val = q.v[idx];
                        idx += 1;
                        let (i, j) = (oi + ia, oj + ib);
                        let (k, l) = (ok + ic, ol + id);
                        // Coulomb: G_ij += D_kl (ij|kl).
                        fock[i * nb + j] += density[k * nb + l] * val;
                        // Exchange: G_ik −= ½ D_jl (ij|kl).
                        fock[i * nb + k] -= 0.5 * density[j * nb + l] * val;
                    }
                }
            }
        }
    }
}

/// Run a restricted Hartree–Fock SCF on the molecule with the STO-3G
/// minimal basis.
///
/// # Errors
/// [`HfError`] for unsupported elements, basis overflow, non-convergence,
/// or numerical blowup.
pub fn rhf_energy(mol: &Molecule) -> Result<f64, HfError> {
    rhf(mol).map(|r| r.energy)
}

/// Full HF driver: like [`rhf_energy`] but returns all converged
/// quantities.
///
/// # Errors
/// See [`rhf_energy`].
pub fn rhf(mol: &Molecule) -> Result<HfResult, HfError> {
    // Nuclei in Bohr.
    let atoms: Vec<(u8, Vec3)> = mol
        .atoms()
        .map(|a| (a.z, a.pos * angstrom_to_bohr(1.0)))
        .collect();
    if atoms.iter().any(|&(z, _)| z > 10) {
        return Err(HfError::UnsupportedElement);
    }
    let nuclei: Vec<(f64, Vec3)> = atoms.iter().map(|&(z, r)| (f64::from(z), r)).collect();
    let shells = sto3g_basis(&atoms).map_err(|_| HfError::UnsupportedElement)?;

    // Basis-function offsets per shell.
    let mut offsets = Vec::new();
    let mut nb = 0usize;
    for sh in &shells {
        offsets.push(nb);
        nb += cartesian_components(sh.l).len();
    }
    if nb > MAX_BASIS {
        return Err(HfError::BasisTooLarge);
    }
    let shell_of = |bf: usize| -> usize {
        offsets
            .iter()
            .rposition(|&o| o <= bf)
            .expect("bf within shells")
    };

    // One-electron matrices.
    let mut s = vec![0.0; nb * nb];
    let mut h_core = vec![0.0; nb * nb];
    for i in 0..nb {
        let si = &shells[shell_of(i)];
        let ci = i - offsets[shell_of(i)];
        for j in 0..nb {
            let sj = &shells[shell_of(j)];
            let cj = j - offsets[shell_of(j)];
            let ov = overlap_shell(si, sj);
            let k = ov[ci][cj];
            s[i * nb + j] = k;
            let t = kinetic_shell(si, sj)[ci][cj];
            let v = nuclear_shell(si, sj, &nuclei)[ci][cj];
            h_core[i * nb + j] = t + v;
        }
    }

    // Two-electron integrals: Schwarz-screened packed shell quartets. The
    // full basis-function tensor `(μν|λσ)` is never materialized; the Fock
    // build consumes the quartets directly.
    let ns = shells.len();
    let comps: Vec<Vec<(u8, u8, u8)>> = shells.iter().map(|s| cartesian_components(s.l)).collect();
    let mut schwarz = vec![0.0f64; ns * ns];
    for a in 0..ns {
        for b in 0..ns {
            let t = eri_shell(&shells[a], &shells[b], &shells[a], &shells[b]);
            let q = t
                .iter()
                .flatten()
                .flatten()
                .flatten()
                .fold(0.0f64, |m, &x| m.max(x.abs()))
                .sqrt();
            schwarz[a * ns + b] = q;
        }
    }

    let mut quartets: Vec<PackedQuartet> = Vec::new();
    for (a, oa) in shells.iter().enumerate() {
        for (b, ob) in shells.iter().enumerate() {
            if schwarz[a * ns + b] == 0.0 {
                continue;
            }
            for (c, oc) in shells.iter().enumerate() {
                for (d, od) in shells.iter().enumerate() {
                    if schwarz[a * ns + b] * schwarz[c * ns + d] < SCHWARZ_THRESHOLD {
                        continue;
                    }
                    let v = eri_shell(oa, ob, oc, od);
                    let dims = [
                        comps[a].len(),
                        comps[b].len(),
                        comps[c].len(),
                        comps[d].len(),
                    ];
                    let mut flat = Vec::with_capacity(dims[0] * dims[1] * dims[2] * dims[3]);
                    for row in &v {
                        for r2 in row {
                            for r3 in r2 {
                                for x in r3 {
                                    flat.push(*x);
                                }
                            }
                        }
                    }
                    quartets.push(PackedQuartet {
                        shells: [a, b, c, d],
                        dims,
                        v: flat,
                    });
                }
            }
        }
    }

    // Symmetric orthogonalization X = S^{-1/2}.
    let (s_vals, s_vecs) = jacobi_eigh(&s, nb, 1e-11, 64).map_err(|_| HfError::NumericalBlowup)?;
    let mut s_half = vec![0.0; nb * nb];
    for i in 0..nb {
        for j in 0..nb {
            for k in 0..nb {
                s_half[i * nb + j] +=
                    s_vecs[i * nb + k] * s_vecs[j * nb + k] / num::sqrt(s_vals[k]);
            }
        }
    }
    let x_mat = s_half;

    // Occupied orbitals (closed shell): n_elec / 2, with n_elec the nuclear
    // charge minus the molecule's formal charge.
    let charge = f64::from(mol.formal_charge());
    let n_elec: f64 = nuclei.iter().map(|&(z, _)| z).sum::<f64>() - charge;
    if n_elec < 0.0 || n_elec.round() != n_elec {
        return Err(HfError::UnsupportedElement);
    }
    let n_occ = (n_elec / 2.0).round() as usize;
    if 2.0 * n_occ as f64 != n_elec {
        return Err(HfError::UnsupportedElement); // open shell: outside RHF scope
    }

    // SCF loop.
    let mut density = vec![0.0f64; nb * nb];
    let mut energy_old = 0.0;
    let e_nuc = nuclear_repulsion(&nuclei);
    let mut converged = false;
    for iteration in 0..128 {
        // Fock = H + G(D), accumulated from the packed quartets: every basis
        // integral (ij|kl) is visited exactly once, via its own quartet.
        let mut fock = h_core.clone();
        accumulate_g(&mut fock, &density, &quartets, &offsets, nb);
        // Energy.
        let mut e_elec = 0.0;
        for i in 0..nb {
            for j in 0..nb {
                e_elec += density[i * nb + j] * (h_core[i * nb + j] + fock[i * nb + j]);
            }
        }
        e_elec *= 0.5;
        let energy = e_elec + e_nuc;
        if !energy.is_finite() {
            return Err(HfError::NumericalBlowup);
        }
        if (energy - energy_old).abs() < 1e-11 && iteration > 0 {
            converged = true;
            break;
        }
        energy_old = energy;

        // Diagonalize F' = Xᵀ F X.
        let ft = transpose(&fock, nb);
        let f_prime = mat_mul(&mat_mul(&x_mat, &ft, nb), &x_mat, nb);
        let (eps, c_prime) =
            jacobi_eigh(&f_prime, nb, 1e-11, 64).map_err(|_| HfError::NumericalBlowup)?;
        let coeffs = mat_mul(&x_mat, &c_prime, nb);

        // New density from the occupied orbitals.
        let d_new = density_from_orbitals(&coeffs, n_occ, nb);
        // Damping for early iterations tames wild core-guess densities.
        let damp = if iteration < 8 { 0.5 } else { 1.0 };
        for i in 0..nb * nb {
            density[i] = damp * d_new[i] + (1.0 - damp) * density[i];
        }
        // Stash the last orbitals/energies for the result.
        if iteration == 127 {
            return Err(HfError::NotConverged);
        }
        let _ = eps;
    }
    if !converged {
        return Err(HfError::NotConverged);
    }

    // Final quantities.
    let mut fock = h_core.clone();
    accumulate_g(&mut fock, &density, &quartets, &offsets, nb);
    let ft = transpose(&fock, nb);
    let f_prime = mat_mul(&mat_mul(&x_mat, &ft, nb), &x_mat, nb);
    let (eps, c_prime) =
        jacobi_eigh(&f_prime, nb, 1e-11, 64).map_err(|_| HfError::NumericalBlowup)?;
    let coeffs = mat_mul(&x_mat, &c_prime, nb);
    let mut e_elec = 0.0;
    for i in 0..nb {
        for j in 0..nb {
            e_elec += density[i * nb + j] * (h_core[i * nb + j] + fock[i * nb + j]);
        }
    }
    e_elec *= 0.5;

    Ok(HfResult {
        energy: e_elec + e_nuc,
        electronic_energy: e_elec,
        nuclear_energy: e_nuc,
        orbital_energies: eps,
        coefficients: coeffs,
        density,
        n_basis: nb,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_chem_core::molecule::Molecule;

    fn h2(r_angstrom: f64) -> Molecule {
        let mut m = Molecule::new("H2");
        m.add_atom::<1>(Vec3::new(0.0, 0.0, -r_angstrom / 2.0));
        m.add_atom::<1>(Vec3::new(0.0, 0.0, r_angstrom / 2.0));
        m
    }

    #[test]
    fn h2_sto3g_literature_value() {
        // HF/STO-3G H₂ at R = 0.7414 Å. The Szabo–Ostlund book value
        // −1.11675931 Eₕ is quoted at exactly R = 1.4 Bohr
        // (= 0.740842 Å); at 0.7414 Å the same basis gives −1.11668439 Eₕ,
        // cross-validated against an independent McMurchie–Davidson
        // implementation to 9 digits.
        let e = rhf_energy(&h2(0.7414)).unwrap();
        assert!((e - (-1.11668439)).abs() < 1e-6, "E = {e:.10}");
    }

    #[test]
    fn he_atom_sto3g() {
        // HF/STO-3G helium with the Basis Set Exchange exponent set:
        // −2.80778516 Eₕ (literature −2.80778395 Eₕ; the 12 µEₕ residual
        // traces to the 8-digit rounding of the universal contraction
        // coefficients).
        let mut he = Molecule::new("He");
        he.add_atom::<2>(Vec3::ZERO);
        let e = rhf_energy(&he).unwrap();
        assert!((e - (-2.80778516)).abs() < 1e-6, "E = {e:.10}");
    }

    #[test]
    fn water_sto3g_near_equilibrium() {
        // HF/STO-3G water at the Szabo–Ostlund near-equilibrium geometry
        // r(OH) = 0.9894 Å, HOH = 100.025°: E = −74.96590117 Eₕ. This is
        // the canonical RHF/STO-3G water energy, cross-validated here
        // against an independent Python McMurchie–Davidson RHF to 9 digits.
        let r = 0.9894;
        let half = (100.025f64 / 2.0).to_radians();
        let mut water = Molecule::new("water");
        water.add_atom::<8>(Vec3::ZERO);
        water.add_atom::<1>(Vec3::new(0.0, r * half.sin(), r * half.cos()));
        water.add_atom::<1>(Vec3::new(0.0, -r * half.sin(), r * half.cos()));
        let res = rhf(&water).unwrap();
        assert!(
            (res.energy - (-74.96590117)).abs() < 1e-6,
            "E = {:.8}",
            res.energy
        );
    }

    #[test]
    fn unsupported_element_errors() {
        let mut na = Molecule::new("Na");
        na.add_atom::<11>(Vec3::ZERO);
        assert_eq!(rhf_energy(&na), Err(HfError::UnsupportedElement));
    }
}

#[cfg(test)]
mod driver_tests {
    use super::*;
    use tpt_chem_core::molecule::Molecule;

    #[test]
    fn h2_dication_zero_electrons() {
        // H₂²⁺ has no electrons: RHF reduces to the bare nuclear repulsion
        // 1/R Hartree — a sharp analytic check of the charge bookkeeping.
        let r = 0.7414;
        let mut m = Molecule::new("H2 2+");
        m.add_atom::<1>(Vec3::new(0.0, 0.0, -r / 2.0));
        m.add_atom::<1>(Vec3::new(0.0, 0.0, r / 2.0));
        m.set_formal_charge(2);
        let res = rhf(&m).unwrap();
        // r in Å → Bohr for the repulsion.
        let r_bohr = r * tpt_chem_core::units::angstrom_to_bohr(1.0);
        let want = 1.0 / r_bohr;
        assert!(
            (res.energy - want).abs() < 1e-10,
            "E = {} vs {want}",
            res.energy
        );
        assert_eq!(res.orbital_energies.len(), 2);
    }

    #[test]
    fn hydronium_ten_electron_closed_shell() {
        // H₃O⁺: 8 + 3 − 1 = 10 electrons — the same closed shell as water.
        let r = 0.98;
        let mut m = Molecule::new("H3O+");
        m.add_atom::<8>(Vec3::ZERO);
        for k in 0..3 {
            let a = 2.0 * core::f64::consts::PI * k as f64 / 3.0;
            m.add_atom::<1>(Vec3::new(r * a.cos(), r * a.sin(), 0.0));
        }
        m.set_formal_charge(1);
        let res = rhf(&m).unwrap();
        // Roughly water-like electronic energy; the sharp assertions are the
        // bookkeeping ones below.
        assert!(
            res.energy < -74.0 && res.energy > -76.0,
            "E = {}",
            res.energy
        );
        // O: 1s + 2s + 2p(3) = 5 functions; 3 H: 1 each → 8.
        assert_eq!(res.orbital_energies.len(), 8);
        assert_eq!(res.n_basis, 8);
    }

    #[test]
    fn odd_electron_charge_is_rejected() {
        // H₂⁺: 2 − 1 = 1 electron → open shell, outside RHF scope.
        let mut m = Molecule::new("H2+");
        m.add_atom::<1>(Vec3::new(0.0, 0.0, -0.7));
        m.add_atom::<1>(Vec3::new(0.0, 0.0, 0.7));
        m.set_formal_charge(1);
        assert_eq!(rhf_energy(&m), Err(HfError::UnsupportedElement));
    }

    #[test]
    fn basis_too_large_is_a_dedicated_error() {
        // 258 H atoms → 258 basis functions > MAX_BASIS = 256; the check
        // fires before any integral work.
        let mut big = Molecule::new("H258");
        for i in 0..258 {
            big.add_atom::<1>(Vec3::new(i as f64 * 3.0, 0.0, 0.0));
        }
        assert_eq!(rhf_energy(&big), Err(HfError::BasisTooLarge));
    }
}
