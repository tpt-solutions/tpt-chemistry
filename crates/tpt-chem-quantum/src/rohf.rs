//! Restricted open-shell Hartree–Fock (ROHF): one orbital set shared by the
//! alpha and beta spin, with doubly occupied, singly occupied ("open") and
//! virtual orbital classes.
//!
//! Each iteration builds the spin Fock matrices `F^α` and `F^β` from the
//! current densities, then combines them into a single *effective* Fock
//! matrix in the MO basis of the current orbitals (Roothaan's rules):
//!
//! | block          | closed | open | virtual |
//! |----------------|--------|------|---------|
//! | **closed**     | `F̄`    | `F^β`| `F̄`    |
//! | **open**       | `F^β`  | `F̄` | `F^α`   |
//! | **virtual**    | `F̄`    | `F^α`| `F̄`    |
//!
//! with `F̄ = (F^α + F^β) / 2`. The effective matrix is diagonalized and the
//! new orbitals are its eigenvectors, so the lowest `n_β` are doubly occupied
//! and the next `n_α − n_β` singly occupied in each following iteration.
//!
//! For a closed-shell system ROHF reduces to RHF. For a single open-shell
//! electron ROHF equals UHF. In general `E_ROHF ≥ E_UHF`, because ROHF
//! restricts the spatial orbitals to be shared by both spins.

use tpt_chem_core::molecule::Molecule;

use crate::basis::BasisSet;
use crate::hf::{scf_setup, HfError, ScfSetup};
use crate::linalg::{jacobi_eigh, mat_mul, transpose};
use crate::uhf::{accumulate_jk, diag, spin_density, spin_occupations};

/// Result of a converged ROHF calculation.
#[derive(Clone, Debug)]
pub struct RohfResult {
    /// Total energy (Hartree).
    pub energy: f64,
    /// Nuclear repulsion energy (Hartree).
    pub nuclear_energy: f64,
    /// Orbital energies (eigenvalues of the effective Fock matrix), ascending.
    pub orbital_energies: Vec<f64>,
    /// MO coefficients, row-major (columns = orbitals), shared by both spins.
    pub coefficients: Vec<f64>,
    /// Alpha density matrix (no factor 2).
    pub density_alpha: Vec<f64>,
    /// Beta density matrix (no factor 2).
    pub density_beta: Vec<f64>,
    /// Number of alpha electrons.
    pub n_alpha: usize,
    /// Number of beta electrons.
    pub n_beta: usize,
    /// Number of basis functions.
    pub n_basis: usize,
}

/// Roothaan effective Fock matrix in the MO basis. `n_closed` orbitals are
/// doubly occupied and the orbitals up to `n_occ` are singly occupied.
fn effective_fock(fa: &[f64], fb: &[f64], nb: usize, n_closed: usize, n_occ: usize) -> Vec<f64> {
    // 0 = closed, 1 = open, 2 = virtual.
    let class = |p: usize| {
        if p < n_closed {
            0
        } else if p < n_occ {
            1
        } else {
            2
        }
    };
    let mut f = vec![0.0; nb * nb];
    for p in 0..nb {
        for q in 0..nb {
            let (lo, hi) = {
                let (cp, cq) = (class(p), class(q));
                (cp.min(cq), cp.max(cq))
            };
            f[p * nb + q] = match (lo, hi) {
                (0, 1) => fb[p * nb + q],
                (1, 2) => fa[p * nb + q],
                _ => 0.5 * (fa[p * nb + q] + fb[p * nb + q]),
            };
        }
    }
    f
}

/// ROHF energy with the given spin multiplicity `2S + 1` (STO-3G basis).
///
/// # Errors
/// As [`crate::uhf::uhf`]: an incompatible multiplicity is rejected.
pub fn rohf(mol: &Molecule, multiplicity: u32) -> Result<RohfResult, HfError> {
    rohf_with_basis(mol, multiplicity, BasisSet::Sto3g)
}

/// [`rohf`] with an explicit basis set.
///
/// # Errors
/// See [`rohf`].
pub fn rohf_with_basis(
    mol: &Molecule,
    multiplicity: u32,
    basis: BasisSet,
) -> Result<RohfResult, HfError> {
    let ScfSetup {
        nb,
        offsets,
        h_core,
        quartets,
        x_mat,
        e_nuc,
        n_elec,
        ..
    } = scf_setup(mol, basis)?;
    let (n_alpha, n_beta) = spin_occupations(n_elec, multiplicity, nb)?;

    // Start from the core-Hamiltonian guess.
    let (_, mut c) = diag(&h_core, &x_mat, nb)?;
    let mut e_old = 0.0;
    let mut last = None;
    let mut converged = false;
    let (mut da_old, mut db_old) = (vec![0.0; nb * nb], vec![0.0; nb * nb]);
    for it in 0..256 {
        let da = spin_density(&c, n_alpha, nb);
        let db = spin_density(&c, n_beta, nb);
        let dt: Vec<f64> = da.iter().zip(&db).map(|(a, b)| a + b).collect();
        let (j, ka, kb) = accumulate_jk(&quartets, &offsets, nb, &dt, &da, &db);
        let fa: Vec<f64> = (0..nb * nb).map(|i| h_core[i] + j[i] - ka[i]).collect();
        let fb: Vec<f64> = (0..nb * nb).map(|i| h_core[i] + j[i] - kb[i]).collect();
        let e_elec: f64 = 0.5
            * (0..nb * nb)
                .map(|i| dt[i] * h_core[i] + da[i] * fa[i] + db[i] * fb[i])
                .sum::<f64>();
        let energy = e_elec + e_nuc;
        if !energy.is_finite() {
            return Err(HfError::NumericalBlowup);
        }

        // The orbitals are S-orthonormal, so `Cᵀ F C` is the MO-basis Fock.
        let ct = transpose(&c, nb);
        let fa_mo = mat_mul(&mat_mul(&ct, &fa, nb), &c, nb);
        let fb_mo = mat_mul(&mat_mul(&ct, &fb, nb), &c, nb);
        let feff = effective_fock(&fa_mo, &fb_mo, nb, n_beta, n_alpha);
        let (eps, u) = jacobi_eigh(&feff, nb, 1e-11, 64).map_err(|_| HfError::NumericalBlowup)?;
        let c_new = mat_mul(&c, &u, nb);

        let da_new = spin_density(&c_new, n_alpha, nb);
        let db_new = spin_density(&c_new, n_beta, nb);
        let dd = da_new
            .iter()
            .zip(&da_old)
            .chain(db_new.iter().zip(&db_old))
            .map(|(n, o)| (n - o).abs())
            .fold(0.0f64, f64::max);
        last = Some((eps, c_new.clone(), energy));
        if it > 0 && (energy - e_old).abs() < 1e-10 && dd < 1e-7 {
            converged = true;
            break;
        }
        e_old = energy;
        da_old = da;
        db_old = db;
        c = c_new;
    }
    if !converged {
        return Err(HfError::NotConverged);
    }
    let (orbital_energies, coefficients, energy) = last.ok_or(HfError::NotConverged)?;
    Ok(RohfResult {
        energy,
        nuclear_energy: e_nuc,
        orbital_energies,
        density_alpha: spin_density(&coefficients, n_alpha, nb),
        density_beta: spin_density(&coefficients, n_beta, nb),
        coefficients,
        n_alpha,
        n_beta,
        n_basis: nb,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hf::rhf_energy;
    use crate::uhf::uhf;
    use tpt_chem_core::vec3::Vec3;

    fn h2(r: f64) -> Molecule {
        let mut m = Molecule::new("H2");
        m.add_atom::<1>(Vec3::new(0.0, 0.0, -r / 2.0));
        m.add_atom::<1>(Vec3::new(0.0, 0.0, r / 2.0));
        m
    }

    #[test]
    fn hydrogen_atom_matches_uhf() {
        // One open-shell electron: ROHF and UHF are the same calculation.
        let mut h = Molecule::new("H");
        h.add_atom::<1>(Vec3::ZERO);
        let r = rohf(&h, 2).unwrap();
        assert!((r.energy - (-0.466_581_85)).abs() < 1e-6, "{}", r.energy);
    }

    #[test]
    fn closed_shell_matches_rhf() {
        let m = h2(0.74);
        let r = rohf(&m, 1).unwrap();
        assert!((r.energy - rhf_energy(&m).unwrap()).abs() < 1e-8);
    }

    #[test]
    fn rohf_is_above_uhf_for_triplet() {
        // The restricted orbitals cannot be lower than the unrestricted ones.
        let m = h2(0.74);
        let r = rohf(&m, 3).unwrap();
        let u = uhf(&m, 3).unwrap();
        assert!(r.energy >= u.energy - 1e-8, "{} < {}", r.energy, u.energy);
        assert_eq!(r.n_alpha, 2);
        assert_eq!(r.n_beta, 0);
    }

    #[test]
    fn oh_radical_is_above_uhf_and_close_to_it() {
        // OH radical (doublet) at a typical bond length. The restricted
        // solution is above UHF by a few millihartree, not more.
        let mut m = Molecule::new("OH");
        m.add_atom::<8>(Vec3::ZERO);
        m.add_atom::<1>(Vec3::new(0.0, 0.0, 0.97));
        let r = rohf(&m, 2).unwrap();
        let u = uhf(&m, 2).unwrap();
        assert!(r.energy >= u.energy - 1e-8, "{} < {}", r.energy, u.energy);
        assert!(r.energy - u.energy < 0.02, "{} vs {}", r.energy, u.energy);
    }

    #[test]
    fn effective_fock_blocks_follow_roothaan_rules() {
        // Two orbitals closed, one open, one virtual (nb = 4) with distinct
        // Fα and Fβ so each block's source is visible.
        let nb = 4;
        let fa: Vec<f64> = (0..16).map(|i| 1.0 + i as f64).collect();
        let fb: Vec<f64> = (0..16).map(|i| 100.0 + i as f64).collect();
        let f = effective_fock(&fa, &fb, nb, 2, 3);
        let avg = |p: usize, q: usize| 0.5 * (fa[p * nb + q] + fb[p * nb + q]);
        assert_eq!(f[1], avg(0, 1)); // closed–closed
        assert_eq!(f[2], fb[2]); // closed–open
        assert_eq!(f[2 * nb + 3], fa[2 * nb + 3]); // open–virtual
        assert_eq!(f[3 * nb + 3], avg(3, 3)); // virtual–virtual
        assert_eq!(f[2 * nb + 2], avg(2, 2)); // open–open
    }
}
