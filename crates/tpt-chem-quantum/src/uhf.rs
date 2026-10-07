//! Unrestricted Hartree–Fock for open-shell systems (radicals, ions,
//! triplets). Separate alpha and beta orbitals:
//! `F^s = H + J[Da + Db] - K[D^s]`, `E = 1/2 sum (D.H + Da.Fa + Db.Fb) + E_nuc`.

use tpt_chem_core::molecule::Molecule;

use crate::basis::BasisSet;
use crate::hf::{scf_setup, HfError, PackedQuartet, ScfSetup};
use crate::linalg::{jacobi_eigh, mat_mul, transpose};

/// Result of a converged UHF calculation.
#[derive(Clone, Debug)]
pub struct UhfResult {
    /// Total energy (Hartree).
    pub energy: f64,
    /// Nuclear repulsion energy (Hartree).
    pub nuclear_energy: f64,
    /// Alpha orbital energies, ascending.
    pub orbital_energies_alpha: Vec<f64>,
    /// Beta orbital energies, ascending.
    pub orbital_energies_beta: Vec<f64>,
    /// Alpha MO coefficients, row-major (columns = orbitals).
    pub coefficients_alpha: Vec<f64>,
    /// Beta MO coefficients, row-major (columns = orbitals).
    pub coefficients_beta: Vec<f64>,
    /// Alpha density matrix (no factor 2).
    pub density_alpha: Vec<f64>,
    /// Beta density matrix (no factor 2).
    pub density_beta: Vec<f64>,
    /// Expectation value of S^2 (exact value is S(S+1); the excess is spin
    /// contamination).
    pub s_squared: f64,
    /// Number of alpha electrons.
    pub n_alpha: usize,
    /// Number of beta electrons.
    pub n_beta: usize,
    /// Number of basis functions.
    pub n_basis: usize,
}

fn spin_density(coeffs: &[f64], n_occ: usize, nb: usize) -> Vec<f64> {
    let mut d = vec![0.0; nb * nb];
    for i in 0..nb {
        for j in 0..nb {
            d[i * nb + j] = (0..n_occ)
                .map(|o| coeffs[i * nb + o] * coeffs[j * nb + o])
                .sum();
        }
    }
    d
}

/// Accumulate `J[Dt]`, `K[Da]` and `K[Db]` from the packed quartets.
fn accumulate_jk(
    quartets: &[PackedQuartet],
    offsets: &[usize],
    nb: usize,
    dt: &[f64],
    da: &[f64],
    db: &[f64],
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut j = vec![0.0; nb * nb];
    let mut ka = vec![0.0; nb * nb];
    let mut kb = vec![0.0; nb * nb];
    for q in quartets {
        let [a, b, c, d] = q.shells;
        let [na, nb2, nc, nd] = q.dims;
        let (oi, oj, ok, ol) = (offsets[a], offsets[b], offsets[c], offsets[d]);
        let mut idx = 0;
        for ia in 0..na {
            for ib in 0..nb2 {
                for ic in 0..nc {
                    for id in 0..nd {
                        let v = q.v[idx];
                        idx += 1;
                        let (i, jj, k, l) = (oi + ia, oj + ib, ok + ic, ol + id);
                        j[i * nb + jj] += dt[k * nb + l] * v;
                        ka[i * nb + k] += da[jj * nb + l] * v;
                        kb[i * nb + k] += db[jj * nb + l] * v;
                    }
                }
            }
        }
    }
    (j, ka, kb)
}

fn diag(fock: &[f64], x: &[f64], nb: usize) -> Result<(Vec<f64>, Vec<f64>), HfError> {
    let ft = transpose(fock, nb);
    let fp = mat_mul(&mat_mul(x, &ft, nb), x, nb);
    let (eps, cp) = jacobi_eigh(&fp, nb, 1e-11, 64).map_err(|_| HfError::NumericalBlowup)?;
    Ok((eps, mat_mul(x, &cp, nb)))
}

/// UHF energy with the given spin multiplicity `2S + 1` (STO-3G basis).
///
/// # Errors
/// [`HfError::UnsupportedElement`] if the multiplicity is incompatible with
/// the electron count (zero, wrong parity, or more unpaired electrons than
/// the basis can hold); otherwise see [`crate::hf::rhf_energy`].
pub fn uhf(mol: &Molecule, multiplicity: u32) -> Result<UhfResult, HfError> {
    uhf_with_basis(mol, multiplicity, BasisSet::Sto3g)
}

/// [`uhf`] with an explicit basis set.
///
/// # Errors
/// See [`uhf`].
pub fn uhf_with_basis(
    mol: &Molecule,
    multiplicity: u32,
    basis: BasisSet,
) -> Result<UhfResult, HfError> {
    let ScfSetup {
        nb,
        offsets,
        s,
        h_core,
        quartets,
        x_mat,
        e_nuc,
        n_elec,
    } = scf_setup(mol, basis)?;
    if multiplicity == 0 {
        return Err(HfError::UnsupportedElement);
    }
    let unpaired = (multiplicity - 1) as usize;
    if unpaired > n_elec || (n_elec - unpaired) % 2 != 0 {
        return Err(HfError::UnsupportedElement);
    }
    let n_beta = (n_elec - unpaired) / 2;
    let n_alpha = n_beta + unpaired;
    if n_alpha > nb {
        return Err(HfError::UnsupportedElement);
    }

    let zeros = vec![0.0; nb * nb];
    let (mut da, mut db) = (zeros.clone(), zeros);
    let mut e_old = 0.0;
    let mut last = None;
    let mut converged = false;
    for it in 0..256 {
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
        let (ea, ca) = diag(&fa, &x_mat, nb)?;
        let (eb, cb) = diag(&fb, &x_mat, nb)?;
        let da_new = spin_density(&ca, n_alpha, nb);
        let db_new = spin_density(&cb, n_beta, nb);
        let dd = da_new
            .iter()
            .zip(&da)
            .chain(db_new.iter().zip(&db))
            .map(|(n, o)| (n - o).abs())
            .fold(0.0f64, f64::max);
        last = Some((ea, ca, eb, cb, energy));
        if it > 0 && (energy - e_old).abs() < 1e-10 && dd < 1e-7 {
            converged = true;
            break;
        }
        e_old = energy;
        let damp = if it < 10 { 0.5 } else { 0.8 };
        for i in 0..nb * nb {
            da[i] = damp * da_new[i] + (1.0 - damp) * da[i];
            db[i] = damp * db_new[i] + (1.0 - damp) * db[i];
        }
    }
    if !converged {
        return Err(HfError::NotConverged);
    }
    let (ea, ca, eb, cb, energy) = last.ok_or(HfError::NotConverged)?;
    let da = spin_density(&ca, n_alpha, nb);
    let db = spin_density(&cb, n_beta, nb);

    // <S^2> = Sz(Sz+1) + n_beta - sum_ij |<psi_a,i | psi_b,j>|^2.
    let sz = 0.5 * (n_alpha as f64 - n_beta as f64);
    let sc = mat_mul(&s, &cb, nb);
    let mut overlap_sum = 0.0;
    for i in 0..n_alpha {
        for j in 0..n_beta {
            let o: f64 = (0..nb).map(|m| ca[m * nb + i] * sc[m * nb + j]).sum();
            overlap_sum += o * o;
        }
    }
    let s_squared = sz * (sz + 1.0) + n_beta as f64 - overlap_sum;

    Ok(UhfResult {
        energy,
        nuclear_energy: e_nuc,
        orbital_energies_alpha: ea,
        orbital_energies_beta: eb,
        coefficients_alpha: ca,
        coefficients_beta: cb,
        density_alpha: da,
        density_beta: db,
        s_squared,
        n_alpha,
        n_beta,
        n_basis: nb,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hf::rhf_energy;
    use tpt_chem_core::vec3::Vec3;

    fn h2(r: f64) -> Molecule {
        let mut m = Molecule::new("H2");
        m.add_atom::<1>(Vec3::new(0.0, 0.0, -r / 2.0));
        m.add_atom::<1>(Vec3::new(0.0, 0.0, r / 2.0));
        m
    }

    #[test]
    fn hydrogen_atom_doublet() {
        let mut h = Molecule::new("H");
        h.add_atom::<1>(Vec3::new(0.0, 0.0, 0.0));
        let r = uhf(&h, 2).unwrap();
        // HF/STO-3G hydrogen atom: -0.46658185 Eh.
        assert!((r.energy - (-0.466_581_85)).abs() < 1e-6, "{}", r.energy);
        assert!((r.s_squared - 0.75).abs() < 1e-9);
    }

    #[test]
    fn uhf_singlet_matches_rhf() {
        let m = h2(0.74);
        let u = uhf(&m, 1).unwrap();
        assert!((u.energy - rhf_energy(&m).unwrap()).abs() < 1e-8);
        assert!(u.s_squared.abs() < 1e-8);
    }

    #[test]
    fn h2_triplet_is_above_singlet() {
        let m = h2(0.74);
        let t = uhf(&m, 3).unwrap();
        assert!(t.energy > uhf(&m, 1).unwrap().energy + 0.3);
        assert!((t.s_squared - 2.0).abs() < 1e-9);
    }

    #[test]
    fn bad_multiplicity_rejected() {
        assert!(uhf(&h2(0.74), 2).is_err());
        assert!(uhf(&h2(0.74), 0).is_err());
        assert!(uhf(&h2(0.74), 5).is_err());
    }

    #[test]
    fn helium_cation_doublet() {
        let mut he = Molecule::new("He+");
        he.add_atom::<2>(Vec3::new(0.0, 0.0, 0.0));
        he.set_formal_charge(1);
        let r = uhf(&he, 2).unwrap();
        // One electron in a 1s STO-3G orbital on Z = 2: close to -1.99.
        assert!(r.energy < -1.9 && r.energy > -2.1, "{}", r.energy);
    }
}
