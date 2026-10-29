//! STO-3G basis set data for hydrogen through neon.
//!
//! The standard STO-3G parameterization: every atomic orbital is fit by 3
//! contracted primitives; the contraction coefficients are universal for a
//! given shell class (inner s: [0.15432897, 0.53532814, 0.44463454];
//! valence s: [−0.09996723, 0.39951283, 0.70011547]; valence p:
//! [0.15591627, 0.60768372, 0.39195739]) — only the exponents vary per
//! element (Hehre, Stewart & Pople 1969).
//!
//! Values are the Basis Set Exchange "STO-3G" release numbers.

use tpt_chem_core::vec3::Vec3;

use crate::gaussian::Shell;

/// Contraction coefficients shared by all STO-3G shells of a class.
const C_INNER_S: [f64; 3] = [0.15432897, 0.53532814, 0.44463454];
const C_VALENCE_S: [f64; 3] = [-0.09996723, 0.39951283, 0.70011547];
const C_VALENCE_P: [f64; 3] = [0.15591627, 0.60768372, 0.39195739];

/// Inner-shell (1s for Li–Ne) exponents.
const INNER_S: [(u8, [f64; 3]); 8] = [
    (3, [16.1195750, 2.8698150, 0.7657880]),
    (4, [27.1045800, 4.8224560, 1.2895730]),
    (5, [45.6970220, 8.5182460, 2.2918220]),
    (6, [71.6168370, 13.0450960, 3.5305120]),
    (7, [99.1061690, 18.0527880, 4.8862190]),
    (8, [130.7093200, 23.8088610, 6.4436080]),
    (9, [166.6791300, 30.3608120, 8.2168200]),
    (10, [208.9629200, 38.0581570, 10.2055430]),
];

/// Valence (2s2p) SP-shell exponents for Li–Ne.
const VALENCE_SP: [(u8, [f64; 3]); 8] = [
    (3, [0.6362890, 0.1478600, 0.0488040]),
    (4, [1.0986560, 0.2554210, 0.0844180]),
    (5, [1.9420880, 0.4595040, 0.1506770]),
    (6, [2.9412490, 0.6834830, 0.2222899]),
    (7, [3.7804550, 0.8784970, 0.2857140]),
    (8, [5.0331510, 1.1695960, 0.3803890]),
    (9, [6.4648030, 1.5022810, 0.4885880]),
    (10, [8.2168200, 1.9181380, 0.6224140]),
];

/// STO-3G exponents for H (1s) and He (1s). He uses its own Basis Set
/// Exchange STO-3G exponent set (not the naive H-rescale), reproducing the
/// literature HF/STO-3G energy −2.8077840 Eₕ.
const H_HE_1S: [(u8, [f64; 3]); 2] = [
    (1, [3.42525091, 0.62391373, 0.16885540]),
    (2, [6.36242139, 1.15887060, 0.31364979]),
];

/// Errors from basis-set construction.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BasisError {
    /// The element is outside the parameterized set (H–Ne).
    UnsupportedElement(u8),
}

impl core::fmt::Display for BasisError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            BasisError::UnsupportedElement(z) => {
                write!(f, "STO-3G not parameterized for Z = {z} (supported: 1–10)")
            }
        }
    }
}

impl std::error::Error for BasisError {}

/// Build the minimal STO-3G basis for one atom of element `z` at `center`
/// (Bohr): H/He get a single 1s shell; Li–Ne get an inner 1s shell plus an
/// SP valence shell (emitted here as one s and one p shell sharing
/// exponents).
///
/// # Errors
/// [`BasisError::UnsupportedElement`] for Z > 10.
pub fn sto3g_shells(z: u8, center: Vec3) -> Result<Vec<Shell>, BasisError> {
    let mk = |l: u8, exps: [f64; 3], coefs: &[f64; 3]| {
        let mut sh = Shell::new(
            l,
            center,
            exps.iter().copied().zip(coefs.iter().copied()).collect(),
        );
        sh.normalize();
        sh
    };
    match z {
        1 | 2 => {
            let (_, exps) = H_HE_1S[(z - 1) as usize];
            Ok(vec![mk(0, exps, &C_INNER_S)])
        }
        3..=10 => {
            let idx = (z - 3) as usize;
            let (_, inner) = INNER_S[idx];
            let (_, valence) = VALENCE_SP[idx];
            Ok(vec![
                mk(0, inner, &C_INNER_S),
                mk(0, valence, &C_VALENCE_S),
                mk(1, valence, &C_VALENCE_P),
            ])
        }
        other => Err(BasisError::UnsupportedElement(other)),
    }
}

/// Build the STO-3G basis for a whole geometry: `atoms` is a list of
/// `(atomic number, position in Bohr)`.
///
/// The returned shells are in the canonical order: all shells of atom 0,
/// then atom 1, …
///
/// # Errors
/// [`BasisError::UnsupportedElement`] if any atom is outside H–Ne.
pub fn sto3g_basis(atoms: &[(u8, Vec3)]) -> Result<Vec<Shell>, BasisError> {
    let mut out = Vec::new();
    for &(z, r) in atoms {
        out.extend(sto3g_shells(z, r)?);
    }
    Ok(out)
}

/// 6-31G shells for element `z` (H–Ar) at `center` (Bohr).
///
/// # Errors
/// [`BasisError::UnsupportedElement`] for Z > 18.
pub fn pople_631g_shells(z: u8, center: Vec3) -> Result<Vec<Shell>, BasisError> {
    match z {
        1..=10 => Ok(crate::basis_631g::POPLE_631G[(z - 1) as usize]
            .iter()
            .map(|&(l, prims)| {
                let mut sh = Shell::new(l, center, prims.to_vec());
                sh.normalize();
                sh
            })
            .collect()),
        11..=18 => {
            let parsed = parse_element_shells(POPLE_631G_NA_AR_DATA, z)
                .ok_or(BasisError::UnsupportedElement(z))?;
            Ok(build_shells(parsed, center))
        }
        _ => Err(BasisError::UnsupportedElement(z)),
    }
}

/// 6-31G data for Na–Ar, generated from Basis Set Exchange by
/// `tools/gen_631g_na_ar.py`.
const POPLE_631G_NA_AR_DATA: &str = include_str!("../data/631g-na-ar.txt");

/// cc-pVDZ data for H–Ar, generated from Basis Set Exchange by
/// `tools/gen_cc_pvdz.py`.
const CC_PVDZ_DATA: &str = include_str!("../data/cc-pvdz.txt");

/// A shell read from the cc-pVDZ data: `(l, spherical, primitives)`.
type ParsedShell = (u8, bool, Vec<(f64, f64)>);

/// cc-pVDZ shells for element `z` (H–Ar) at `center` (Bohr). d shells are
/// spherical (five functions), the rest are Cartesian.
///
/// # Errors
/// [`BasisError::UnsupportedElement`] outside H–Ar.
pub fn cc_pvdz_shells(z: u8, center: Vec3) -> Result<Vec<Shell>, BasisError> {
    let parsed = if (1..=18).contains(&z) {
        parse_element_shells(CC_PVDZ_DATA, z)
    } else {
        None
    };
    parsed
        .map(|parsed| build_shells(parsed, center))
        .ok_or(BasisError::UnsupportedElement(z))
}

/// Parses the shells of element `z` from a data file in the
/// `@ <Z>` / `shell <l> <spherical>` / `<exponent> <coefficient>` format.
/// Returns `None` if the element is absent or the data is malformed.
fn parse_element_shells(data: &str, z: u8) -> Option<Vec<ParsedShell>> {
    let mut parsed: Vec<ParsedShell> = Vec::new();
    let mut inside = false;
    for line in data.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if let Some(element) = line.strip_prefix("@ ") {
            if inside {
                break;
            }
            inside = element.trim() == z.to_string();
            continue;
        }
        if !inside {
            continue;
        }
        if let Some(header) = line.strip_prefix("shell ") {
            let mut fields = header.split_whitespace();
            let l = fields.next().and_then(|x| x.parse::<u8>().ok())?;
            let sph = fields.next().and_then(|x| x.parse::<u8>().ok())?;
            parsed.push((l, sph == 1, Vec::new()));
        } else {
            let mut fields = line.split_whitespace();
            let alpha = fields.next().and_then(|x| x.parse::<f64>().ok())?;
            let coef = fields.next().and_then(|x| x.parse::<f64>().ok())?;
            parsed.last_mut()?.2.push((alpha, coef));
        }
    }
    if parsed.is_empty() {
        None
    } else {
        Some(parsed)
    }
}

/// Normalized shells from parsed data, with spherical flags applied.
fn build_shells(parsed: Vec<ParsedShell>, center: Vec3) -> Vec<Shell> {
    parsed
        .into_iter()
        .map(|(l, sph, prims)| {
            let mut sh = Shell::new(l, center, prims);
            if sph {
                sh = sh.with_spherical();
            }
            sh.normalize();
            sh
        })
        .collect()
}

/// Available Gaussian basis sets.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum BasisSet {
    /// Minimal STO-3G (H–Ne).
    #[default]
    Sto3g,
    /// Split-valence 6-31G (H–Ne); Cartesian s and p functions only.
    Pople631g,
    /// Dunning cc-pVDZ (H–Ne), with spherical d functions. Energies only:
    /// gradients, cube and Molden output are not implemented for it yet.
    CcPvdz,
}

impl BasisSet {
    /// Shells of one element at `center` (Bohr).
    ///
    /// # Errors
    /// [`BasisError::UnsupportedElement`] outside H–Ne.
    pub fn shells(self, z: u8, center: Vec3) -> Result<Vec<Shell>, BasisError> {
        match self {
            BasisSet::Sto3g => sto3g_shells(z, center),
            BasisSet::Pople631g => pople_631g_shells(z, center),
            BasisSet::CcPvdz => cc_pvdz_shells(z, center),
        }
    }

    /// Whether any shell of this basis uses spherical harmonics.
    pub fn is_spherical(self) -> bool {
        matches!(self, BasisSet::CcPvdz)
    }

    /// Shells for a whole geometry (`(Z, position in Bohr)`), atom by atom.
    ///
    /// # Errors
    /// See [`BasisSet::shells`].
    pub fn build(self, atoms: &[(u8, Vec3)]) -> Result<Vec<Shell>, BasisError> {
        let mut out = Vec::new();
        for &(z, r) in atoms {
            out.extend(self.shells(z, r)?);
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hydrogen_single_shell() {
        let shells = sto3g_shells(1, Vec3::ZERO).unwrap();
        assert_eq!(shells.len(), 1);
        assert_eq!(shells[0].l, 0);
        assert_eq!(shells[0].primitives.len(), 3);
        // Normalized to unit self-overlap.
        let s_ov = crate::integrals::overlap_shell(&shells[0], &shells[0])[0][0];
        assert!((s_ov - 1.0).abs() < 1e-12, "S = {s_ov}");
    }

    #[test]
    fn carbon_three_shells() {
        let shells = sto3g_shells(6, Vec3::new(0.1, 0.2, 0.3)).unwrap();
        assert_eq!(shells.len(), 3);
        assert_eq!(shells[0].l, 0); // 1s
        assert_eq!(shells[1].l, 0); // 2s
        assert_eq!(shells[2].l, 1); // 2p
                                    // The 2s and 2p shells share exponents.
        let exps = |sh: &Shell| sh.primitives.iter().map(|p| p.0).collect::<Vec<_>>();
        assert_eq!(exps(&shells[1]), exps(&shells[2]));
    }

    #[test]
    fn heavy_element_rejected() {
        assert_eq!(
            sto3g_shells(11, Vec3::ZERO),
            Err(BasisError::UnsupportedElement(11))
        );
    }

    #[test]
    fn pople_631g_na_ar_function_counts() {
        // Na and Ar in BSE 6-31G: 1s, 2sp, 3sp, 4sp (13 functions each).
        let na = pople_631g_shells(11, Vec3::ZERO).unwrap();
        let ar = pople_631g_shells(18, Vec3::ZERO).unwrap();
        let count = |shells: &[Shell]| {
            shells
                .iter()
                .map(|s| crate::gaussian::n_cartesian(s.l))
                .sum::<usize>()
        };
        assert_eq!(count(&na), 13);
        assert_eq!(count(&ar), 13);
    }

    #[test]
    fn pople_631g_na_ar_shells_are_normalized() {
        use crate::integrals::overlap_shell;
        for z in 11..=18 {
            for sh in pople_631g_shells(z, Vec3::ZERO).unwrap() {
                let s = overlap_shell(&sh, &sh);
                for (i, row) in s.iter().enumerate() {
                    assert!(
                        (row[i] - 1.0).abs() < 1e-8,
                        "Z={z} l={} diag {}",
                        sh.l,
                        row[i]
                    );
                }
            }
        }
    }

    #[test]
    fn pople_631g_rejects_beyond_argon() {
        assert!(pople_631g_shells(19, Vec3::ZERO).is_err());
    }
}
