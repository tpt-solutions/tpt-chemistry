//! Pitzer ion-interaction model (Pitzer 1973, 1991 multicomponent form)
//! for electrolyte mixtures.
//!
//! Supported terms: Debye–Hückel long-range term, cation–anion `β⁰, β¹, β²,
//! Cᵠ`, like-charge `θ` and triplet `ψ`. Unsymmetric electrostatic mixing
//! terms (`ᴱθ`) and neutral-solute `λ` terms are not implemented; neutral
//! solutes receive `γ = 1`. Parameters not supplied default to zero, so an
//! empty parameter set reduces to the Debye–Hückel term alone.
//!
//! Parameter temperature dependence is linear about 298.15 K, `p(T) = p₂₉₈ +
//! (dp/dT)(T − 298.15)`.

use super::{check_state, ActivityModel};
use crate::constants::T_REF;
use crate::solvent::debye_huckel_a_phi;
use crate::species::{ionic_strength, Species};
use crate::Result;

/// Pitzer cation–anion parameters.
#[derive(Debug, Clone, PartialEq)]
pub struct PitzerPair {
    /// Cation name (must match [`Species::name`]).
    pub cation: String,
    /// Anion name.
    pub anion: String,
    /// `β⁰` at 298.15 K.
    pub beta0: f64,
    /// `β¹` at 298.15 K.
    pub beta1: f64,
    /// `β²` at 298.15 K (used for 2:2 and higher-valence-type pairs).
    pub beta2: f64,
    /// `Cᵠ` at 298.15 K.
    pub c_phi: f64,
    /// `α₁` (2.0 for most pairs, 1.4 for 2:2).
    pub alpha1: f64,
    /// `α₂` (12.0 where `β²` is used).
    pub alpha2: f64,
    /// `dβ⁰/dT`, K⁻¹.
    pub d_beta0: f64,
    /// `dβ¹/dT`, K⁻¹.
    pub d_beta1: f64,
    /// `dβ²/dT`, K⁻¹.
    pub d_beta2: f64,
    /// `dCᵠ/dT`, K⁻¹.
    pub d_c_phi: f64,
}

impl PitzerPair {
    /// Creates a pair with `α₁ = 2`, `α₂ = 12` and no temperature dependence.
    pub fn new(cation: &str, anion: &str, beta0: f64, beta1: f64, beta2: f64, c_phi: f64) -> Self {
        Self {
            cation: cation.to_string(),
            anion: anion.to_string(),
            beta0,
            beta1,
            beta2,
            c_phi,
            alpha1: 2.0,
            alpha2: 12.0,
            d_beta0: 0.0,
            d_beta1: 0.0,
            d_beta2: 0.0,
            d_c_phi: 0.0,
        }
    }
}

/// Like-charge mixing parameter `θ` between two ions of the same sign.
#[derive(Debug, Clone, PartialEq)]
pub struct PitzerTheta {
    /// First ion name.
    pub a: String,
    /// Second ion name.
    pub b: String,
    /// `θ` at 298.15 K.
    pub theta: f64,
    /// `dθ/dT`, K⁻¹.
    pub d_theta: f64,
}

/// Triplet parameter `ψ` for two like-charged ions and one counter-ion.
#[derive(Debug, Clone, PartialEq)]
pub struct PitzerPsi {
    /// The two like-charged ions.
    pub like: (String, String),
    /// The counter-ion of opposite sign.
    pub other: String,
    /// `ψ` (temperature independent).
    pub psi: f64,
}

/// Pitzer ion-interaction model.
#[derive(Debug, Clone, PartialEq)]
pub struct Pitzer {
    pairs: Vec<PitzerPair>,
    thetas: Vec<PitzerTheta>,
    psis: Vec<PitzerPsi>,
    /// Universal parameter `b` (1.2 (kg/mol)^½).
    pub b: f64,
}

impl Default for Pitzer {
    fn default() -> Self {
        Self::new()
    }
}

/// `g(x) = 2[1 − (1+x)e^{−x}]/x²`.
fn g(x: f64) -> f64 {
    if x < 1e-6 {
        1.0 - x + 0.75 * x * x
    } else {
        2.0 * (1.0 - (1.0 + x) * (-x).exp()) / (x * x)
    }
}

/// `g'(x) = −2[1 − (1 + x + x²/2)e^{−x}]/x²`.
fn gp(x: f64) -> f64 {
    if x < 1e-6 {
        -0.5 * x
    } else {
        -2.0 * (1.0 - (1.0 + x + 0.5 * x * x) * (-x).exp()) / (x * x)
    }
}

/// Parameter values for one cation–anion pair at the working temperature.
#[derive(Clone, Copy, Default)]
struct PairAt {
    b0: f64,
    b1: f64,
    b2: f64,
    cphi: f64,
    a1: f64,
    a2: f64,
}

impl PairAt {
    /// Second virial coefficient `B_ca` (activity form).
    fn b(&self, si: f64) -> f64 {
        self.b0 + self.b1 * g(self.a1 * si) + self.b2 * g(self.a2 * si)
    }
    /// `B'_ca = dB/dI`.
    fn b_prime(&self, si: f64, i: f64) -> f64 {
        (self.b1 * gp(self.a1 * si) + self.b2 * gp(self.a2 * si)) / i
    }
    /// `B^φ_ca` (osmotic form).
    fn b_phi(&self, si: f64) -> f64 {
        self.b0 + self.b1 * (-self.a1 * si).exp() + self.b2 * (-self.a2 * si).exp()
    }
}

struct Working<'a> {
    cats: Vec<usize>,
    ans: Vec<usize>,
    sp: &'a [Species],
    pairs: Vec<Vec<PairAt>>, // [cat][an]
    cc: Vec<Vec<f64>>,       // θ among cations
    aa: Vec<Vec<f64>>,       // θ among anions
    i: f64,
    si: f64,
    z_sum: f64,
}

impl Pitzer {
    /// Creates an empty parameter set (Debye–Hückel term only).
    pub fn new() -> Self {
        Self {
            pairs: Vec::new(),
            thetas: Vec::new(),
            psis: Vec::new(),
            b: 1.2,
        }
    }

    /// Adds (or replaces) cation–anion parameters.
    pub fn add_pair(&mut self, pair: PitzerPair) {
        self.pairs
            .retain(|p| !(p.cation == pair.cation && p.anion == pair.anion));
        self.pairs.push(pair);
    }

    /// Adds like-charge `θ` parameters.
    pub fn add_theta(&mut self, theta: PitzerTheta) {
        self.thetas.push(theta);
    }

    /// Adds triplet `ψ` parameters.
    pub fn add_psi(&mut self, psi: PitzerPsi) {
        self.psis.push(psi);
    }

    fn pair_at(&self, c: &str, a: &str, t: f64) -> PairAt {
        match self.pairs.iter().find(|p| p.cation == c && p.anion == a) {
            Some(p) => {
                let dt = t - T_REF;
                PairAt {
                    b0: p.beta0 + p.d_beta0 * dt,
                    b1: p.beta1 + p.d_beta1 * dt,
                    b2: p.beta2 + p.d_beta2 * dt,
                    cphi: p.c_phi + p.d_c_phi * dt,
                    a1: p.alpha1,
                    a2: p.alpha2,
                }
            }
            None => PairAt {
                a1: 2.0,
                a2: 12.0,
                ..PairAt::default()
            },
        }
    }

    fn theta_at(&self, x: &str, y: &str, t: f64) -> f64 {
        self.thetas
            .iter()
            .find(|p| (p.a == x && p.b == y) || (p.a == y && p.b == x))
            .map_or(0.0, |p| p.theta + p.d_theta * (t - T_REF))
    }

    fn psi_at(&self, x: &str, y: &str, other: &str) -> f64 {
        self.psis
            .iter()
            .find(|p| {
                p.other == other
                    && ((p.like.0 == x && p.like.1 == y) || (p.like.0 == y && p.like.1 == x))
            })
            .map_or(0.0, |p| p.psi)
    }

    fn prepare<'a>(&self, sp: &'a [Species], t: f64) -> Working<'a> {
        let cats: Vec<usize> = (0..sp.len()).filter(|&k| sp[k].charge > 0).collect();
        let ans: Vec<usize> = (0..sp.len()).filter(|&k| sp[k].charge < 0).collect();
        let pairs = cats
            .iter()
            .map(|&c| {
                ans.iter()
                    .map(|&a| self.pair_at(&sp[c].name, &sp[a].name, t))
                    .collect()
            })
            .collect();
        let cc = cats
            .iter()
            .map(|&c| {
                cats.iter()
                    .map(|&d| self.theta_at(&sp[c].name, &sp[d].name, t))
                    .collect()
            })
            .collect();
        let aa = ans
            .iter()
            .map(|&a| {
                ans.iter()
                    .map(|&b| self.theta_at(&sp[a].name, &sp[b].name, t))
                    .collect()
            })
            .collect();
        let i = ionic_strength(sp);
        let z_sum = sp
            .iter()
            .map(|s| s.molality * f64::from(s.charge.abs()))
            .sum();
        Working {
            cats,
            ans,
            sp,
            pairs,
            cc,
            aa,
            i,
            si: i.sqrt(),
            z_sum,
        }
    }

    /// `C_ca = Cᵠ / (2 √|z_c z_a|)`.
    fn c_ca(w: &Working<'_>, ci: usize, ai: usize) -> f64 {
        let zc = f64::from(w.sp[w.cats[ci]].charge);
        let za = f64::from(w.sp[w.ans[ai]].charge);
        w.pairs[ci][ai].cphi / (2.0 * (zc * za.abs()).sqrt())
    }

    fn f_gamma(&self, w: &Working<'_>, a_phi: f64) -> f64 {
        let b = self.b;
        let mut f = -a_phi * (w.si / (1.0 + b * w.si) + 2.0 / b * (1.0 + b * w.si).ln());
        for (ci, &c) in w.cats.iter().enumerate() {
            for (ai, &a) in w.ans.iter().enumerate() {
                f += w.sp[c].molality * w.sp[a].molality * w.pairs[ci][ai].b_prime(w.si, w.i);
            }
        }
        f
    }
}

impl ActivityModel for Pitzer {
    fn ln_gamma(&self, species: &[Species], t: f64) -> Result<Vec<f64>> {
        check_state(species, t)?;
        let mut out = vec![0.0; species.len()];
        let w = self.prepare(species, t);
        if w.i < 1e-14 || w.cats.is_empty() || w.ans.is_empty() {
            return Ok(out);
        }
        let sp = species;
        let f = self.f_gamma(&w, debye_huckel_a_phi(t));
        // Σ_c Σ_a m_c m_a C_ca
        let mut cc_sum = 0.0;
        for ci in 0..w.cats.len() {
            for ai in 0..w.ans.len() {
                cc_sum += sp[w.cats[ci]].molality * sp[w.ans[ai]].molality * Self::c_ca(&w, ci, ai);
            }
        }
        for (mi, &m) in w.cats.iter().enumerate() {
            let zm = f64::from(sp[m].charge);
            let mut s = zm * zm * f + zm.abs() * cc_sum;
            for (ai, &a) in w.ans.iter().enumerate() {
                s += sp[a].molality
                    * (2.0 * w.pairs[mi][ai].b(w.si) + w.z_sum * Self::c_ca(&w, mi, ai));
            }
            for (ci, &c) in w.cats.iter().enumerate() {
                if ci == mi {
                    continue;
                }
                let mut inner = 2.0 * w.cc[mi][ci];
                for &a in &w.ans {
                    inner += sp[a].molality * self.psi_at(&sp[m].name, &sp[c].name, &sp[a].name);
                }
                s += sp[c].molality * inner;
            }
            for (ai, &a) in w.ans.iter().enumerate() {
                for &a2 in &w.ans[ai + 1..] {
                    s += sp[a].molality
                        * sp[a2].molality
                        * self.psi_at(&sp[a].name, &sp[a2].name, &sp[m].name);
                }
            }
            out[m] = s;
        }
        for (xi, &x) in w.ans.iter().enumerate() {
            let zx = f64::from(sp[x].charge);
            let mut s = zx * zx * f + zx.abs() * cc_sum;
            for (ci, &c) in w.cats.iter().enumerate() {
                s += sp[c].molality
                    * (2.0 * w.pairs[ci][xi].b(w.si) + w.z_sum * Self::c_ca(&w, ci, xi));
            }
            for (ai, &a) in w.ans.iter().enumerate() {
                if ai == xi {
                    continue;
                }
                let mut inner = 2.0 * w.aa[xi][ai];
                for &c in &w.cats {
                    inner += sp[c].molality * self.psi_at(&sp[x].name, &sp[a].name, &sp[c].name);
                }
                s += sp[a].molality * inner;
            }
            for (ci, &c) in w.cats.iter().enumerate() {
                for &c2 in &w.cats[ci + 1..] {
                    s += sp[c].molality
                        * sp[c2].molality
                        * self.psi_at(&sp[c].name, &sp[c2].name, &sp[x].name);
                }
            }
            out[x] = s;
        }
        Ok(out)
    }

    fn osmotic_coefficient(&self, species: &[Species], t: f64) -> Result<f64> {
        check_state(species, t)?;
        let w = self.prepare(species, t);
        let sp = species;
        let total_m: f64 = w.cats.iter().chain(&w.ans).map(|&k| sp[k].molality).sum();
        if w.i < 1e-14 || total_m <= 0.0 {
            return Ok(1.0);
        }
        let b = self.b;
        let mut s = -debye_huckel_a_phi(t) * w.i * w.si / (1.0 + b * w.si);
        for (ci, &c) in w.cats.iter().enumerate() {
            for (ai, &a) in w.ans.iter().enumerate() {
                let p = &w.pairs[ci][ai];
                s += sp[c].molality
                    * sp[a].molality
                    * (p.b_phi(w.si) + w.z_sum * Self::c_ca(&w, ci, ai));
            }
        }
        for (ci, &c) in w.cats.iter().enumerate() {
            for (cj, &c2) in w.cats.iter().enumerate().skip(ci + 1) {
                let mut inner = w.cc[ci][cj];
                for &a in &w.ans {
                    inner += sp[a].molality * self.psi_at(&sp[c].name, &sp[c2].name, &sp[a].name);
                }
                s += sp[c].molality * sp[c2].molality * inner;
            }
        }
        for (ai, &a) in w.ans.iter().enumerate() {
            for (aj, &a2) in w.ans.iter().enumerate().skip(ai + 1) {
                let mut inner = w.aa[ai][aj];
                for &c in &w.cats {
                    inner += sp[c].molality * self.psi_at(&sp[a].name, &sp[a2].name, &sp[c].name);
                }
                s += sp[a].molality * sp[a2].molality * inner;
            }
        }
        Ok(1.0 + 2.0 * s / total_m)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn salt(c: &str, zc: i32, a: &str, za: i32, m: f64) -> Vec<Species> {
        // electroneutral 1 formula unit per `m`
        vec![
            Species::new(c, zc, m * f64::from(za.abs())),
            Species::new(a, za, m * f64::from(zc)),
        ]
    }

    fn mean_gamma(model: &Pitzer, sol: &[Species], nu: (f64, f64)) -> f64 {
        let g = model.ln_gamma(sol, 298.15).unwrap();
        ((nu.0 * g[0] + nu.1 * g[1]) / (nu.0 + nu.1)).exp()
    }

    fn nacl_model() -> Pitzer {
        let mut p = Pitzer::new();
        p.add_pair(PitzerPair::new("Na+", "Cl-", 0.0765, 0.2664, 0.0, 0.00127));
        p
    }

    #[test]
    fn nacl_one_molal() {
        let model = nacl_model();
        let sol = salt("Na+", 1, "Cl-", -1, 1.0);
        assert!((mean_gamma(&model, &sol, (1.0, 1.0)) - 0.657).abs() < 0.005);
        let phi = model.osmotic_coefficient(&sol, 298.15).unwrap();
        assert!((phi - 0.937).abs() < 0.004, "phi = {phi}");
    }

    #[test]
    fn kcl_one_molal() {
        let mut model = Pitzer::new();
        model.add_pair(PitzerPair::new("K+", "Cl-", 0.04835, 0.2122, 0.0, -0.00084));
        let sol = salt("K+", 1, "Cl-", -1, 1.0);
        assert!((mean_gamma(&model, &sol, (1.0, 1.0)) - 0.604).abs() < 0.006);
    }

    #[test]
    fn cacl2_one_molal() {
        let mut model = Pitzer::new();
        model.add_pair(PitzerPair::new("Ca2+", "Cl-", 0.3159, 1.614, 0.0, -0.00034));
        let sol = salt("Ca2+", 2, "Cl-", -1, 1.0);
        let g = mean_gamma(&model, &sol, (1.0, 2.0));
        assert!((g - 0.50).abs() < 0.03, "gamma = {g}");
    }

    #[test]
    fn gamma_tends_to_one_when_dilute() {
        let g = nacl_model()
            .ln_gamma(&salt("Na+", 1, "Cl-", -1, 1e-8), 298.15)
            .unwrap();
        assert!(g[0].abs() < 1e-3 && g[1].abs() < 1e-3);
    }

    #[test]
    fn gibbs_duhem_consistency_between_gamma_and_phi() {
        // ln γ± = (φ − 1) + ∫₀ᵐ (φ − 1)/m dm
        let model = nacl_model();
        let m_end = 1.5_f64;
        let phi = |m: f64| {
            model
                .osmotic_coefficient(&salt("Na+", 1, "Cl-", -1, m), 298.15)
                .unwrap()
        };
        // Substitute m = u² so the integrand (φ−1)/m dm = 2(φ−1)/u du is smooth.
        let n = 4000;
        let u_end = m_end.sqrt();
        let h = u_end / f64::from(n);
        let mut integral = 0.0;
        for k in 0..n {
            let u = (f64::from(k) + 0.5) * h;
            integral += 2.0 * (phi(u * u) - 1.0) / u * h;
        }
        let lhs = {
            let g = model
                .ln_gamma(&salt("Na+", 1, "Cl-", -1, m_end), 298.15)
                .unwrap();
            0.5 * (g[0] + g[1])
        };
        let rhs = (phi(m_end) - 1.0) + integral;
        assert!((lhs - rhs).abs() < 2e-3, "{lhs} vs {rhs}");
    }

    #[test]
    fn mixing_theta_changes_mixture_gamma() {
        let mut with = nacl_model();
        with.add_pair(PitzerPair::new("K+", "Cl-", 0.04835, 0.2122, 0.0, -0.00084));
        let mut without = with.clone();
        with.add_theta(PitzerTheta {
            a: "Na+".into(),
            b: "K+".into(),
            theta: -0.012,
            d_theta: 0.0,
        });
        let sol = [
            Species::new("Na+", 1, 0.5),
            Species::new("K+", 1, 0.5),
            Species::new("Cl-", -1, 1.0),
        ];
        let a = with.ln_gamma(&sol, 298.15).unwrap();
        let b = without.ln_gamma(&sol, 298.15).unwrap();
        assert!((a[0] - b[0]).abs() > 1e-4);
        without.add_theta(PitzerTheta {
            a: "Na+".into(),
            b: "K+".into(),
            theta: -0.012,
            d_theta: 0.0,
        });
        assert_eq!(without.ln_gamma(&sol, 298.15).unwrap(), a);
    }

    #[test]
    fn empty_parameter_set_gives_debye_huckel_term_only() {
        let model = Pitzer::new();
        let g = model
            .ln_gamma(&salt("Na+", 1, "Cl-", -1, 0.001), 298.15)
            .unwrap();
        let limiting = -crate::solvent::debye_huckel_a_gamma(298.15) * 0.001_f64.sqrt();
        assert!((g[0] - limiting).abs() / limiting.abs() < 0.08);
    }
}

/// Validation against measured mean activity coefficients (25 °C).
///
/// The data are mean activity coefficients from the CRC Handbook of Chemistry
/// and Physics (97th ed.), entered as numbers. They were checked against a
/// compiled copy of the table, so see `docs/third-party.md` for provenance.
/// Parameters are the binary Pitzer set used by the tests above (NaCl, KCl,
/// CaCl₂), and PHREEQC's `pitzer.dat` (25 °C column) for MgCl₂.
#[cfg(test)]
mod crc_validation {
    use super::*;

    fn gamma_mean(model: &Pitzer, cation: &str, z: i32, anion: &str, za: i32, m: f64) -> f64 {
        let sol = vec![
            Species::new(cation, z, m * f64::from(za.abs())),
            Species::new(anion, za, m * f64::from(z)),
        ];
        let g = model.ln_gamma(&sol, 298.15).unwrap();
        // ν₊ and ν₋ for a 1:1 salt (1,1) or a 2:1 salt (2,1) with cation charge 2.
        let (nu_c, nu_a) = if z.abs() == za.abs() {
            (1.0, 1.0)
        } else if z == 2 {
            (1.0, 2.0)
        } else {
            (2.0, 1.0)
        };
        ((nu_c * g[0] + nu_a * g[1]) / (nu_c + nu_a)).exp()
    }

    /// Largest |predicted − measured| over the points with `m ≤ m_max`.
    fn worst_deviation(
        model: &Pitzer,
        data: &[(f64, f64)],
        salt: (&str, i32, &str, i32),
        m_max: f64,
    ) -> f64 {
        data.iter()
            .filter(|(m, _)| *m <= m_max + 1e-12)
            .map(|&(m, measured)| {
                let predicted = gamma_mean(model, salt.0, salt.1, salt.2, salt.3, m);
                (predicted - measured).abs()
            })
            .fold(0.0, f64::max)
    }

    const NACL: [(f64, f64); 12] = [
        (0.1, 0.778),
        (0.2, 0.735),
        (0.3, 0.71),
        (0.4, 0.693),
        (0.5, 0.681),
        (0.6, 0.673),
        (0.7, 0.667),
        (0.8, 0.662),
        (0.9, 0.659),
        (1.0, 0.657),
        (2.0, 0.668),
        (5.0, 0.874),
    ];
    const KCL: [(f64, f64); 12] = [
        (0.1, 0.77),
        (0.2, 0.718),
        (0.3, 0.688),
        (0.4, 0.666),
        (0.5, 0.649),
        (0.6, 0.637),
        (0.7, 0.626),
        (0.8, 0.618),
        (0.9, 0.61),
        (1.0, 0.604),
        (2.0, 0.573),
        (5.0, 0.593),
    ];
    const CACL2: [(f64, f64); 11] = [
        (0.1, 0.518),
        (0.2, 0.472),
        (0.3, 0.455),
        (0.4, 0.448),
        (0.5, 0.448),
        (0.6, 0.453),
        (0.7, 0.46),
        (0.8, 0.47),
        (0.9, 0.484),
        (1.0, 0.5),
        (2.0, 0.784),
    ];
    const MGCL2: [(f64, f64); 10] = [
        (0.1, 0.529),
        (0.2, 0.489),
        (0.3, 0.477),
        (0.4, 0.475),
        (0.5, 0.481),
        (0.6, 0.491),
        (0.7, 0.506),
        (0.8, 0.522),
        (0.9, 0.544),
        (1.0, 0.57),
    ];

    #[test]
    fn nacl_matches_measured_series_to_5_molal() {
        let mut model = Pitzer::new();
        model.add_pair(PitzerPair::new("Na+", "Cl-", 0.0765, 0.2664, 0.0, 0.00127));
        let worst = worst_deviation(&model, &NACL, ("Na+", 1, "Cl-", -1), 5.0);
        assert!(worst < 0.005, "NaCl worst |Δγ±| = {worst:.4}");
    }

    #[test]
    fn kcl_matches_measured_series_to_5_molal() {
        let mut model = Pitzer::new();
        model.add_pair(PitzerPair::new("K+", "Cl-", 0.04835, 0.2122, 0.0, -0.00084));
        let worst = worst_deviation(&model, &KCL, ("K+", 1, "Cl-", -1), 5.0);
        assert!(worst < 0.005, "KCl worst |Δγ±| = {worst:.4}");
    }

    #[test]
    fn cacl2_matches_measured_series_to_1_molal() {
        // A binary-only model is not accurate above 1 mol/kg for CaCl₂ (the
        // 2 mol/kg point deviates by about 0.016), so the test stops at 1.
        let mut model = Pitzer::new();
        model.add_pair(PitzerPair::new("Ca2+", "Cl-", 0.3159, 1.614, 0.0, -0.00034));
        let worst = worst_deviation(&model, &CACL2, ("Ca2+", 2, "Cl-", -1), 1.0);
        assert!(worst < 0.003, "CaCl2 worst |Δγ±| = {worst:.4}");
    }

    #[test]
    fn mgcl2_matches_measured_series_to_1_molal() {
        let mut model = Pitzer::new();
        model.add_pair(PitzerPair::new("Mg+2", "Cl-", 0.351, 1.65, 0.0, 0.00651));
        let worst = worst_deviation(&model, &MGCL2, ("Mg+2", 2, "Cl-", -1), 1.0);
        assert!(worst < 0.008, "MgCl2 worst |Δγ±| = {worst:.4}");
    }
}
