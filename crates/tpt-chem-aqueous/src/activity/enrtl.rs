//! Single-solvent electrolyte NRTL (eNRTL) in the Chen–Evans local
//! composition form with a Pitzer–Debye–Hückel long-range term.
//!
//! The excess Gibbs energy `nGᴱˣ/RT` is evaluated as a function of mole
//! numbers and activity coefficients are obtained as its derivatives
//! (central differences), so the model is thermodynamically consistent by
//! construction (Gibbs–Duhem holds identically).
//!
//! Standard simplifications are used: water is the only solvent, ions of the
//! same sign have no direct interaction (`G = 1, τ = 0`), and the binary
//! parameters are the ion-pair/solvent set `τ_ca,m` and `τ_m,ca`. Ion
//! activity coefficients are converted to the unsymmetric (infinite-dilution)
//! convention and the molality scale. **No parameter sets are shipped**;
//! `τ` must be regressed or supplied by the user. Neutral solutes are
//! ignored (`γ = 1`).
//!
//! The infinite-dilution reference for each ion is evaluated with all solutes
//! scaled to a vanishing concentration at fixed relative proportions.

use super::{check_state, ActivityModel};
use crate::constants::M_WATER;
use crate::solvent::debye_huckel_a_phi;
use crate::species::Species;
use crate::{AqueousError, Result};

/// Binary eNRTL parameters for a cation–anion pair in water.
#[derive(Debug, Clone, PartialEq)]
pub struct EnrtlPair {
    /// Cation name.
    pub cation: String,
    /// Anion name.
    pub anion: String,
    /// `τ_ca,m` (energy parameter of the ion pair in the solvent cell).
    pub tau_ca_m: f64,
    /// `τ_m,ca` (energy parameter of the solvent about the ion-pair cell).
    pub tau_m_ca: f64,
}

/// Single-solvent eNRTL model.
#[derive(Debug, Clone, PartialEq)]
pub struct Enrtl {
    /// Non-randomness factor α (0.2 is standard for water).
    pub alpha: f64,
    /// Closest-approach parameter ρ of the Pitzer–Debye–Hückel term (14.9).
    pub rho: f64,
    pairs: Vec<EnrtlPair>,
}

impl Default for Enrtl {
    fn default() -> Self {
        Self {
            alpha: 0.2,
            rho: 14.9,
            pairs: Vec::new(),
        }
    }
}

struct Layout {
    z: Vec<f64>,
    n: Vec<f64>,
    cats: Vec<usize>,
    ans: Vec<usize>,
    tau1: Vec<Vec<f64>>, // [cat][an]
    tau2: Vec<Vec<f64>>,
    idx: Vec<usize>, // species index for each ion slot
}

impl Enrtl {
    /// Creates a model with no pair parameters (all `τ = 0`; PDH term only).
    pub fn new() -> Self {
        Self::default()
    }

    /// Adds (or replaces) the parameters of a cation–anion pair.
    pub fn add_pair(&mut self, pair: EnrtlPair) {
        self.pairs
            .retain(|p| !(p.cation == pair.cation && p.anion == pair.anion));
        self.pairs.push(pair);
    }

    fn layout(&self, species: &[Species]) -> Layout {
        let idx: Vec<usize> = (0..species.len())
            .filter(|&k| species[k].charge != 0)
            .collect();
        let z: Vec<f64> = idx.iter().map(|&k| f64::from(species[k].charge)).collect();
        let n: Vec<f64> = idx.iter().map(|&k| species[k].molality).collect();
        let cats: Vec<usize> = (0..idx.len()).filter(|&s| z[s] > 0.0).collect();
        let ans: Vec<usize> = (0..idx.len()).filter(|&s| z[s] < 0.0).collect();
        let lookup = |c: usize, a: usize, first: bool| {
            self.pairs
                .iter()
                .find(|p| p.cation == species[idx[c]].name && p.anion == species[idx[a]].name)
                .map_or(0.0, |p| if first { p.tau_ca_m } else { p.tau_m_ca })
        };
        let tau1 = cats
            .iter()
            .map(|&c| ans.iter().map(|&a| lookup(c, a, true)).collect())
            .collect();
        let tau2 = cats
            .iter()
            .map(|&c| ans.iter().map(|&a| lookup(c, a, false)).collect())
            .collect();
        Layout {
            z,
            n,
            cats,
            ans,
            tau1,
            tau2,
            idx,
        }
    }

    /// `n_total · Gᴱˣ/(n_total R T)` for `n_w` moles of water and ion mole
    /// numbers `ion_n` (same slots as `lay`).
    fn g_excess(&self, lay: &Layout, n_w: f64, ion_n: &[f64], a_phi: f64) -> f64 {
        let n_tot = n_w + ion_n.iter().sum::<f64>();
        let xm = n_w / n_tot;
        let xc: Vec<f64> = lay
            .cats
            .iter()
            .map(|&s| lay.z[s].abs() * ion_n[s] / n_tot)
            .collect();
        let xa: Vec<f64> = lay
            .ans
            .iter()
            .map(|&s| lay.z[s].abs() * ion_n[s] / n_tot)
            .collect();
        let sxc: f64 = xc.iter().sum();
        let sxa: f64 = xa.iter().sum();

        // Long-range Pitzer–Debye–Hückel term.
        let ix: f64 = 0.5
            * lay
                .z
                .iter()
                .zip(ion_n)
                .map(|(z, n)| z * z * n / n_tot)
                .sum::<f64>();
        let pdh = if ix > 0.0 {
            -n_tot * (1000.0 / M_WATER).sqrt() * 4.0 * a_phi * ix / self.rho
                * (1.0 + self.rho * ix.sqrt()).ln()
        } else {
            0.0
        };
        if sxc <= 0.0 || sxa <= 0.0 {
            return pdh;
        }

        let alpha = self.alpha;
        let g_ca_m = |ci: usize, ai: usize| (-alpha * lay.tau1[ci][ai]).exp();
        // Effective cation/anion–solvent parameters.
        let g_cm: Vec<f64> = (0..xc.len())
            .map(|ci| (0..xa.len()).map(|ai| xa[ai] / sxa * g_ca_m(ci, ai)).sum())
            .collect();
        let g_am: Vec<f64> = (0..xa.len())
            .map(|ai| (0..xc.len()).map(|ci| xc[ci] / sxc * g_ca_m(ci, ai)).sum())
            .collect();
        let tau_cm: Vec<f64> = g_cm.iter().map(|g: &f64| -g.ln() / alpha).collect();
        let tau_am: Vec<f64> = g_am.iter().map(|g: &f64| -g.ln() / alpha).collect();

        // Solvent-centred cell.
        let mut num = 0.0;
        let mut den = xm;
        for ci in 0..xc.len() {
            num += xc[ci] * g_cm[ci] * tau_cm[ci];
            den += xc[ci] * g_cm[ci];
        }
        for ai in 0..xa.len() {
            num += xa[ai] * g_am[ai] * tau_am[ai];
            den += xa[ai] * g_am[ai];
        }
        let mut lc = xm * num / den;

        // Cation-centred cells.
        for (ci, &xci) in xc.iter().enumerate() {
            let mut cell = 0.0;
            for (ai, &xai) in xa.iter().enumerate() {
                let gm = (-alpha * lay.tau2[ci][ai]).exp();
                cell += xai / sxa * xm * gm * lay.tau2[ci][ai] / (xm * gm + sxa);
            }
            lc += xci * cell;
        }
        // Anion-centred cells.
        for (ai, &xai) in xa.iter().enumerate() {
            let mut cell = 0.0;
            for (ci, &xci) in xc.iter().enumerate() {
                let gm = (-alpha * lay.tau2[ci][ai]).exp();
                cell += xci / sxc * xm * gm * lay.tau2[ci][ai] / (xm * gm + sxc);
            }
            lc += xai * cell;
        }
        pdh + n_tot * lc
    }

    /// Rational (mole-fraction) activity coefficients `ln γ` of water and of
    /// each ion slot, symmetric convention (`∂(nGᴱˣ/RT)/∂nᵢ`), for the
    /// per-kilogram-of-water mole numbers in `lay`.
    fn ln_gamma_x(&self, lay: &Layout, scale: f64, a_phi: f64) -> (f64, Vec<f64>) {
        let n_w = 1000.0 / M_WATER;
        let ion_n: Vec<f64> = lay.n.iter().map(|n| n * scale).collect();
        let hw = 1e-6 * n_w;
        let mut plus = ion_n.clone();
        let g_w = (self.g_excess(lay, n_w + hw, &plus, a_phi)
            - self.g_excess(lay, n_w - hw, &plus, a_phi))
            / (2.0 * hw);
        let mut out = vec![0.0; ion_n.len()];
        for s in 0..ion_n.len() {
            let h = 1e-5 * ion_n[s].max(1e-300);
            plus[s] = ion_n[s] + h;
            let up = self.g_excess(lay, n_w, &plus, a_phi);
            plus[s] = ion_n[s] - h;
            let dn = self.g_excess(lay, n_w, &plus, a_phi);
            plus[s] = ion_n[s];
            out[s] = (up - dn) / (2.0 * h);
        }
        (g_w, out)
    }

    /// Rational `ln γ_w` and `ln x_w` for the actual solution.
    fn water_terms(&self, lay: &Layout, a_phi: f64) -> (f64, f64) {
        let (g_w, _) = self.ln_gamma_x(lay, 1.0, a_phi);
        let n_w = 1000.0 / M_WATER;
        let x_w = n_w / (n_w + lay.n.iter().sum::<f64>());
        (g_w, x_w.ln())
    }

    /// Rational `ln γ_w` (see [`Enrtl::water_activity_terms`]) for external
    /// consistency checks: returns `(ln x_w, ln γ_w)` for the solution.
    pub fn water_activity_terms(&self, species: &[Species], t: f64) -> Result<(f64, f64)> {
        check_state(species, t)?;
        let lay = self.layout(species);
        let (gw, lnx) = self.water_terms(&lay, debye_huckel_a_phi(t));
        Ok((lnx, gw))
    }
}

impl ActivityModel for Enrtl {
    fn ln_gamma(&self, species: &[Species], t: f64) -> Result<Vec<f64>> {
        check_state(species, t)?;
        let lay = self.layout(species);
        let mut out = vec![0.0; species.len()];
        if lay.cats.is_empty() || lay.ans.is_empty() {
            return Ok(out);
        }
        let a_phi = debye_huckel_a_phi(t);
        let (_, g_sym) = self.ln_gamma_x(&lay, 1.0, a_phi);
        let (_, g_inf) = self.ln_gamma_x(&lay, 1e-8, a_phi);
        let total_m: f64 = lay.n.iter().sum();
        let to_molal = (1.0 + 1e-3 * M_WATER * total_m).ln();
        for (slot, &k) in lay.idx.iter().enumerate() {
            let v = g_sym[slot] - g_inf[slot] - to_molal;
            if !v.is_finite() {
                return Err(AqueousError::NonFinite("eNRTL activity coefficient"));
            }
            out[k] = v;
        }
        Ok(out)
    }

    fn osmotic_coefficient(&self, species: &[Species], t: f64) -> Result<f64> {
        check_state(species, t)?;
        let lay = self.layout(species);
        let total_m: f64 = lay.n.iter().sum();
        if total_m <= 0.0 {
            return Ok(1.0);
        }
        let (g_w, ln_x) = self.water_terms(&lay, debye_huckel_a_phi(t));
        Ok(-(ln_x + g_w) / (1e-3 * M_WATER * total_m))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nacl(m: f64) -> Vec<Species> {
        vec![Species::new("Na+", 1, m), Species::new("Cl-", -1, m)]
    }

    fn model() -> Enrtl {
        let mut e = Enrtl::new();
        e.add_pair(EnrtlPair {
            cation: "Na+".into(),
            anion: "Cl-".into(),
            tau_ca_m: 8.885,
            tau_m_ca: -4.549,
        });
        e
    }

    #[test]
    fn dilute_limit_follows_debye_huckel() {
        let g = Enrtl::new().ln_gamma(&nacl(1e-4), 298.15).unwrap();
        let limiting = -crate::solvent::debye_huckel_a_gamma(298.15) * 1e-2;
        assert!(
            (g[0] - limiting).abs() / limiting.abs() < 0.06,
            "{} vs {limiting}",
            g[0]
        );
    }

    #[test]
    fn zero_taus_have_no_local_composition_contribution() {
        // With τ = 0, only the PDH term acts: γ < 1 and decreasing with m.
        let a = Enrtl::new().ln_gamma(&nacl(0.1), 298.15).unwrap()[0];
        let b = Enrtl::new().ln_gamma(&nacl(0.5), 298.15).unwrap()[0];
        assert!(a < 0.0 && b < a);
    }

    #[test]
    fn infinite_dilution_reference_is_unity() {
        let g = model().ln_gamma(&nacl(1e-7), 298.15).unwrap();
        assert!(g[0].abs() < 1e-2 && g[1].abs() < 1e-2, "{g:?}");
    }

    #[test]
    fn gibbs_duhem_holds_numerically() {
        // n_w Δlnγ_w + Σ n_i Δlnγ_i = 0 for the rational coefficients when
        // the solute amount changes (symmetric convention).
        let m = model();
        let t = 298.15;
        let a_phi = debye_huckel_a_phi(t);
        let n_w = 1000.0 / M_WATER;
        let eval = |c: f64| {
            let lay = m.layout(&nacl(c));
            m.ln_gamma_x(&lay, 1.0, a_phi)
        };
        let (c0, dc) = (0.8, 1e-4);
        let (w0, i0) = eval(c0);
        let (w1, i1) = eval(c0 + dc);
        let lhs = n_w * (w1 - w0) + c0 * ((i1[0] - i0[0]) + (i1[1] - i0[1]));
        let scale = n_w * (w1 - w0).abs();
        assert!(
            lhs.abs() < 1e-3 * scale.max(1e-9),
            "GD residual {lhs}, scale {scale}"
        );
    }

    #[test]
    fn osmotic_coefficient_is_finite_and_unity_when_dilute() {
        let phi = model().osmotic_coefficient(&nacl(1e-6), 298.15).unwrap();
        assert!((phi - 1.0).abs() < 1e-2);
        assert!(model()
            .osmotic_coefficient(&nacl(1.0), 298.15)
            .unwrap()
            .is_finite());
    }
}
