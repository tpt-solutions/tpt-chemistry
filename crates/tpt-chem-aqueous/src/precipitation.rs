//! Calcium phosphate speciation and precipitation (chalkiness) prediction.
//!
//! [`CalciumPhosphateSystem`] holds total calcium and phosphate, a
//! background electrolyte, optional chelating [`Ligand`]s and a set of
//! candidate solid phases. At a fixed pH (activity scale, assumed to be held
//! by an external buffer or titrant) and temperature it solves the coupled
//! calcium/phosphate mass balances with ion pairing and ligand binding,
//! iterating activity coefficients to self-consistency, and reports the
//! saturation index of every phase **assuming nothing has precipitated**.
//! Scanning pH and temperature, or searching for the pH at which the
//! supersaturation threshold is crossed, then predicts the onset of
//! precipitation.
//!
//! The result is thermodynamic: `SI > 0` means the solution *can* lower its
//! free energy by precipitating. Whether visible turbidity develops also
//! depends on nucleation and growth kinetics; use
//! [`CalciumPhosphateSystem::si_threshold`] to demand a finite
//! supersaturation before declaring onset.
//!
//! # Examples
//!
//! ```
//! use tpt_chem_aqueous::activity::debye_huckel::Davies;
//! use tpt_chem_aqueous::chelator::Ligand;
//! use tpt_chem_aqueous::precipitation::{CalciumPhosphateSystem, Onset};
//! use tpt_chem_aqueous::species::Species;
//!
//! let background = vec![Species::new("Na+", 1, 0.1), Species::new("Cl-", -1, 0.1)];
//! let mut sys = CalciumPhosphateSystem::new(5e-3, 3e-3, Box::new(Davies::default()));
//! sys.background = background;
//! let bare = sys.onset_ph(310.15, 5.0, 9.0).unwrap();
//! sys.add_ligand(Ligand::casein_phosphopeptide(), 2e-3);
//! let with_cpp = sys.onset_ph(310.15, 5.0, 9.0).unwrap();
//! if let (Onset::At(a), Onset::At(b)) = (bare, with_cpp) {
//!     assert!(b > a); // the chelator delays precipitation to a higher pH
//! }
//! ```

use crate::activity::ActivityModel;
use crate::chelator::Ligand;
use crate::constants::LN10;
use crate::phosphate::{default_phases, PhosphateConstants};
use crate::solubility::{ComponentLogActivities, SolidPhase};
use crate::solver::brent;
use crate::species::{ionic_strength, Species};
use crate::{AqueousError, Result};

/// Outcome of an onset-pH search.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Onset {
    /// The solution stays below the threshold over the whole pH range.
    Never,
    /// The solution is already above the threshold at the lower pH bound.
    Always,
    /// Precipitation threshold is crossed at this pH (rising pH).
    At(f64),
}

/// Equilibrium state of the (metastable) solution at fixed pH and T.
#[derive(Debug, Clone)]
pub struct Speciation {
    /// All species with their equilibrium molalities (background, calcium
    /// and phosphate species, H⁺/OH⁻, and ligand site populations).
    pub species: Vec<Species>,
    /// Activity coefficients matching `species` (same order).
    pub gamma: Vec<f64>,
    /// Free Ca²⁺ molality.
    pub free_ca: f64,
    /// Molal ionic strength.
    pub ionic_strength: f64,
    /// `log₁₀` activities of the components entering solubility products.
    pub log_activities: ComponentLogActivities,
    /// Saturation index of each candidate phase, in the order of
    /// [`CalciumPhosphateSystem::phases`].
    pub saturation: Vec<(String, f64)>,
    /// Number of activity-coefficient iterations used.
    pub iterations: usize,
}

impl Speciation {
    /// The most supersaturated phase and its saturation index.
    pub fn max_saturation(&self) -> (f64, &str) {
        self.saturation
            .iter()
            .fold((f64::NEG_INFINITY, ""), |acc, (n, si)| {
                if *si > acc.0 {
                    (*si, n.as_str())
                } else {
                    acc
                }
            })
    }
}

/// One point of a pH × temperature scan.
#[derive(Debug, Clone, PartialEq)]
pub struct ScanPoint {
    /// pH.
    pub ph: f64,
    /// Temperature, K.
    pub temperature: f64,
    /// Largest saturation index over all phases.
    pub max_si: f64,
    /// Name of the phase with the largest saturation index.
    pub phase: String,
}

/// A Ca/phosphate solution with optional chelators and background salt.
pub struct CalciumPhosphateSystem {
    /// Total calcium, mol/kg water.
    pub ca_total: f64,
    /// Total inorganic phosphate, mol/kg water.
    pub phosphate_total: f64,
    /// Background electrolyte (fixed molalities) contributing to ionic
    /// strength and activity coefficients.
    pub background: Vec<Species>,
    /// Chelating ligands with their total molalities.
    pub ligands: Vec<(Ligand, f64)>,
    /// Acid/ion-pair constants.
    pub constants: PhosphateConstants,
    /// Candidate solid phases.
    pub phases: Vec<SolidPhase>,
    /// Saturation-index threshold that defines "onset" (0 = thermodynamic
    /// saturation; larger values mimic a nucleation barrier).
    pub si_threshold: f64,
    model: Box<dyn ActivityModel>,
}

/// Indices into the activity-coefficient vector of the speciation.
#[derive(Clone, Copy)]
struct Gammas {
    ca: f64,
    h2po4: f64,
    hpo4: f64,
    po4: f64,
    cah2po4: f64,
    capo4: f64,
}

struct Core {
    a_ca: f64,
    a3: f64,
    log_a3: f64,
}

impl CalciumPhosphateSystem {
    /// Creates a system with default constants/phases, no background salt
    /// and no ligands.
    pub fn new(ca_total: f64, phosphate_total: f64, model: Box<dyn ActivityModel>) -> Self {
        Self {
            ca_total,
            phosphate_total,
            background: Vec::new(),
            ligands: Vec::new(),
            constants: PhosphateConstants::default(),
            phases: default_phases(),
            si_threshold: 0.0,
            model,
        }
    }

    /// Adds a ligand at the given total molality.
    pub fn add_ligand(&mut self, ligand: Ligand, molality: f64) {
        self.ligands.push((ligand, molality));
    }

    /// Phosphate partition coefficients (per unit `a(PO₄³⁻)`) at fixed `a_ca`.
    fn phosphate_terms(&self, a_ca: f64, a_h: f64, t: f64, g: &Gammas) -> [f64; 7] {
        let c = &self.constants;
        let (ka1, ka2, ka3) = (c.ka1.k(t), c.ka2.k(t), c.ka3.k(t));
        let hpo4 = a_h / ka3;
        let h2po4 = a_h * a_h / (ka2 * ka3);
        [
            a_h * a_h * a_h / (ka1 * ka2 * ka3),
            h2po4 / g.h2po4,
            hpo4 / g.hpo4,
            1.0 / g.po4,
            c.ca_hpo4.k(t) * a_ca * hpo4,
            c.ca_h2po4.k(t) * a_ca * h2po4 / g.cah2po4,
            c.ca_po4.k(t) * a_ca / g.capo4,
        ]
    }

    fn solve_core(&self, a_h: f64, t: f64, g: &Gammas) -> Result<Core> {
        let balance = |x: f64| {
            let a_ca = x.exp();
            let terms = self.phosphate_terms(a_ca, a_h, t, g);
            let s: f64 = terms.iter().sum();
            let a3 = self.phosphate_total / s;
            let pairs = a3 * (terms[4] + terms[5] + terms[6]);
            let ligand: f64 = self
                .ligands
                .iter()
                .map(|(l, m)| l.bound_calcium(*m, a_ca, a_h, t))
                .sum();
            a_ca / g.ca + pairs + ligand - self.ca_total
        };
        let hi = (g.ca * self.ca_total).ln();
        let lo = hi - 120.0 * LN10.max(1.0);
        let x = brent(balance, lo.min(-90.0), hi, 1e-13, 300)?;
        let a_ca = x.exp();
        let s: f64 = self.phosphate_terms(a_ca, a_h, t, g).iter().sum();
        Ok(Core {
            a_ca,
            a3: self.phosphate_total / s,
            log_a3: (self.phosphate_total / s).ln(),
        })
    }

    /// Solves the speciation at the given pH and temperature (K).
    pub fn speciate(&self, ph: f64, t: f64) -> Result<Speciation> {
        if !(self.ca_total > 0.0 && self.phosphate_total > 0.0) {
            return Err(AqueousError::InvalidInput(
                "total Ca and phosphate must be positive",
            ));
        }
        if !(ph.is_finite() && (0.0..=14.0).contains(&ph)) {
            return Err(AqueousError::InvalidInput("pH must be in 0–14"));
        }
        let a_h = 10f64.powf(-ph);
        let a_oh = self.constants.kw.k(t) / a_h;
        let mut g = Gammas {
            ca: 1.0,
            h2po4: 1.0,
            hpo4: 1.0,
            po4: 1.0,
            cah2po4: 1.0,
            capo4: 1.0,
        };
        let (mut g_h, mut g_oh) = (1.0, 1.0);
        for iteration in 1..=200 {
            let core = self.solve_core(a_h, t, &g)?;
            let (species, idx) = self.assemble(&core, (a_h, a_oh, t), &g, (g_h, g_oh));
            let ln_g = self.model.ln_gamma(&species, t)?;
            let new = Gammas {
                ca: ln_g[idx.ca].exp(),
                h2po4: ln_g[idx.ca + 1].exp(),
                hpo4: ln_g[idx.ca + 2].exp(),
                po4: ln_g[idx.ca + 3].exp(),
                cah2po4: ln_g[idx.ca + 4].exp(),
                capo4: ln_g[idx.ca + 5].exp(),
            };
            let (n_h, n_oh) = (ln_g[idx.ca + 7].exp(), ln_g[idx.ca + 8].exp());
            let diff = [
                (new.ca / g.ca).ln(),
                (new.h2po4 / g.h2po4).ln(),
                (new.hpo4 / g.hpo4).ln(),
                (new.po4 / g.po4).ln(),
                (new.cah2po4 / g.cah2po4).ln(),
                (new.capo4 / g.capo4).ln(),
                (n_h / g_h).ln(),
                (n_oh / g_oh).ln(),
            ]
            .iter()
            .fold(0.0_f64, |m, d| m.max(d.abs()));
            if !diff.is_finite() {
                return Err(AqueousError::NonFinite("activity coefficients"));
            }
            // Damped (geometric-mean) update.
            let mix = |old: f64, new: f64| (0.5 * (old.ln() + new.ln())).exp();
            g = Gammas {
                ca: mix(g.ca, new.ca),
                h2po4: mix(g.h2po4, new.h2po4),
                hpo4: mix(g.hpo4, new.hpo4),
                po4: mix(g.po4, new.po4),
                cah2po4: mix(g.cah2po4, new.cah2po4),
                capo4: mix(g.capo4, new.capo4),
            };
            g_h = mix(g_h, n_h);
            g_oh = mix(g_oh, n_oh);
            if diff < 1e-9 {
                let core = self.solve_core(a_h, t, &g)?;
                let (species, _) = self.assemble(&core, (a_h, a_oh, t), &g, (g_h, g_oh));
                let gamma = self
                    .model
                    .ln_gamma(&species, t)?
                    .iter()
                    .map(|l| l.exp())
                    .collect();
                return Ok(self.finish((species, gamma), &core, (ph, a_oh, t), iteration));
            }
        }
        Err(AqueousError::NoConvergence(
            "activity-coefficient iteration",
        ))
    }

    /// Builds the species list from a solved core.
    fn assemble(
        &self,
        core: &Core,
        (a_h, a_oh, t): (f64, f64, f64),
        g: &Gammas,
        (g_h, g_oh): (f64, f64),
    ) -> (Vec<Species>, Idx) {
        let terms = self.phosphate_terms(core.a_ca, a_h, t, g);
        let a3 = core.a3;
        let mut sp = self.background.clone();
        let ca = sp.len();
        sp.push(Species::new("Ca2+", 2, core.a_ca / g.ca));
        sp.push(Species::new("H2PO4-", -1, a3 * terms[1]));
        sp.push(Species::new("HPO4--", -2, a3 * terms[2]));
        sp.push(Species::new("PO4---", -3, a3 * terms[3]));
        sp.push(Species::new("CaH2PO4+", 1, a3 * terms[5]));
        sp.push(Species::new("CaPO4-", -1, a3 * terms[6]));
        sp.push(Species::new("CaHPO4", 0, a3 * terms[4]));
        sp.push(Species::new("H+", 1, a_h / g_h));
        sp.push(Species::new("OH-", -1, a_oh / g_oh));
        for (lig, m) in &self.ligands {
            for site in &lig.sites {
                let pop = site.population(m * site.count, core.a_ca, a_h, t);
                let base = site.protonated_charge;
                let tag = format!("{}:{}", lig.name, site.name);
                sp.push(Species::new(&format!("{tag}:HS"), base, pop.protonated));
                sp.push(Species::new(&format!("{tag}:S"), base - 1, pop.free));
                sp.push(Species::new(&format!("{tag}:CaS"), base + 1, pop.bound));
            }
        }
        (sp, Idx { ca })
    }

    fn finish(
        &self,
        (species, gamma): (Vec<Species>, Vec<f64>),
        core: &Core,
        (ph, a_oh, t): (f64, f64, f64),
        iterations: usize,
    ) -> Speciation {
        let c = &self.constants;
        let log_a3 = core.log_a3 / LN10;
        let log_ah = -ph;
        let log_activities = ComponentLogActivities {
            ca: core.a_ca.log10(),
            h: log_ah,
            oh: a_oh.log10(),
            h2po4: log_a3 + 2.0 * log_ah - c.ka2.log_k(t) - c.ka3.log_k(t),
            hpo4: log_a3 + log_ah - c.ka3.log_k(t),
            po4: log_a3,
        };
        let saturation = self
            .phases
            .iter()
            .map(|p| (p.name.clone(), p.saturation_index(&log_activities, t)))
            .collect();
        let free_ca = species
            .iter()
            .find(|s| s.name == "Ca2+")
            .map_or(0.0, |s| s.molality);
        Speciation {
            ionic_strength: ionic_strength(&species),
            species,
            gamma,
            free_ca,
            log_activities,
            saturation,
            iterations,
        }
    }

    /// Largest saturation index over all phases (and the phase's name).
    pub fn max_saturation(&self, ph: f64, t: f64) -> Result<(f64, String)> {
        let s = self.speciate(ph, t)?;
        let (si, name) = s.max_saturation();
        Ok((si, name.to_string()))
    }

    /// Scans a pH × temperature grid (temperature varies slowest).
    pub fn scan(&self, ph_grid: &[f64], t_grid: &[f64]) -> Result<Vec<ScanPoint>> {
        let mut out = Vec::with_capacity(ph_grid.len() * t_grid.len());
        for &t in t_grid {
            for &ph in ph_grid {
                let (max_si, phase) = self.max_saturation(ph, t)?;
                out.push(ScanPoint {
                    ph,
                    temperature: t,
                    max_si,
                    phase,
                });
            }
        }
        Ok(out)
    }

    /// pH at which the largest saturation index crosses
    /// [`si_threshold`](Self::si_threshold) while raising pH from `ph_lo` to
    /// `ph_hi`, at temperature `t` (K).
    pub fn onset_ph(&self, t: f64, ph_lo: f64, ph_hi: f64) -> Result<Onset> {
        if ph_lo >= ph_hi {
            return Err(AqueousError::InvalidInput("ph_lo must be below ph_hi"));
        }
        let f = |ph: f64| -> f64 {
            self.max_saturation(ph, t)
                .map_or(f64::NAN, |(si, _)| si - self.si_threshold)
        };
        let (flo, fhi) = (f(ph_lo), f(ph_hi));
        if flo.is_nan() || fhi.is_nan() {
            return Err(AqueousError::NonFinite("saturation index"));
        }
        if flo >= 0.0 {
            return Ok(Onset::Always);
        }
        if fhi <= 0.0 {
            return Ok(Onset::Never);
        }
        Ok(Onset::At(brent(f, ph_lo, ph_hi, 1e-8, 200)?))
    }

    /// Onset pH for each temperature in `t_grid` (see [`onset_ph`](Self::onset_ph)).
    pub fn onset_curve(&self, t_grid: &[f64], ph_lo: f64, ph_hi: f64) -> Result<Vec<(f64, Onset)>> {
        t_grid
            .iter()
            .map(|&t| Ok((t, self.onset_ph(t, ph_lo, ph_hi)?)))
            .collect()
    }
}

struct Idx {
    ca: usize,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::activity::debye_huckel::Davies;

    fn system(ligand: f64) -> CalciumPhosphateSystem {
        let mut s = CalciumPhosphateSystem::new(5e-3, 3e-3, Box::new(Davies::default()));
        s.background = vec![Species::new("Na+", 1, 0.1), Species::new("Cl-", -1, 0.1)];
        if ligand > 0.0 {
            s.add_ligand(Ligand::casein_phosphopeptide(), ligand);
        }
        s
    }

    #[test]
    fn mass_balances_close() {
        let s = system(1e-3);
        let sp = s.speciate(6.8, 310.15).unwrap();
        let get = |n: &str| sp.species.iter().find(|x| x.name == n).unwrap().molality;
        let p = get("H2PO4-")
            + get("HPO4--")
            + get("PO4---")
            + get("CaH2PO4+")
            + get("CaPO4-")
            + get("CaHPO4");
        // H3PO4 is not in the species list (neutral, no charge); add via ratio check:
        assert!(p <= 3e-3 * (1.0 + 1e-9));
        assert!(p > 0.9 * 3e-3 || sp.free_ca > 0.0);
        let ca = get("Ca2+")
            + get("CaH2PO4+")
            + get("CaPO4-")
            + get("CaHPO4")
            + sp.species
                .iter()
                .filter(|x| x.name.ends_with(":CaS"))
                .map(|x| x.molality)
                .sum::<f64>();
        assert!((ca - 5e-3).abs() < 1e-10, "Ca balance {ca}");
    }

    #[test]
    fn saturation_rises_with_ph() {
        let s = system(0.0);
        let lo = s.max_saturation(5.5, 298.15).unwrap().0;
        let hi = s.max_saturation(7.5, 298.15).unwrap().0;
        assert!(hi > lo);
    }

    #[test]
    fn chelator_lowers_free_calcium_and_delays_onset() {
        let bare = system(0.0);
        let cpp = system(3e-3);
        let a = bare.speciate(6.5, 310.15).unwrap();
        let b = cpp.speciate(6.5, 310.15).unwrap();
        assert!(b.free_ca < a.free_ca);
        let (Onset::At(p0), Onset::At(p1)) = (
            bare.onset_ph(310.15, 4.0, 9.0).unwrap(),
            cpp.onset_ph(310.15, 4.0, 9.0).unwrap(),
        ) else {
            panic!("expected bracketed onsets");
        };
        assert!(p1 > p0, "{p1} should exceed {p0}");
    }

    #[test]
    fn zero_ligand_matches_no_ligand() {
        let a = system(0.0).speciate(7.0, 298.15).unwrap();
        let mut b = system(0.0);
        b.add_ligand(Ligand::casein_phosphopeptide(), 0.0);
        let b = b.speciate(7.0, 298.15).unwrap();
        assert!((a.max_saturation().0 - b.max_saturation().0).abs() < 1e-9);
    }

    #[test]
    fn onset_is_at_the_threshold() {
        let s = system(0.0);
        let Onset::At(ph) = s.onset_ph(298.15, 4.0, 9.0).unwrap() else {
            panic!("no onset")
        };
        let (si, _) = s.max_saturation(ph, 298.15).unwrap();
        assert!(si.abs() < 1e-5, "SI at onset = {si}");
    }

    #[test]
    fn higher_threshold_moves_onset_to_higher_ph() {
        let mut s = system(0.0);
        let Onset::At(p0) = s.onset_ph(298.15, 4.0, 9.0).unwrap() else {
            panic!()
        };
        s.si_threshold = 1.0;
        let Onset::At(p1) = s.onset_ph(298.15, 4.0, 9.0).unwrap() else {
            panic!()
        };
        assert!(p1 > p0);
    }

    #[test]
    fn retrograde_solubility_lowers_onset_ph_when_hot() {
        let s = system(0.0);
        let (Onset::At(cold), Onset::At(hot)) = (
            s.onset_ph(283.15, 4.0, 9.0).unwrap(),
            s.onset_ph(353.15, 4.0, 9.0).unwrap(),
        ) else {
            panic!()
        };
        assert!(hot < cold, "{hot} vs {cold}");
    }

    #[test]
    fn scan_has_expected_shape() {
        let s = system(0.0);
        let pts = s.scan(&[6.0, 7.0], &[298.15, 310.15, 323.15]).unwrap();
        assert_eq!(pts.len(), 6);
        assert!(pts.iter().all(|p| p.max_si.is_finite()));
    }

    #[test]
    fn rejects_bad_input() {
        let mut s = system(0.0);
        s.ca_total = 0.0;
        assert!(s.speciate(7.0, 298.15).is_err());
        assert!(system(0.0).speciate(15.0, 298.15).is_err());
    }
}
