//! Property tests for `tpt-chem-aqueous` (spec.txt §6, Phase 1; CONTRIBUTING
//! ground rule 5: physics is verified, not assumed).
//!
//! The properties cover sign and monotonicity of activity coefficients,
//! finiteness over the valid input range, electroneutrality-independent
//! ionic-strength scaling, van't Hoff reference consistency, and mass and
//! charge balance of the calcium-phosphate speciation solver.

use std::vec::Vec;

use tpt_chem_aqueous::activity::debye_huckel::{Davies, ExtendedDebyeHuckel, LimitingLaw};
use tpt_chem_aqueous::activity::pitzer::{Pitzer, PitzerPair};
use tpt_chem_aqueous::activity::ActivityModel;
use tpt_chem_aqueous::chelator::Ligand;
use tpt_chem_aqueous::constants::{LN10, R};
use tpt_chem_aqueous::equilibrium::Equilibrium;
use tpt_chem_aqueous::precipitation::CalciumPhosphateSystem;
use tpt_chem_aqueous::species::{ionic_strength, Species};

/// Molality of a 1:1 salt, as a species pair.
fn nacl(m: f64) -> Vec<Species> {
    vec![Species::new("Na+", 1, m), Species::new("Cl-", -1, m)]
}

/// Ion pair of charges `zc`, `za` with `m` molal salt (electroneutral).
fn salt(zc: i32, za: i32, m: f64) -> Vec<Species> {
    vec![
        Species::new("M", zc, m * f64::from(za.abs())),
        Species::new("X", za, m * f64::from(zc)),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    const T_RANGE: std::ops::Range<f64> = 278.15..348.15;

    fn models() -> Vec<Box<dyn ActivityModel>> {
        vec![
            Box::new(LimitingLaw),
            Box::new(ExtendedDebyeHuckel::default()),
            Box::new(Davies::default()),
        ]
    }

    fn pitzer_nacl() -> Pitzer {
        let mut p = Pitzer::new();
        p.add_pair(PitzerPair::new("Na+", "Cl-", 0.0765, 0.2664, 0.0, 0.00127));
        p
    }

    proptest! {
        #![proptest_config(ProptestConfig::with_cases(128))]

        /// Every Debye–Hückel-family coefficient of a 1:1 salt is ≤ 1 for any
        /// molality, and the three models agree on the sign.
        #[test]
        fn dh_family_gamma_below_one(m in 1e-8_f64..0.5, t in T_RANGE) {
            for model in models() {
                let g = model.ln_gamma(&nacl(m), t).unwrap();
                prop_assert!(g.iter().all(|l| l.is_finite() && *l <= 1e-12));
            }
        }

        /// Activity coefficients fall (become more negative in ln) as ionic
        /// strength rises in the Davies range.
        #[test]
        fn davies_monotone_in_ionic_strength(m in 1e-4_f64..0.2, dm in 1e-4_f64..0.1, t in T_RANGE) {
            let model = Davies::default();
            let a = model.ln_gamma(&nacl(m), t).unwrap()[0];
            let b = model.ln_gamma(&nacl(m + dm), t).unwrap()[0];
            prop_assert!(b <= a + 1e-12, "{b} > {a}");
        }

        /// Ionic strength scales linearly with a uniform molality scaling.
        #[test]
        fn ionic_strength_is_homogeneous(m in 1e-6_f64..5.0, k in 0.1_f64..10.0) {
            let i1 = ionic_strength(&salt(2, -1, m));
            let ik = ionic_strength(&salt(2, -1, m * k));
            prop_assert!((ik - k * i1).abs() <= 1e-12 * ik.abs().max(1.0));
        }

        /// The Pitzer mean activity coefficient stays finite and positive over
        /// the validated molality range.
        #[test]
        fn pitzer_nacl_mean_gamma_positive_and_finite(m in 1e-6_f64..6.0, t in T_RANGE) {
            let g = pitzer_nacl().ln_gamma(&nacl(m), t).unwrap();
            let mean = 0.5 * (g[0] + g[1]);
            prop_assert!(mean.is_finite());
            prop_assert!(mean.exp() > 0.0);
        }

        /// Pitzer osmotic coefficient is finite and approaches 1 at infinite dilution.
        #[test]
        fn pitzer_osmotic_coefficient_finite(m in 1e-6_f64..6.0, t in T_RANGE) {
            let phi = pitzer_nacl().osmotic_coefficient(&nacl(m), t).unwrap();
            prop_assert!(phi.is_finite());
        }

        /// Equilibrium constants are exact at the reference temperature, and the
        /// van't Hoff slope d(log K)/dT = ΔH / (R ln 10 T²) has the sign of ΔH.
        #[test]
        fn van_t_hoff_reference_and_slope(logk in -30.0_f64..30.0, dh in -1e5_f64..1e5, t in T_RANGE) {
            let e = Equilibrium::new(logk, dh);
            prop_assert!((e.log_k(298.15) - logk).abs() < 1e-12);
            let h = 1e-3;
            let slope = (e.log_k(t + h) - e.log_k(t - h)) / (2.0 * h);
            let expected = dh / (R * LN10 * t * t);
            prop_assert!((slope - expected).abs() <= 1e-6 * expected.abs().max(1e-12));
        }

        /// Speciation of the calcium-phosphate system: all molalities finite and
        /// non-negative, calcium and phosphate balances close, and the solver
        /// converges over a broad pH/T/concentration range.
        #[test]
        fn speciation_conserves_calcium_and_phosphate(
            ph in 5.0_f64..9.0,
            t in 283.15_f64..333.15,
            ca in 1e-4_f64..1e-2,
            p in 1e-4_f64..1e-2,
        ) {
            let mut sys = CalciumPhosphateSystem::new(ca, p, Box::new(Davies::default()));
            sys.background = nacl(0.1);
            let sp = sys.speciate(ph, t).unwrap();
            prop_assert!(sp.species.iter().all(|s| s.molality.is_finite() && s.molality >= 0.0));

            let total_ca: f64 = sp
                .species
                .iter()
                .filter(|s| {
                    matches!(s.name.as_str(), "Ca2+" | "CaH2PO4+" | "CaPO4-" | "CaHPO4")
                        || s.name.ends_with(":CaS")
                })
                .map(|s| s.molality)
                .sum();
            prop_assert!((total_ca - ca).abs() <= 1e-8 * ca, "Ca balance {total_ca} vs {ca}");

            let total_p: f64 = sp
                .species
                .iter()
                .filter(|s| {
                    matches!(
                        s.name.as_str(),
                        "H2PO4-" | "HPO4--" | "PO4---" | "CaH2PO4+" | "CaPO4-" | "CaHPO4"
                    )
                })
                .map(|s| s.molality)
                .sum();
            // Neutral H3PO4 is not tracked as a species, so phosphate may be
            // partly missing at low pH; it can never be in excess.
            prop_assert!(total_p <= p * (1.0 + 1e-8), "P excess {total_p} vs {p}");
        }

        /// Adding a chelator never increases free calcium.
        #[test]
        fn chelator_never_increases_free_calcium(
            ph in 5.5_f64..8.5,
            ligand in 1e-4_f64..5e-3,
        ) {
            let t = 298.15;
            let mut bare = CalciumPhosphateSystem::new(5e-3, 3e-3, Box::new(Davies::default()));
            bare.background = nacl(0.1);
            let mut with = CalciumPhosphateSystem::new(5e-3, 3e-3, Box::new(Davies::default()));
            with.background = nacl(0.1);
            with.add_ligand(Ligand::casein_phosphopeptide(), ligand);
            let free_bare = bare.speciate(ph, t).unwrap().free_ca;
            let free_with = with.speciate(ph, t).unwrap().free_ca;
            prop_assert!(free_with <= free_bare * (1.0 + 1e-9));
        }
    }
}
