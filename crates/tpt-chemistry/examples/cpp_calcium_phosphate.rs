//! Calcium-phosphate precipitation onset (chalkiness) with and without a
//! casein-phosphopeptide-like chelator, across temperature.
//!
//! Uses illustrative constants (see `tpt_chem_aqueous::phosphate` and
//! `tpt_chem_aqueous::chelator`): the output shows the trend, not a
//! calibrated prediction for any particular milk or formulation.
//!
//! Run: `cargo run -p tpt-chemistry --features aqueous --example cpp_calcium_phosphate`

use tpt_chemistry::aqueous::activity::debye_huckel::Davies;
use tpt_chemistry::aqueous::chelator::Ligand;
use tpt_chemistry::aqueous::precipitation::{CalciumPhosphateSystem, Onset};
use tpt_chemistry::aqueous::species::Species;

fn system(cpp_molal: Option<f64>) -> CalciumPhosphateSystem {
    // Illustrative totals: 5 mM Ca and 3 mM inorganic phosphate in 0.1 m NaCl.
    // Real milk is far more concentrated and needs measured constants.
    let mut sys = CalciumPhosphateSystem::new(5e-3, 3e-3, Box::new(Davies::default()));
    sys.background = vec![Species::new("Na+", 1, 0.1), Species::new("Cl-", -1, 0.1)];
    if let Some(m) = cpp_molal {
        sys.add_ligand(Ligand::casein_phosphopeptide(), m);
    }
    sys
}

fn describe(onset: Onset) -> String {
    match onset {
        Onset::At(ph) => format!("{ph:.2}"),
        Onset::Always => "precipitates at all pH".to_string(),
        Onset::Never => "no precipitation in pH 4–9".to_string(),
    }
}

fn main() {
    let temperatures_c = [4.0, 25.0, 45.0, 60.0, 75.0];
    let t_grid: Vec<f64> = temperatures_c.iter().map(|c| c + 273.15).collect();

    let bare = system(None).onset_curve(&t_grid, 4.0, 9.0).unwrap();
    let with_cpp = system(Some(3e-3)).onset_curve(&t_grid, 4.0, 9.0).unwrap();

    println!("Onset pH of calcium phosphate supersaturation (SI = 0), Ca 5 mM, P 3 mM, pH 4–9:");
    println!("{:>8}  {:>12}  {:>22}", "T / °C", "no chelator", "with CPP (3 mM)");
    for ((c, (_, b)), (_, w)) in temperatures_c.iter().zip(&bare).zip(&with_cpp) {
        println!("{c:>8.0}  {:>12}  {:>22}", describe(*b), describe(*w));
    }
    println!(
        "\n(the chelator sequesters Ca²⁺, so precipitation is pushed to higher pH; \
         a higher si_threshold would model a nucleation barrier)"
    );
}
