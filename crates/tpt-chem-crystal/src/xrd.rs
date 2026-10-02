//! Powder X-ray diffraction simulation from structure factors.
//!
//! Given a lattice, atomic positions (fractional) with atomic numbers, and
//! X-ray wavelength, the simulator:
//! 1. enumerates Miller indices up to a d-spacing cutoff,
//! 2. computes the structure factor `F(hkl) = Σ_j f_j(G) exp(2πi g·r_j)`
//!    with Cromer-Mann atomic form factors,
//! 3. keeps allowed reflections (`|F|² > cutoff`), and reports the 2θ
//!    peak positions with intensities `|F|² × Lorentz-polarization`.

use std::vec::Vec;

use crate::bravais::Lattice;
use crate::miller::Miller;
use tpt_chem_core::num;
use tpt_chem_core::vec3::Vec3;

/// One simulated powder-diffraction peak.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Reflection {
    /// Miller index.
    pub miller: Miller,
    /// d-spacing (Å).
    pub d: f64,
    /// 2θ angle in degrees.
    pub two_theta: f64,
    /// Relative intensity normalized to the strongest peak (0..1).
    pub intensity: f64,
}

/// Cromer-Mann 9-parameter form factor coefficients per element.
struct CromerMann {
    /// a₁..a₄.
    a: [f64; 4],
    /// b₁..b₄.
    b: [f64; 4],
    /// c constant.
    c: f64,
}

/// Cromer-Mann coefficients for common elements (International Tables).
fn cromer_mann(z: u8) -> Option<CromerMann> {
    let cm = |a: [f64; 4], b: [f64; 4], c: f64| Some(CromerMann { a, b, c });
    match z {
        1 => cm(
            [0.493002, 0.322912, 0.140191, 0.04081],
            [10.5109, 26.1257, 3.14236, 57.7997],
            0.003038,
        ),
        6 => cm(
            [2.31, 1.0202, 1.5886, 0.865],
            [20.8439, 10.2075, 0.5687, 51.6512],
            0.2156,
        ),
        7 => cm(
            [12.2126, 3.1322, 2.0125, 1.1663],
            [0.0057, 9.8933, 28.9975, 0.5826],
            -11.529,
        ),
        8 => cm(
            [3.0485, 2.2868, 1.5463, 0.867],
            [13.2771, 5.7011, 0.3239, 32.9089],
            0.2508,
        ),
        14 => cm(
            [6.2915, 3.0353, 1.9891, 1.541],
            [2.4386, 32.3337, 0.6785, 81.6937],
            1.1407,
        ),
        _ => None,
    }
}


/// Atomic form factor `f(G)` for an element at scattering magnitude
/// `s = sinθ/λ` (Cromer-Mann): `f = Σ aᵢ exp(−bᵢ s²) + c`.
pub fn form_factor(z: u8, s: f64) -> f64 {
    match cromer_mann(z) {
        Some(cm) => cm
            .a
            .iter()
            .zip(&cm.b)
            .map(|(&a, &b)| a * num::exp(-b * s * s))
            .sum::<f64>()
            + cm.c,
        None => z as f64, // fallback: atomic number
    }
}

/// Atomic positions in a structure: atomic numbers + fractional coords.
#[derive(Clone, Debug)]
pub struct CrystalStructure {
    /// Atomic numbers per site.
    pub atomic_numbers: Vec<u8>,
    /// Fractional coordinates per site.
    pub frac_positions: Vec<Vec3>,
    /// Debye-Waller isotropic B factor (default 1.0 if zero).
    pub b_factor: f64,
}

impl CrystalStructure {
    /// A structure from parallel arrays.
    pub fn new(atomic_numbers: Vec<u8>, frac_positions: Vec<Vec3>) -> Self {
        CrystalStructure {
            atomic_numbers,
            frac_positions,
            b_factor: 1.0,
        }
    }
}

/// Powder XRD simulator.
#[derive(Clone, Debug)]
pub struct XrdSimulator {
    /// X-ray wavelength (Å).
    pub wavelength: f64,
    /// Minimum |F|² to keep a reflection (background threshold).
    pub min_intensity: f64,
    /// Maximum number of Miller indices scanned per axis.
    pub max_index: i32,
}

impl Default for XrdSimulator {
    fn default() -> Self {
        XrdSimulator {
            wavelength: 1.5406, // Cu Kα
            min_intensity: 1e-4,
            max_index: 10,
        }
    }
}

impl XrdSimulator {
    /// Simulate the powder pattern; reflections sorted by 2θ with
    /// intensities normalized to the strongest.
    pub fn simulate(&self, lat: &Lattice, structure: &CrystalStructure) -> Vec<Reflection> {
        let recip = lat.reciprocal_no_2pi();
        let mut raw: Vec<(Miller, f64, f64, f64)> = Vec::new(); // (miller, d, 2θ, |F|²)

        for h in -self.max_index..=self.max_index {
            for k in -self.max_index..=self.max_index {
                for l in -self.max_index..=self.max_index {
                    if h == 0 && k == 0 && l == 0 {
                        continue;
                    }
                    let m = Miller::new(h, k, l);
                    let g = m.reciprocal_vector(lat);
                    let d = 1.0 / g.norm();
                    if d < self.wavelength / 2.0 {
                        continue; // beyond the Bragg limit
                    }
                    let two_theta = 2.0 * num::asin(self.wavelength / (2.0 * d)).to_degrees();

                    // Structure factor.
                    let s = num::sin(theta_from_two_theta(two_theta)) / self.wavelength;
                    let mut f_re = 0.0;
                    let mut f_im = 0.0;
                    for (z, r_frac) in structure.atomic_numbers.iter().zip(&structure.frac_positions)
                    {
                        let pos = lat.to_cartesian(*r_frac);
                        // Phase: 2π g·r with g in cycles/Å (no 2π convention).
                        let phase = 2.0 * core::f64::consts::PI * pos.dot(g);
                        let ff = form_factor(*z, s);
                        let dw = num::exp(-structure.b_factor * s * s / 2.0);
                        f_re += ff * dw * phase.cos();
                        f_im += ff * dw * phase.sin();
                    }
                    // Lorentz-polarization factor for an unpolarized
                    // source: (1 + cos^2 2θ)/(sin^2 θ cos θ).
                    let th = theta_from_two_theta(two_theta);
                    let lp = (1.0 + num::powf(2.0 * th.cos(), 2.0))
                        / (num::powf(th.sin(), 2.0) * th.cos());
                    let intensity = (f_re * f_re + f_im * f_im) * lp;
                    if intensity > self.min_intensity {
                        raw.push((m, d, two_theta, intensity));
                    }
                }
            }
        }

        if raw.is_empty() {
            return Vec::new();
        }

        // Merge symmetry-equivalent reflections (same 2θ within tolerance):
        // powder rings add intensities; the reported Miller index is the
        // family member with the largest h+k+l (prefers all-positive).
        raw.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());
        let mut merged: Vec<(Miller, f64, f64, f64)> = Vec::new();
        for (m, d, tt, inten) in raw {
            match merged.last_mut() {
                Some(last) if (last.2 - tt).abs() < 1e-4 => {
                    last.3 += inten;
                    if m.h + m.k + m.l > last.0.h + last.0.k + last.0.l {
                        last.0 = m;
                    }
                }
                _ => merged.push((m, d, tt, inten)),
            }
        }

        let max_int = merged.iter().map(|r| r.3).fold(0.0f64, f64::max);
        merged
            .into_iter()
            .map(|(m, d, tt, inten)| Reflection {
                miller: m,
                d,
                two_theta: tt,
                intensity: inten / max_int,
            })
            .collect()
    }
}

fn theta_from_two_theta(tt: f64) -> f64 {
    tt.to_radians() / 2.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn form_factor_carbon_at_zero_is_z() {
        // Cromer-Mann at s = 0 sums to Z (with the constant c).
        assert!((form_factor(6, 0.0) - 6.0).abs() < 1e-3);
        assert!((form_factor(8, 0.0) - 8.0).abs() < 1e-3);
    }

    #[test]
    fn form_factor_decays() {
        assert!(form_factor(6, 0.5) < form_factor(6, 0.1));
    }

    #[test]
    fn si_diamond_pattern() {
        // Diamond-cubic silicon: a = 5.431, atoms at (0,0,0)+fcc and
        // (¼,¼,¼)+fcc. Cu Kα.
        let lat = Lattice::cubic(5.431);
        let mut structure = CrystalStructure::new(
            vec![14; 8],
            vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(0.0, 0.5, 0.5),
                Vec3::new(0.5, 0.0, 0.5),
                Vec3::new(0.5, 0.5, 0.0),
                Vec3::new(0.25, 0.25, 0.25),
                Vec3::new(0.25, 0.75, 0.75),
                Vec3::new(0.75, 0.25, 0.75),
                Vec3::new(0.75, 0.75, 0.25),
            ],
        );
        structure.b_factor = 0.5;
        let sim = XrdSimulator {
            wavelength: 1.5406,
            min_intensity: 1e-3,
            max_index: 6,
        };
        let peaks = sim.simulate(&lat, &structure);

        // (111) ≈ 28.4° 2θ is the strongest Si peak.
        let p111 = peaks
            .iter()
            .find(|p| p.miller == Miller::new(1, 1, 1))
            .expect("111 present");
        assert!((p111.two_theta - 28.44).abs() < 0.3, "2θ = {}", p111.two_theta);

        // Extinction rules for diamond-cubic (fcc + glide):
        // (100) forbidden, (110) forbidden, (200) forbidden.
        for forbidden in [Miller::new(1, 0, 0), Miller::new(1, 1, 0), Miller::new(2, 0, 0)] {
            assert!(
                !peaks.iter().any(|p| p.miller == forbidden),
                "{forbidden:?} must be forbidden"
            );
        }

        // (220) ≈ 47.3° present and strong.
        let p220 = peaks
            .iter()
            .find(|p| p.miller == Miller::new(2, 2, 0))
            .expect("220 present");
        assert!((p220.two_theta - 47.3).abs() < 0.4, "2θ = {}", p220.two_theta);
        assert!(p220.intensity > 0.4, "220 strong: {}", p220.intensity);

        // (111) is the strongest peak.
        let strongest = peaks.iter().map(|p| p.intensity).fold(0.0f64, f64::max);
        assert!((strongest - 1.0).abs() < 1e-9);
        assert!((p111.intensity - 1.0).abs() < 1e-9);
    }

    #[test]
    fn nacl_rock_salt_extinctions() {
        // NaCl: fcc with Na at (0,0,0) and Cl at (½,0,0) (both fcc sites).
        let lat = Lattice::cubic(5.64);
        let mut structure = CrystalStructure::new(
            vec![11, 11, 11, 11, 17, 17, 17, 17],
            vec![
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(0.0, 0.5, 0.5),
                Vec3::new(0.5, 0.0, 0.5),
                Vec3::new(0.5, 0.5, 0.0),
                Vec3::new(0.5, 0.0, 0.0),
                Vec3::new(0.5, 0.5, 0.5),
                Vec3::new(0.0, 0.5, 0.0),
                Vec3::new(0.0, 0.0, 0.5),
            ],
        );
        structure.b_factor = 1.0;
        let sim = XrdSimulator {
            wavelength: 1.5406,
            min_intensity: 1e-5,
            max_index: 6,
        };
        let peaks = sim.simulate(&lat, &structure);
        // fcc: mixed-index reflections forbidden.
        for forbidden in [Miller::new(1, 0, 0), Miller::new(1, 1, 0), Miller::new(2, 1, 0)] {
            assert!(!peaks.iter().any(|p| p.miller == forbidden));
        }
        // (111) and (200) present.
        assert!(peaks.iter().any(|p| p.miller == Miller::new(1, 1, 1)));
        assert!(peaks.iter().any(|p| p.miller == Miller::new(2, 0, 0)));
    }

    #[test]
    fn peaks_sorted_by_two_theta() {
        let lat = Lattice::cubic(5.431);
        let structure = CrystalStructure::new(vec![14], vec![Vec3::ZERO]);
        let sim = XrdSimulator::default();
        let peaks = sim.simulate(&lat, &structure);
        for w in peaks.windows(2) {
            assert!(w[0].two_theta <= w[1].two_theta);
        }
    }
}
