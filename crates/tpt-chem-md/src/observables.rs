//! Structural and dynamical observables from trajectories: the radial
//! distribution function g(r), mean-squared displacement, the velocity
//! autocorrelation function, and Einstein/Green–Kubo diffusion
//! coefficients.
//!
//! All functions take frame data in MD units (positions Å, velocities
//! Å·fs⁻¹); diffusion comes out in Å²·fs⁻¹ (`1 Å²·fs⁻¹ = 10⁻⁴ cm²·s⁻¹`).

use tpt_chem_core::vec3::Vec3;

use crate::box3::Box3;

/// Radial distribution function g(r) from a periodic trajectory.
///
/// `frames` are sampled configurations (wrapped positions), all with the
/// same particle count inside `box_`. Returns `r` bin centers and the
/// normalized g(r) per bin; each pair is counted once, and the shell
/// counts are normalized by the ideal-gas expectation
/// `2π N ρ r² Δr` per frame (ρ = N/V), so a uniform fluid gives g ≈ 1.
pub fn rdf(frames: &[Vec<Vec3>], box_: &Box3, r_max: f64, bins: usize) -> (Vec<f64>, Vec<f64>) {
    assert!(bins > 0, "at least one bin");
    assert!(!frames.is_empty(), "at least one frame");
    let n = frames[0].len();
    let dr = r_max / bins as f64;
    let volume = box_.volume();
    let rho = n as f64 / volume;

    let mut hist = vec![0.0f64; bins];
    for frame in frames {
        debug_assert_eq!(frame.len(), n, "all frames must have equal length");
        for a in 0..n {
            for b in (a + 1)..n {
                let d = box_.min_image(frame[b] - frame[a]);
                let r = d.norm();
                if r < r_max {
                    hist[(r / dr) as usize] += 2.0;
                }
            }
        }
    }
    let n_frames = frames.len() as f64;
    let mut centers = Vec::with_capacity(bins);
    let mut g = Vec::with_capacity(bins);
    for (i, &count) in hist.iter().enumerate() {
        let r_mid = (i as f64 + 0.5) * dr;
        // hist counts each pair in both orderings (+= 2), so the ideal
        // expectation is N·ρ·4π r² Δr per frame.
        let ideal = 4.0 * core::f64::consts::PI * n as f64 * rho * r_mid * r_mid * dr * n_frames;
        centers.push(r_mid);
        g.push(count / ideal);
    }
    (centers, g)
}

/// Mean-squared displacement as a function of lag, averaged over all time
/// origins: `MSD(k·Δt) = ⟨|r(t+k) − r(t)|²⟩_t`.
///
/// Positions must be **unwrapped** (continuous); wrapped coordinates fold
/// across the periodic boundary and corrupt large lags.
/// Returns MSD for lags `0..=max_lag` (Å²).
pub fn msd(frames: &[Vec<Vec3>], max_lag: usize) -> Vec<f64> {
    assert!(frames.len() > 1, "at least two frames");
    let n_frames = frames.len();
    let max_lag = max_lag.min(n_frames - 1);
    let mut out = Vec::with_capacity(max_lag + 1);
    for lag in 0..=max_lag {
        let mut acc = 0.0f64;
        let mut count = 0u64;
        for (t, cur) in frames.iter().take(n_frames - lag).enumerate() {
            for (i, p) in cur.iter().enumerate() {
                acc += (frames[t + lag][i] - *p).norm_sq();
            }
            count += cur.len() as u64;
        }
        out.push(if count > 0 { acc / count as f64 } else { 0.0 });
    }
    out
}

/// Velocity autocorrelation function: `C(k·Δt) = ⟨v(t)·v(t+τ)⟩_t`
/// (Å²·fs⁻²), averaged over time origins and particles.
pub fn vacf(vel_frames: &[Vec<Vec3>], max_lag: usize) -> Vec<f64> {
    assert!(vel_frames.len() > 1, "at least two frames");
    let n_frames = vel_frames.len();
    let max_lag = max_lag.min(n_frames - 1);
    let mut out = Vec::with_capacity(max_lag + 1);
    for lag in 0..=max_lag {
        let mut acc = 0.0f64;
        let mut count = 0u64;
        for (t, cur) in vel_frames.iter().take(n_frames - lag).enumerate() {
            for (i, v) in cur.iter().enumerate() {
                acc += v.dot(vel_frames[t + lag][i]);
            }
            count += cur.len() as u64;
        }
        out.push(if count > 0 { acc / count as f64 } else { 0.0 });
    }
    out
}

/// Diffusion coefficient from the Einstein relation, `D = slope(MSD)/(6)`
/// (Å²·fs⁻¹), by least squares over MSD lags `1..=fit_lag`.
///
/// The fit ignores lag 0 (always 0) and the short-lag ballistic regime
/// when `fit_lag` is chosen past it.
pub fn diffusion_from_msd(msd_curve: &[f64], dt: f64, fit_lag: usize) -> f64 {
    assert!(msd_curve.len() > 1, "need at least lag 1");
    let fit_lag = fit_lag.min(msd_curve.len() - 1).max(1);
    let n = fit_lag as f64;
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    for (lag, &y) in msd_curve.iter().enumerate().take(fit_lag + 1).skip(1) {
        let x = lag as f64 * dt;
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
    }
    let slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    slope / 6.0
}

/// Diffusion coefficient from the Green–Kubo relation,
/// `D = (1/3)∫₀^∞ C(τ) dτ` (Å²·fs⁻¹), integrating the VACF by the
/// trapezoid rule up to the series length.
pub fn diffusion_from_vacf(vacf_curve: &[f64], dt: f64) -> f64 {
    if vacf_curve.len() < 2 {
        return 0.0;
    }
    let mut integral = 0.0;
    for pair in vacf_curve.windows(2) {
        integral += (pair[0] + pair[1]) * 0.5 * dt;
    }
    integral / 3.0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn uniform(box_: &Box3, n: usize, seed: u64) -> Vec<Vec3> {
        let mut rng = tpt_chem_core::rng::Rng::new(seed);
        (0..n)
            .map(|_| {
                Vec3::new(
                    rng.range(0.0, box_.length.x),
                    rng.range(0.0, box_.length.y),
                    rng.range(0.0, box_.length.z),
                )
            })
            .collect()
    }

    #[test]
    fn rdf_of_uniform_gas_is_one() {
        let box_ = Box3::cubic(20.0);
        let frames: Vec<Vec<Vec3>> = (0..40).map(|s| uniform(&box_, 300, s)).collect();
        let (centers, g) = rdf(&frames, &box_, 5.0, 25);
        // Interior bins (well inside r_max, away from boundary artifacts):
        // g ≈ 1 within sampling noise.
        for (r, gi) in centers.iter().skip(4).zip(g.iter().skip(4)) {
            assert!(
                (gi - 1.0).abs() < 0.15,
                "g(r = {r}) = {gi} too far from unity"
            );
        }
    }

    #[test]
    fn rdf_of_lattice_peaks_at_shell_distances() {
        // Perfect simple-cubic lattice: g(r) is a comb of shells at
        // a, √2 a, √3 a, ... — with finite bins, peaks land in the right
        // bins and the empty gap bins stay near zero.
        let a = 3.0;
        let box_ = Box3::cubic(4.0 * a);
        let mut frame = Vec::new();
        for i in 1..4 {
            for j in 1..4 {
                for k in 1..4 {
                    frame.push(Vec3::new(i as f64 * a, j as f64 * a, k as f64 * a));
                }
            }
        }
        let (centers, g) = rdf(&[frame], &box_, 6.0, 60);
        let peak_at = |target: f64| -> f64 {
            let bin = (target / 0.1) as usize;
            let lo = bin.saturating_sub(1);
            g[lo..(bin + 2).min(g.len())]
                .iter()
                .fold(0.0f64, |m, v| m.max(*v))
        };
        // Shell normalization for a 3×3×3 interior block is smaller than
        // the ideal-gas estimate (finite-N), so only peak *positions* are
        // asserted here.
        let s1 = peak_at(a);
        let s2 = peak_at(a * core::f64::consts::SQRT_2);
        let gap = g[15];
        assert!(s1 > 1.0, "1st shell peak missing: {s1}");
        assert!(s2 > 0.5, "2nd shell peak missing: {s2}");
        assert!(gap < 0.05, "gap bin must be empty: {gap}");
        let _ = centers;
    }

    #[test]
    fn msd_of_ballistic_motion_is_v2t2() {
        let v = Vec3::new(0.3, -0.1, 0.2);
        let frames: Vec<Vec<Vec3>> = (0..50)
            .map(|t| vec![v * (t as f64), v * (t as f64) + Vec3::new(5.0, 0.0, 0.0)])
            .collect();
        let curve = msd(&frames, 40);
        for (lag, m) in curve.iter().enumerate() {
            let want = v.norm_sq() * (lag * lag) as f64;
            assert!((m - want).abs() < 1e-9, "lag {lag}: {m} vs {want}");
        }
    }

    #[test]
    fn diffusion_recovers_known_d_from_brownian_motion() {
        // Brownian increments with D = 0.05 Å²/fs: each frame step is a
        // Gaussian with variance 6 D dt (3D).
        let d_true: f64 = 0.05;
        let dt = 1.0;
        let n_frames = 2_000;
        let n_parts = 50;
        let mut rng = tpt_chem_core::rng::Rng::new(99);
        let mut frames: Vec<Vec<Vec3>> = Vec::new();
        let sigma = (2.0 * d_true * dt).sqrt();
        let mut cur: Vec<Vec3> = (0..n_parts).map(|_| Vec3::ZERO).collect();
        for _ in 0..n_frames {
            frames.push(cur.clone());
            for p in cur.iter_mut() {
                *p += Vec3::new(
                    rng.normal(0.0, sigma),
                    rng.normal(0.0, sigma),
                    rng.normal(0.0, sigma),
                );
            }
        }
        let curve = msd(&frames, 400);
        let d = diffusion_from_msd(&curve, dt, 300);
        assert!((d - d_true).abs() / d_true < 0.08, "D = {d} vs {d_true}");
    }

    #[test]
    fn vacf_of_constant_velocity_is_flat_and_kubo_matches() {
        let v = Vec3::new(0.2, 0.1, -0.15);
        let frames: Vec<Vec<Vec3>> = (0..100).map(|_| vec![v; 3]).collect();
        let c = vacf(&frames, 50);
        for (lag, val) in c.iter().enumerate() {
            assert!(
                (val - v.norm_sq()).abs() < 1e-12,
                "lag {lag}: {val} vs {}",
                v.norm_sq()
            );
        }
        // D = ∫C/3 over the flat curve: 50 fs of it.
        let d = diffusion_from_vacf(&c, 1.0);
        let want = v.norm_sq() * 50.0 / 3.0; // trapezoid over 50 windows
        assert!((d - want).abs() < 1e-9);
    }
}
