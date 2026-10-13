//! Physical constants (CODATA 2022 / SI 2019 exact values) used by the crate.

/// Molar gas constant, J·mol⁻¹·K⁻¹.
pub const R: f64 = 8.314_462_618_153_24;
/// Avogadro constant, mol⁻¹.
pub const N_A: f64 = 6.022_140_76e23;
/// Elementary charge, C.
pub const E_CHARGE: f64 = 1.602_176_634e-19;
/// Boltzmann constant, J·K⁻¹.
pub const K_B: f64 = 1.380_649e-23;
/// Vacuum permittivity, F·m⁻¹.
pub const EPS0: f64 = 8.854_187_812_8e-12;
/// Reference temperature for tabulated equilibrium constants, K (25 °C).
pub const T_REF: f64 = 298.15;
/// Molar mass of water, g·mol⁻¹.
pub const M_WATER: f64 = 18.015_28;
/// Natural logarithm of 10.
pub const LN10: f64 = std::f64::consts::LN_10;
