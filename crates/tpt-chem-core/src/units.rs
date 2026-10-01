//! Physical constants (NIST CODATA 2022) and chemistry unit conversions.
//!
//! Two unit systems are used across `tpt-chemistry`:
//!
//! * **MD units**: Å, fs, amu, kJ·mol⁻¹ — the default for classical
//!   simulation.
//! * **Atomic units** (a.u.): Bohr, ħ/Eₕ, mₑ, Hartree — the default for
//!   electronic structure.
//!
//! The conversion factors below are `const` so they fold into compile-time
//! expressions; the `uom` re-export (via `tpt-math-units`) is available for
//! fully unit-checked SI APIs at library boundaries.
//!
//! ```
//! use tpt_chem_core::units::*;
//!
//! let r_bohr = 2.0;
//! assert!((r_bohr * BOHR_ANGSTROM - 1.05835442).abs() < 1e-7); // Å
//! assert!((HARTREE_KJ_PER_MOL - 2625.4996).abs() < 1e-3);
//! ```

/// Re-export of the `uom`-backed SI quantity system (`tpt-math-units`) for
/// callers that want fully unit-checked quantities.
pub use tpt_math_units;

/// Bohr radius `a₀` in Ångström (CODATA 2022).
pub const BOHR_ANGSTROM: f64 = 0.529_177_210_544;

/// Hartree energy `Eₕ` in kJ·mol⁻¹.
pub const HARTREE_KJ_PER_MOL: f64 = 2_625.499_639_481;

/// Hartree energy `Eₕ` in kcal·mol⁻¹.
pub const HARTREE_KCAL_PER_MOL: f64 = 627.509_474_063_1;

/// Hartree energy `Eₕ` in eV.
pub const HARTREE_EV: f64 = 27.211_386_245_981;

/// Hartree energy `Eₕ` in joules.
pub const HARTREE_J: f64 = 4.359_744_722_206_1e-18;

/// Bohr radius `a₀` in metres.
pub const BOHR_M: f64 = 5.291_772_105_44e-11;

/// Atomic unit of time `ħ/Eₕ` in seconds.
pub const ATOMIC_TIME_S: f64 = 2.418_884_326_586_4e-17;

/// Atomic unit of time `ħ/Eₕ` in femtoseconds.
pub const ATOMIC_TIME_FS: f64 = 0.024_188_843_265_864;

/// Atomic unit of velocity `a₀·Eₕ/ħ` in m·s⁻¹.
pub const ATOMIC_VELOCITY_M_PER_S: f64 = 2.187_691_262_16e6;

/// Electron mass in amu.
pub const ELECTRON_MASS_AMU: f64 = 5.485_799_090_441e-4;

/// Atomic mass unit in kg.
pub const AMU_KG: f64 = 1.660_539_068_92e-27;

/// Avogadro constant `N_A` in mol⁻¹ (exact, SI 2019).
pub const AVOGADRO: f64 = 6.022_140_76e23;

/// Boltzmann constant `k_B` in J·K⁻¹ (exact, SI 2019).
pub const BOLTZMANN_J_PER_K: f64 = 1.380_649e-23;

/// Molar gas constant `R = k_B·N_A` in J·mol⁻¹·K⁻¹.
pub const GAS_CONSTANT_J_PER_MOL_K: f64 = 8.314_462_618;

/// Boltzmann constant `k_B` in kJ·mol⁻¹·K⁻¹ — the form used by MD code.
pub const BOLTZMANN_KJ_PER_MOL_K: f64 = 8.314_462_618_153_24e-3;

/// Boltzmann constant `k_B` in Hartree·K⁻¹ — the form used by quantum code.
pub const BOLTZMANN_HARTREE_PER_K: f64 = 3.166_811_563_452_4e-6;

/// Planck constant `h` in J·s (exact, SI 2019).
pub const PLANCK_J_S: f64 = 6.626_070_15e-34;

/// Reduced Planck constant `ħ` in J·s (exact, SI 2019).
pub const HBAR_J_S: f64 = 1.054_571_817e-34;

/// Speed of light `c` in m·s⁻¹ (exact, SI 2019).
pub const C_M_PER_S: f64 = 299_792_458.0;

/// Elementary charge in coulombs (exact, SI 2019).
pub const ELEMENTARY_CHARGE_C: f64 = 1.602_176_634e-19;

/// Vacuum permittivity `ε₀` in F·m⁻¹ (CODATA 2022).
pub const VACUUM_PERMITTIVITY: f64 = 8.854_187_818_8e-12;

/// Coulomb constant `1/(4πε₀)` in N·m²·C⁻² (CODATA 2022).
pub const COULOMB_CONSTANT: f64 = 8.987_551_791_9e9;

/// Coulomb's-law prefactor in kJ·mol⁻¹·Å per (e·e): the constant in
/// `E = COULOMB_PREFACTOR_KJ_ANG * qa * qb / r_Å` for charges in units of
/// the elementary charge.
pub const COULOMB_PREFACTOR_KJ_ANG: f64 = 1_389.354_582_655;

/// 1 kcal in kJ (thermochemical calorie, exact: 1 cal = 4.184 J).
pub const KCAL_PER_KJ: f64 = 1.0 / 4.184;

/// Convert a length from Bohr to Ångström.
pub fn bohr_to_angstrom(r_bohr: f64) -> f64 {
    r_bohr * BOHR_ANGSTROM
}

/// Convert a length from Ångström to Bohr.
pub fn angstrom_to_bohr(r_ang: f64) -> f64 {
    r_ang / BOHR_ANGSTROM
}

/// Convert an energy from Hartree to kJ·mol⁻¹.
pub fn hartree_to_kj_per_mol(e_hartree: f64) -> f64 {
    e_hartree * HARTREE_KJ_PER_MOL
}

/// Convert an energy from kJ·mol⁻¹ to Hartree.
pub fn kj_per_mol_to_hartree(e_kj_mol: f64) -> f64 {
    e_kj_mol / HARTREE_KJ_PER_MOL
}

/// Convert an energy from Hartree to eV.
pub fn hartree_to_ev(e_hartree: f64) -> f64 {
    e_hartree * HARTREE_EV
}

/// Convert an energy from kJ·mol⁻¹ to kcal·mol⁻¹.
pub fn kj_per_mol_to_kcal_per_mol(e_kj_mol: f64) -> f64 {
    e_kj_mol * KCAL_PER_KJ
}

/// Convert a time from atomic units to femtoseconds.
pub fn atomic_time_to_fs(t_au: f64) -> f64 {
    t_au * ATOMIC_TIME_FS
}

/// Convert a temperature to thermal energy in kJ·mol⁻¹ (`k_B·T`).
pub fn thermal_energy_kj_per_mol(t_kelvin: f64) -> f64 {
    t_kelvin * BOLTZMANN_KJ_PER_MOL_K
}

/// Convert a temperature to thermal energy in Hartree (`k_B·T`).
pub fn thermal_energy_hartree(t_kelvin: f64) -> f64 {
    t_kelvin * BOLTZMANN_HARTREE_PER_K
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() <= tol
    }

    #[test]
    fn roundtrips() {
        for r in [0.5, 1.0, 2.0, 10.0, 0.74] {
            assert!(close(bohr_to_angstrom(angstrom_to_bohr(r)), r, 1e-13));
        }
        for e in [0.1, 1.0, -74.96, 27.2] {
            assert!(close(
                hartree_to_kj_per_mol(kj_per_mol_to_hartree(e)),
                e,
                1e-10
            ));
        }
    }

    #[test]
    fn reference_conversions() {
        // H₂ equilibrium bond length 0.7414 Å in Bohr.
        assert!(close(angstrom_to_bohr(0.7414), 1.401_035, 1e-5));
        // 1 Hartree ≈ 627.5 kcal/mol.
        assert!(close(
            hartree_to_kj_per_mol(1.0) * KCAL_PER_KJ,
            627.5095,
            1e-3
        ));
        // 1 eV in kJ/mol: 96.485.
        assert!(close(HARTREE_KJ_PER_MOL / HARTREE_EV, 96.485332, 1e-5));
    }

    #[test]
    fn thermal_energies() {
        // k_B·T at 298.15 K ≈ 2.479 kJ/mol.
        assert!(close(thermal_energy_kj_per_mol(298.15), 2.4790, 1e-3));
        // k_B·T at 298.15 K ≈ 9.4410e-4 Hartree.
        assert!(close(thermal_energy_hartree(298.15), 9.441_0e-4, 1e-7));
    }
}
