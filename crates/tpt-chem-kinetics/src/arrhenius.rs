//! Temperature-dependent rate constants: Arrhenius and Eyring laws
//! (spec.txt §4, "Temperature-dependent rate constants").

use tpt_chem_core::units::{
    BOLTZMANN_J_PER_K, BOLTZMANN_KJ_PER_MOL_K, GAS_CONSTANT_J_PER_MOL_K, PLANCK_J_S,
};

/// Arrhenius law: `k(T) = A·exp(−Ea/(RT))`.
///
/// `A` in the rate's native units (same units as `k`), `Ea` in
/// kJ·mol⁻¹, `T` in K. The temperature-independent fallback
/// [`RateLaw::rate_at_300k`] evaluates at `T = 298.15 K`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Arrhenius {
    /// Pre-exponential factor `A`.
    pub a: f64,
    /// Activation energy `Ea` (kJ·mol⁻¹).
    pub ea: f64,
}

impl Arrhenius {
    /// Evaluate `k(T)`.
    pub fn rate(&self, t_kelvin: f64) -> f64 {
        self.a * (-self.ea / (BOLTZMANN_KJ_PER_MOL_K * t_kelvin)).exp()
    }

    /// The constant-rate view at the reference temperature 298.15 K.
    pub fn at_300k(&self) -> f64 {
        self.rate(298.15)
    }

    /// As a generic [`RateLaw`].
    pub fn as_law(&self) -> RateLaw {
        RateLaw::Arrhenius(*self)
    }
}

/// Eyring law (transition-state theory):
/// `k(T) = (kB·T/h)·exp(ΔS‡/R)·exp(−ΔH‡/(RT))`.
///
/// `ΔH‡` in kJ·mol⁻¹, `ΔS‡` in J·mol⁻¹·K⁻¹; `k(T)` in s⁻¹ (first-order
/// units, from the `kB·T/h` prefactor).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Eyring {
    /// Enthalpy of activation `ΔH‡` (kJ·mol⁻¹).
    pub dh: f64,
    /// Entropy of activation `ΔS‡` (J·mol⁻¹·K⁻¹).
    pub ds: f64,
}

impl Eyring {
    /// Evaluate `k(T)` in s⁻¹.
    pub fn rate(&self, t_kelvin: f64) -> f64 {
        let exp_factor = (self.ds / GAS_CONSTANT_J_PER_MOL_K
            - self.dh * 1000.0 / (GAS_CONSTANT_J_PER_MOL_K * t_kelvin))
            .exp();
        (BOLTZMANN_J_PER_K * t_kelvin / PLANCK_J_S) * exp_factor
    }

    /// The constant-rate view at the reference temperature 298.15 K.
    pub fn at_300k(&self) -> f64 {
        self.rate(298.15)
    }

    /// As a generic [`RateLaw`].
    pub fn as_law(&self) -> RateLaw {
        RateLaw::Eyring(*self)
    }
}

/// Generic temperature-dependent rate law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum RateLaw {
    /// Arrhenius form.
    Arrhenius(Arrhenius),
    /// Eyring form.
    Eyring(Eyring),
}

impl RateLaw {
    /// Evaluate at a temperature.
    pub fn rate_at(&self, t_kelvin: f64) -> f64 {
        match self {
            RateLaw::Arrhenius(a) => a.rate(t_kelvin),
            RateLaw::Eyring(e) => e.rate(t_kelvin),
        }
    }

    /// The value used when the network is integrated without a temperature:
    /// the law evaluated at 298.15 K (see `ReactionInput::propensity`).
    pub fn rate_at_300k(&self) -> f64 {
        self.rate_at(298.15)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arrhenius_doubling_rule() {
        // Ea = 50 kJ/mol: k doubles roughly every ~10 K near 300 K.
        let a = Arrhenius { a: 1e13, ea: 50.0 };
        let k300 = a.rate(300.0);
        let k310 = a.rate(310.0);
        assert!(k310 / k300 > 1.5 && k310 / k300 < 2.5);
        // Exact value at 298.15 K.
        let k = a.rate(298.15);
        let expect: f64 = 1e13 * (-50.0f64 / (8.314462618e-3 * 298.15)).exp();
        assert!((k - expect).abs() < 1e-3 * expect);
    }

    #[test]
    fn arrhenius_zero_ea_is_constant() {
        let a = Arrhenius { a: 5.0, ea: 0.0 };
        assert_eq!(a.rate(100.0), 5.0);
        assert_eq!(a.rate(1000.0), 5.0);
    }

    #[test]
    fn eyring_first_order_prefactor() {
        // Diffusion-limited-ish: ΔH‡ = 10 kJ/mol, ΔS‡ = −20 J/mol/K.
        let e = Eyring {
            dh: 10.0,
            ds: -20.0,
        };
        let k298 = e.rate(298.15);
        // k = (kBT/h) exp(ΔS/R) exp(−ΔH/RT): kB T/h ≈ 6.2e12 s⁻¹ at 298 K.
        let expect: f64 = (1.380649e-23 * 298.15 / 6.62607015e-34)
            * ((-20.0f64 / 8.314462618) - 10000.0f64 / (8.314462618 * 298.15)).exp();
        assert!((k298 - expect).abs() < 1e-3 * expect);
        // k grows with T.
        assert!(e.rate(350.0) > k298);
    }

    #[test]
    fn rate_law_dispatch() {
        let law = RateLaw::Arrhenius(Arrhenius { a: 2.0, ea: 0.0 });
        assert_eq!(law.rate_at(500.0), 2.0);
        let law = RateLaw::Eyring(Eyring { dh: 0.0, ds: 0.0 });
        // kB T / h at 300 K ≈ 6.25e12.
        assert!(law.rate_at(300.0) > 6e12 && law.rate_at(300.0) < 7e12);
    }
}
