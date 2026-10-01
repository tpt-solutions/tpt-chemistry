//! Float math shims so the crate compiles with and without `std`.
//!
//! With the `std` feature the inherent `f64` methods are used; without it the
//! calls route through `libm` (MIT OR Apache-2.0), keeping the crate
//! `no_std + alloc` compatible per spec.txt §4.

/// √x.
pub fn sqrt(x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.sqrt()
    }
    #[cfg(not(feature = "std"))]
    {
        libm::sqrt(x)
    }
}

/// eˣ.
pub fn exp(x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.exp()
    }
    #[cfg(not(feature = "std"))]
    {
        libm::exp(x)
    }
}

/// ln x.
pub fn ln(x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.ln()
    }
    #[cfg(not(feature = "std"))]
    {
        libm::log(x)
    }
}

/// log₁₀ x.
pub fn log10(x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.log10()
    }
    #[cfg(not(feature = "std"))]
    {
        libm::log10(x)
    }
}

/// xʸ.
pub fn powf(x: f64, y: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.powf(y)
    }
    #[cfg(not(feature = "std"))]
    {
        libm::pow(x, y)
    }
}

/// xⁿ for integer n.
pub fn powi(x: f64, n: i32) -> f64 {
    #[cfg(feature = "std")]
    {
        x.powi(n)
    }
    #[cfg(not(feature = "std"))]
    {
        // Exponentiation by squaring; exact for the small |n| used here.
        match n {
            0 => 1.0,
            n if n < 0 => 1.0 / powi(x, -n),
            _ => {
                let (mut base, mut acc, mut e) = (x, 1.0, n);
                while e > 0 {
                    if e & 1 == 1 {
                        acc *= base;
                    }
                    base *= base;
                    e >>= 1;
                }
                acc
            }
        }
    }
}

/// cos x.
pub fn cos(x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.cos()
    }
    #[cfg(not(feature = "std"))]
    {
        libm::cos(x)
    }
}

/// sin x.
pub fn sin(x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.sin()
    }
    #[cfg(not(feature = "std"))]
    {
        libm::sin(x)
    }
}

/// tan x.
pub fn tan(x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.tan()
    }
    #[cfg(not(feature = "std"))]
    {
        libm::tan(x)
    }
}

/// atan2(y, x) — the angle of the point (x, y).
pub fn atan2(y: f64, x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        y.atan2(x)
    }
    #[cfg(not(feature = "std"))]
    {
        libm::atan2(y, x)
    }
}

/// arccos x, clamped to the valid domain.
pub fn acos(x: f64) -> f64 {
    #[cfg(feature = "std")]
    {
        x.clamp(-1.0, 1.0).acos()
    }
    #[cfg(not(feature = "std"))]
    {
        libm::acos(x.clamp(-1.0, 1.0))
    }
}

/// Fused multiply-add `x * a + b`.
pub fn fma(x: f64, a: f64, b: f64) -> f64 {
    x * a + b
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basics() {
        assert!((sqrt(4.0) - 2.0).abs() < 1e-15);
        assert!((exp(0.0) - 1.0).abs() < 1e-15);
        assert!((ln(1.0)).abs() < 1e-15);
        assert!((powf(2.0, 10.0) - 1024.0).abs() < 1e-9);
        assert!((powi(3.0, 4) - 81.0).abs() < 1e-12);
        assert!((cos(0.0) - 1.0).abs() < 1e-15);
        assert!((sin(core::f64::consts::FRAC_PI_2) - 1.0).abs() < 1e-15);
        assert!((atan2(1.0, 1.0) - core::f64::consts::FRAC_PI_4).abs() < 1e-15);
        assert!((acos(1.0)).abs() < 1e-15);
        assert!((fma(3.0, 2.0, 1.0) - 7.0).abs() < 1e-15);
    }
}
