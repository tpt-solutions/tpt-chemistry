//! From-scratch complex radix-2 FFT (iterative Cooley–Tukey) and a 3D
//! wrapper over row-column passes.
//!
//! This is the transform backend for [`crate::pme`]; it lives here rather
//! than in a `tpt-dsp` dependency so the whole long-range electrostatics
//! stack is self-contained. Sizes must be powers of two.
//!
//! Conventions (match the PME derivation):
//! * forward: `F[m] = Σ_x f[x]·e^{−2πi m x/N}`
//! * inverse: `F[x] = Σ_m f[m]·e^{+2πi m x/N}` (unnormalized; the 1/N, where
//!   needed, is applied by the caller).

/// A minimal complex number (no external numerics dependency).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Complex {
    /// Real part.
    pub re: f64,
    /// Imaginary part.
    pub im: f64,
}

impl Complex {
    /// A complex number from its parts.
    pub const fn new(re: f64, im: f64) -> Self {
        Complex { re, im }
    }
}

/// In-place forward FFT, `len` a power of two.
///
/// # Panics
/// Panics if `data.len()` is not a power of two.
pub fn fft_forward(data: &mut [Complex]) {
    let n = data.len();
    assert!(n.is_power_of_two(), "FFT length must be a power of two");
    if n <= 1 {
        return;
    }
    bit_reverse_permute(data);
    let mut len = 2;
    while len <= n {
        let half = len / 2;
        let step = -2.0 * core::f64::consts::PI / len as f64;
        for start in (0..n).step_by(len) {
            for k in 0..half {
                let w = Complex::cis(step * k as f64);
                let t = w.mul(data[start + k + half]);
                let u = data[start + k];
                data[start + k] = Complex::new(u.re + t.re, u.im + t.im);
                data[start + k + half] = Complex::new(u.re - t.re, u.im - t.im);
            }
        }
        len *= 2;
    }
}

/// In-place inverse FFT (unnormalized): `F[x] = Σ_m f[m]·e^{+2πi m x/N}`.
///
/// # Panics
/// Panics if `data.len()` is not a power of two.
pub fn fft_inverse(data: &mut [Complex]) {
    let n = data.len();
    assert!(n.is_power_of_two(), "FFT length must be a power of two");
    if n <= 1 {
        return;
    }
    bit_reverse_permute(data);
    let mut len = 2;
    while len <= n {
        let half = len / 2;
        let step = 2.0 * core::f64::consts::PI / len as f64;
        for start in (0..n).step_by(len) {
            for k in 0..half {
                let w = Complex::cis(step * k as f64);
                let t = w.mul(data[start + k + half]);
                let u = data[start + k];
                data[start + k] = Complex::new(u.re + t.re, u.im + t.im);
                data[start + k + half] = Complex::new(u.re - t.re, u.im - t.im);
            }
        }
        len *= 2;
    }
}

/// Naive DFT, for reference testing only.
pub fn dft_naive(data: &[Complex], inverse: bool) -> Vec<Complex> {
    let n = data.len();
    let sign = if inverse { 1.0 } else { -1.0 };
    (0..n)
        .map(|m| {
            let mut acc = Complex::default();
            for (x, &v) in data.iter().enumerate() {
                let ang = sign * 2.0 * core::f64::consts::PI * (m * x % n) as f64 / n as f64;
                acc = Complex::new(
                    acc.re + v.re * ang.cos() - v.im * ang.sin(),
                    acc.im + v.re * ang.sin() + v.im * ang.cos(),
                );
            }
            acc
        })
        .collect()
}

fn bit_reverse_permute(data: &mut [Complex]) {
    let n = data.len();
    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if i < j {
            data.swap(i, j);
        }
    }
}

impl Complex {
    fn cis(theta: f64) -> Complex {
        Complex::new(theta.cos(), theta.sin())
    }

    fn mul(self, o: Complex) -> Complex {
        Complex::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
}

/// A 3D complex field with power-of-two extents, stored z-fastest
/// (`index = (x·N2 + y)·N3 + z`).
#[derive(Clone)]
pub struct Field3 {
    /// Extents per axis (powers of two).
    pub dims: [usize; 3],
    /// Flat complex data.
    pub data: Vec<Complex>,
}

impl Field3 {
    /// A zero field of the given extents.
    ///
    /// # Panics
    /// Panics if any extent is not a power of two.
    pub fn zeros(dims: [usize; 3]) -> Self {
        for &d in &dims {
            assert!(d.is_power_of_two(), "field extents must be powers of two");
        }
        Field3 {
            dims,
            data: vec![Complex::default(); dims[0] * dims[1] * dims[2]],
        }
    }

    fn index(&self, x: usize, y: usize, z: usize) -> usize {
        (x * self.dims[1] + y) * self.dims[2] + z
    }

    /// Element access.
    pub fn get(&self, x: usize, y: usize, z: usize) -> Complex {
        self.data[self.index(x, y, z)]
    }

    /// Mutable element access.
    pub fn set(&mut self, x: usize, y: usize, z: usize, v: Complex) {
        let i = self.index(x, y, z);
        self.data[i] = v;
    }

    /// In-place 3D forward transform (row-column passes per axis).
    pub fn fft_forward(&mut self) {
        let [nx, ny, nz] = self.dims;
        for x in 0..nx {
            for y in 0..ny {
                let mut line: Vec<Complex> = (0..nz).map(|z| self.get(x, y, z)).collect();
                fft_forward(&mut line);
                for (z, v) in line.into_iter().enumerate() {
                    self.set(x, y, z, v);
                }
            }
        }
        for x in 0..nx {
            for z in 0..nz {
                let mut line: Vec<Complex> = (0..ny).map(|y| self.get(x, y, z)).collect();
                fft_forward(&mut line);
                for (y, v) in line.into_iter().enumerate() {
                    self.set(x, y, z, v);
                }
            }
        }
        for y in 0..ny {
            for z in 0..nz {
                let mut line: Vec<Complex> = (0..nx).map(|x| self.get(x, y, z)).collect();
                fft_forward(&mut line);
                for (x, v) in line.into_iter().enumerate() {
                    self.set(x, y, z, v);
                }
            }
        }
    }

    /// In-place 3D inverse transform (unnormalized).
    pub fn fft_inverse(&mut self) {
        let [nx, ny, nz] = self.dims;
        for x in 0..nx {
            for y in 0..ny {
                let mut line: Vec<Complex> = (0..nz).map(|z| self.get(x, y, z)).collect();
                fft_inverse(&mut line);
                for (z, v) in line.into_iter().enumerate() {
                    self.set(x, y, z, v);
                }
            }
        }
        for x in 0..nx {
            for z in 0..nz {
                let mut line: Vec<Complex> = (0..ny).map(|y| self.get(x, y, z)).collect();
                fft_inverse(&mut line);
                for (y, v) in line.into_iter().enumerate() {
                    self.set(x, y, z, v);
                }
            }
        }
        for y in 0..ny {
            for z in 0..nz {
                let mut line: Vec<Complex> = (0..nx).map(|x| self.get(x, y, z)).collect();
                fft_inverse(&mut line);
                for (x, v) in line.into_iter().enumerate() {
                    self.set(x, y, z, v);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rand_data(n: usize, seed: u64) -> Vec<Complex> {
        let mut s = seed.wrapping_mul(0x9E3779B97F4A7C15);
        let mut next = || {
            s ^= s >> 12;
            s ^= s << 25;
            s ^= s >> 27;
            let v = s.wrapping_mul(0x2545F4914F6CDD1D);
            ((v >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
        };
        (0..n).map(|_| Complex::new(next(), next())).collect()
    }

    #[test]
    fn forward_matches_dft() {
        for &n in &[2usize, 4, 8, 16, 32] {
            let data = rand_data(n, n as u64);
            let mut fast = data.clone();
            fft_forward(&mut fast);
            let slow = dft_naive(&data, false);
            for (f, s) in fast.iter().zip(slow.iter()) {
                assert!(
                    (f.re - s.re).abs() < 1e-10 && (f.im - s.im).abs() < 1e-10,
                    "n {n}: {f:?} vs {s:?}"
                );
            }
        }
    }

    #[test]
    fn inverse_matches_dft() {
        for &n in &[2usize, 8, 16] {
            let data = rand_data(n, (n * 7) as u64);
            let mut fast = data.clone();
            fft_inverse(&mut fast);
            let slow = dft_naive(&data, true);
            for (f, s) in fast.iter().zip(slow.iter()) {
                assert!(
                    (f.re - s.re).abs() < 1e-9 && (f.im - s.im).abs() < 1e-9,
                    "n {n}: {f:?} vs {s:?}"
                );
            }
        }
    }

    #[test]
    fn round_trip_is_identity() {
        let data = rand_data(64, 99);
        let mut buf = data.clone();
        fft_forward(&mut buf);
        fft_inverse(&mut buf);
        let inv_n = 1.0 / 64.0;
        for (f, d) in buf.iter().zip(data.iter()) {
            assert!((f.re * inv_n - d.re).abs() < 1e-12);
            assert!((f.im * inv_n - d.im).abs() < 1e-12);
        }
    }

    #[test]
    fn field3_matches_independent_dft() {
        let dims = [4, 4, 4];
        let mut field = Field3::zeros(dims);
        for (i, v) in field.data.iter_mut().enumerate() {
            *v = Complex::new(((i * 37) % 11) as f64 - 5.0, ((i * 13) % 7) as f64 - 3.0);
        }
        let original = field.clone();
        field.fft_forward();
        // Round trip.
        let mut back = field.clone();
        back.fft_inverse();
        let inv = 1.0 / 64.0;
        for (f, d) in back.data.iter().zip(original.data.iter()) {
            assert!((f.re * inv - d.re).abs() < 1e-12);
            assert!((f.im * inv - d.im).abs() < 1e-12);
        }
        // Independent 3D DFT at every index.
        for x in 0..4 {
            for y in 0..4 {
                for z in 0..4 {
                    let mut acc = Complex::default();
                    for i in 0..4 {
                        for j in 0..4 {
                            for k in 0..4 {
                                let v = original.data[(i * 4 + j) * 4 + k];
                                let ang = -2.0
                                    * core::f64::consts::PI
                                    * ((x * i + y * j + z * k) % 4) as f64
                                    / 4.0;
                                acc = Complex::new(
                                    acc.re + v.re * ang.cos() - v.im * ang.sin(),
                                    acc.im + v.re * ang.sin() + v.im * ang.cos(),
                                );
                            }
                        }
                    }
                    let got = field.get(x, y, z);
                    assert!(
                        (acc.re - got.re).abs() < 1e-9 && (acc.im - got.im).abs() < 1e-9,
                        "({x},{y},{z}): {acc:?} vs {got:?}"
                    );
                }
            }
        }
    }
}
