//! Brillouin-zone sampling: Monkhorst-Pack and Γ-centered k-point grids
//! in fractional reciprocal coordinates.

use std::vec::Vec;

/// A k-point in fractional reciprocal coordinates plus its weight.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KPoint {
    /// Fractional reciprocal coordinates.
    pub frac: [f64; 3],
    /// Statistical weight (sums to 1 over the grid).
    pub weight: f64,
}

/// A generated k-point grid.
#[derive(Clone, Debug)]
pub struct KGrid {
    /// The points with weights summing to 1.
    pub points: Vec<KPoint>,
    /// Grid subdivisions along each reciprocal axis.
    pub divisions: [usize; 3],
}

impl KGrid {
    /// Monkhorst-Pack grid: `k_i = (2r_i − N_i + 1)/(2N_i)` for
    /// r = 1..N (uniform, includes Γ only when N is odd).
    pub fn monkhorst_pack(divisions: [usize; 3]) -> Self {
        let [na, nb, nc] = divisions;
        let total = (na * nb * nc) as f64;
        let mut points = Vec::with_capacity(na * nb * nc);
        for i in 0..na {
            for j in 0..nb {
                for k in 0..nc {
                    points.push(KPoint {
                        frac: [
                            (2 * i as i32 - na as i32 + 1) as f64 / (2.0 * na as f64),
                            (2 * j as i32 - nb as i32 + 1) as f64 / (2.0 * nb as f64),
                            (2 * k as i32 - nc as i32 + 1) as f64 / (2.0 * nc as f64),
                        ],
                        weight: 1.0 / total,
                    });
                }
            }
        }
        KGrid { points, divisions }
    }

    /// Γ-centered grid: k_i = r_i/N_i for r = 0..N (always includes Γ);
    /// the standard choice for supercell-style Brillouin sampling.
    pub fn gamma_centered(divisions: [usize; 3]) -> Self {
        let [na, nb, nc] = divisions;
        let total = (na * nb * nc) as f64;
        let mut points = Vec::with_capacity(na * nb * nc);
        for i in 0..na {
            for j in 0..nb {
                for k in 0..nc {
                    points.push(KPoint {
                        frac: [
                            i as f64 / na as f64,
                            j as f64 / nb as f64,
                            k as f64 / nc as f64,
                        ],
                        weight: 1.0 / total,
                    });
                }
            }
        }
        KGrid { points, divisions }
    }

    /// Total weight (1.0 up to floating point).
    pub fn total_weight(&self) -> f64 {
        self.points.iter().map(|p| p.weight).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn weights_sum_to_one() {
        let g = KGrid::monkhorst_pack([2, 2, 2]);
        assert!((g.total_weight() - 1.0).abs() < 1e-12);
        assert_eq!(g.points.len(), 8);
        let g = KGrid::gamma_centered([3, 3, 3]);
        assert!((g.total_weight() - 1.0).abs() < 1e-12);
        assert_eq!(g.points.len(), 27);
    }

    #[test]
    fn gamma_centered_includes_gamma() {
        let g = KGrid::gamma_centered([2, 2, 2]);
        assert!(g.points.iter().any(|p| p.frac == [0.0; 3]));
    }

    #[test]
    fn monkhorst_pack_2x_shifts_off_gamma() {
        let g = KGrid::monkhorst_pack([2, 2, 2]);
        // 2×2×2 MP points sit at (±1/4, ±1/4, ±1/4).
        assert!(g.points.iter().all(|p| p
            .frac
            .iter()
            .all(|&v| (v.abs() - 0.25).abs() < 1e-12)));
    }

    #[test]
    fn gamma_centered_values() {
        let g = KGrid::gamma_centered([2, 1, 1]);
        let fracs: Vec<[f64; 3]> = g.points.iter().map(|p| p.frac).collect();
        assert!(fracs.contains(&[0.0, 0.0, 0.0]));
        assert!(fracs.contains(&[0.5, 0.0, 0.0]));
        assert_eq!(g.points.len(), 2);
    }
}
