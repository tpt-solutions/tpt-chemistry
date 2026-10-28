//! Gaussian `.cube` output of RHF densities and molecular orbitals.
//!
//! Values are evaluated on an axis-aligned grid of points in Bohr, from the
//! same Cartesian Gaussian shells the SCF used. The cube writer converts
//! atom positions to Bohr, as the format requires, and writes values with z
//! varying fastest.
//!
//! The grid is a plain point set, so the caller can check any integral of a
//! value over it. For example, the RHF density integrates to the electron
//! count, and each normalized MO integrates to one, once the grid spans the
//! functions and the spacing is fine enough.
//!
//! **Spacing must resolve the tightest Gaussian.** The oxygen 1s STO-3G
//! primitives have exponents near 130 Bohr⁻². On water, a 0.25 Bohr grid
//! gives ∫ρ = 10.70 electrons, 0.15 Bohr gives 9.96 and 0.10 Bohr gives
//! 10.005. Use about 0.1 Bohr when the heavy-atom core matters.

use tpt_chem_core::element::symbol;
use tpt_chem_core::molecule::Molecule;
use tpt_chem_core::units::angstrom_to_bohr;
use tpt_chem_core::vec3::Vec3;

use crate::gaussian::{cartesian_components, n_cartesian, primitive_normalization, Shell};

/// Axis-aligned grid of points, in Bohr.
#[derive(Debug, Clone, PartialEq)]
pub struct CubeGrid {
    /// Position of the first grid point (Bohr).
    pub origin: Vec3,
    /// Spacing along x, y and z (Bohr).
    pub spacing: [f64; 3],
    /// Number of points along x, y and z.
    pub n: [usize; 3],
}

impl CubeGrid {
    /// A grid enclosing `points` (Bohr) with `margin` Bohr of padding on every
    /// side and roughly `spacing` Bohr between points.
    ///
    /// # Panics
    /// Panics if `points` is empty or `spacing` is not positive.
    pub fn enclosing(points: &[Vec3], margin: f64, spacing: f64) -> Self {
        assert!(!points.is_empty(), "grid needs at least one point");
        assert!(spacing > 0.0, "spacing must be positive");
        let mut lo = points[0];
        let mut hi = points[0];
        for p in points {
            lo = Vec3::new(lo.x.min(p.x), lo.y.min(p.y), lo.z.min(p.z));
            hi = Vec3::new(hi.x.max(p.x), hi.y.max(p.y), hi.z.max(p.z));
        }
        let lo = lo - Vec3::new(margin, margin, margin);
        let span = |a: f64, b: f64| ((b - a + 2.0 * margin) / spacing).ceil() as usize + 1;
        Self {
            origin: lo,
            spacing: [spacing; 3],
            n: [span(lo.x, hi.x), span(lo.y, hi.y), span(lo.z, hi.z)],
        }
    }

    /// Total number of grid points.
    pub fn len(&self) -> usize {
        self.n[0] * self.n[1] * self.n[2]
    }

    /// `true` if the grid has no points.
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Volume element `dx·dy·dz` (Bohr³), for integrating grid values.
    pub fn voxel_volume(&self) -> f64 {
        self.spacing[0] * self.spacing[1] * self.spacing[2]
    }

    /// Position of grid point `(i, j, k)` (Bohr).
    pub fn point(&self, i: usize, j: usize, k: usize) -> Vec3 {
        self.origin
            + Vec3::new(
                i as f64 * self.spacing[0],
                j as f64 * self.spacing[1],
                k as f64 * self.spacing[2],
            )
    }

    /// Iterates over every grid point in cube-file order (z fastest).
    pub fn points(&self) -> impl Iterator<Item = Vec3> + '_ {
        (0..self.n[0]).flat_map(move |i| {
            (0..self.n[1]).flat_map(move |j| (0..self.n[2]).map(move |k| self.point(i, j, k)))
        })
    }
}

/// Values of every contracted Cartesian Gaussian at `r` (Bohr), in the
/// basis-function order the SCF uses (shell by shell, components in
/// [`cartesian_components`] order).
pub fn ao_values(shells: &[Shell], r: Vec3) -> Vec<f64> {
    let mut out = Vec::with_capacity(shells.iter().map(|s| n_cartesian(s.l)).sum());
    for sh in shells {
        let d = r - sh.center;
        let r2 = d.norm_sq();
        let radial: Vec<(f64, f64)> = sh
            .primitives
            .iter()
            .map(|&(alpha, c)| (c, (-alpha * r2).exp()))
            .collect();
        for (ic, &comp) in cartesian_components(sh.l).iter().enumerate() {
            let poly = d.x.powi(i32::from(comp.0))
                * d.y.powi(i32::from(comp.1))
                * d.z.powi(i32::from(comp.2));
            let mut value = 0.0;
            for (p, &(c, e)) in sh.primitives.iter().zip(&radial) {
                let norm = sh.normalizations[ic] * primitive_normalization(p.0, comp);
                value += c * norm * e;
            }
            out.push(poly * value);
        }
    }
    out
}

/// Electron density `ρ(r) = Σ_ij D_ij φ_i(r) φ_j(r)` at every grid point.
///
/// `density` is the `nb × nb` AO density matrix from [`crate::hf::HfResult`]
/// (row-major, `nb = Σ n_cartesian(l)`).
pub fn density_on_grid(shells: &[Shell], density: &[f64], nb: usize, grid: &CubeGrid) -> Vec<f64> {
    assert_eq!(density.len(), nb * nb, "density must be nb × nb");
    grid.points()
        .map(|r| {
            let phi = ao_values(shells, r);
            let mut rho = 0.0;
            for (i, pi) in phi.iter().enumerate() {
                for (j, pj) in phi.iter().enumerate() {
                    rho += density[i * nb + j] * pi * pj;
                }
            }
            rho
        })
        .collect()
}

/// Value of molecular orbital `mo` at every grid point.
///
/// `coefficients` is the `nb × nb` AO-by-MO matrix from
/// [`crate::hf::HfResult`] (row-major: `coefficients[ao * nb + mo]`).
pub fn orbital_on_grid(
    shells: &[Shell],
    coefficients: &[f64],
    nb: usize,
    mo: usize,
    grid: &CubeGrid,
) -> Vec<f64> {
    assert_eq!(coefficients.len(), nb * nb, "coefficients must be nb × nb");
    assert!(mo < nb, "MO index out of range");
    grid.points()
        .map(|r| {
            ao_values(shells, r)
                .iter()
                .enumerate()
                .map(|(ao, phi)| coefficients[ao * nb + mo] * phi)
                .sum()
        })
        .collect()
}

/// Renders scalar values on `grid` as a Gaussian `.cube` file.
///
/// `mol` carries Å positions; they are converted to Bohr. `values` must have
/// one entry per grid point in [`CubeGrid::points`] order.
///
/// # Panics
/// Panics if `values.len() != grid.len()`, or if an atom has no symbol.
pub fn write_cube(mol: &Molecule, grid: &CubeGrid, values: &[f64], comment: &str) -> String {
    assert_eq!(values.len(), grid.len(), "one value per grid point");
    let to_bohr = angstrom_to_bohr(1.0);
    let mut out = String::new();
    out.push_str(comment.lines().next().unwrap_or(""));
    out.push('\n');
    out.push_str(comment.lines().nth(1).unwrap_or(""));
    out.push('\n');
    out.push_str(&format!(
        "{:>5}{:>12.6}{:>12.6}{:>12.6}\n",
        mol.len(),
        grid.origin.x,
        grid.origin.y,
        grid.origin.z
    ));
    let axes = [
        Vec3::new(grid.spacing[0], 0.0, 0.0),
        Vec3::new(0.0, grid.spacing[1], 0.0),
        Vec3::new(0.0, 0.0, grid.spacing[2]),
    ];
    for (n, a) in grid.n.iter().zip(axes) {
        out.push_str(&format!(
            "{:>5}{:>12.6}{:>12.6}{:>12.6}\n",
            n, a.x, a.y, a.z
        ));
    }
    for atom in mol.atoms() {
        assert!(symbol(atom.z).is_some(), "atom {} has no symbol", atom.z);
        let p = atom.pos * to_bohr;
        out.push_str(&format!(
            "{:>5}{:>12.6}{:>12.6}{:>12.6}{:>12.6}\n",
            atom.z,
            f64::from(atom.z),
            p.x,
            p.y,
            p.z
        ));
    }
    for (idx, v) in values.iter().enumerate() {
        out.push_str(&format!("{v:13.5E}"));
        if (idx + 1) % 6 == 0 || idx + 1 == values.len() {
            out.push('\n');
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::sto3g_basis;
    use crate::hf::rhf;

    fn h2(r_ang: f64) -> Molecule {
        let mut m = Molecule::new("H2");
        m.add_atom_raw(1, Vec3::new(0.0, 0.0, 0.0));
        m.add_atom_raw(1, Vec3::new(0.0, 0.0, r_ang));
        m
    }

    fn shells_for(mol: &Molecule) -> Vec<Shell> {
        let atoms: Vec<(u8, Vec3)> = mol
            .atoms()
            .map(|a| (a.z, a.pos * angstrom_to_bohr(1.0)))
            .collect();
        sto3g_basis(&atoms).unwrap()
    }

    #[test]
    fn density_integrates_to_electron_count() {
        let mol = h2(0.7414);
        let res = rhf(&mol).unwrap();
        let shells = shells_for(&mol);
        let bohr_pts: Vec<Vec3> = mol.atoms().map(|a| a.pos * angstrom_to_bohr(1.0)).collect();
        let grid = CubeGrid::enclosing(&bohr_pts, 7.0, 0.3);
        let rho = density_on_grid(&shells, &res.density, res.n_basis, &grid);
        let total: f64 = rho.iter().sum::<f64>() * grid.voxel_volume();
        assert!((total - 2.0).abs() < 2e-2, "∫ρ = {total}");
    }

    #[test]
    fn occupied_orbital_is_normalized() {
        let mol = h2(0.7414);
        let res = rhf(&mol).unwrap();
        let shells = shells_for(&mol);
        let bohr_pts: Vec<Vec3> = mol.atoms().map(|a| a.pos * angstrom_to_bohr(1.0)).collect();
        let grid = CubeGrid::enclosing(&bohr_pts, 7.0, 0.3);
        let psi = orbital_on_grid(&shells, &res.coefficients, res.n_basis, 0, &grid);
        let norm: f64 = psi.iter().map(|v| v * v).sum::<f64>() * grid.voxel_volume();
        assert!((norm - 1.0).abs() < 2e-2, "∫ψ² = {norm}");
    }

    #[test]
    fn grid_points_follow_cube_order() {
        let grid = CubeGrid {
            origin: Vec3::new(0.0, 0.0, 0.0),
            spacing: [1.0, 2.0, 3.0],
            n: [2, 2, 2],
        };
        let pts: Vec<Vec3> = grid.points().collect();
        assert_eq!(pts.len(), 8);
        // z varies fastest.
        assert_eq!(pts[1], Vec3::new(0.0, 0.0, 3.0));
        assert_eq!(pts[2], Vec3::new(0.0, 2.0, 0.0));
    }

    #[test]
    fn cube_file_has_header_atoms_and_values() {
        let mol = h2(0.74);
        let grid = CubeGrid {
            origin: Vec3::new(-1.0, -1.0, -1.0),
            spacing: [0.5, 0.5, 0.5],
            n: [2, 2, 2],
        };
        let text = write_cube(&mol, &grid, &[0.0; 8], "title\nsubtitle");
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "title");
        assert!(lines[2].trim_start().starts_with('2'), "{}", lines[2]);
        // 2 comment + 1 header + 3 axis + 2 atoms + 8 values (2 lines of 6 and 2)
        assert_eq!(lines.len(), 2 + 1 + 3 + 2 + 2);
    }
}
