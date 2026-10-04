//! Wigner-Seitz cell construction: the first Brillouin zone of the
//! reciprocal lattice.
//!
//! The Wigner-Seitz cell of a lattice is the set of points closer to the
//! origin than to any other lattice point. For the reciprocal lattice it
//! is the first Brillouin zone. Here it is computed as the intersection of
//! half-spaces `x·g ≤ g·g/2` for all nearest reciprocal lattice vectors
//! `g`, represented by its bounding planes plus the origin-containing
//! polyhedron vertices obtained by plane enumeration.

use std::vec::Vec;

use tpt_chem_core::num;
use tpt_chem_core::vec3::Vec3;

use crate::bravais::Lattice;

/// A bounding half-space: points `x` with `normal·x ≤ offset`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HalfSpace {
    /// Outward normal.
    pub normal: Vec3,
    /// Offset: points satisfy `normal·x ≤ offset`.
    pub offset: f64,
}

impl HalfSpace {
    /// The perpendicular bisector half-space of the lattice point `g`:
    /// points closer to the origin than to `g`.
    pub fn bisector(g: Vec3) -> Self {
        HalfSpace {
            normal: g,
            offset: g.norm_sq() / 2.0,
        }
    }

    /// Distance of a point beyond the plane (positive = outside).
    pub fn outside_by(&self, x: Vec3) -> f64 {
        self.normal.dot(x) - self.offset
    }
}

/// A Wigner-Seitz cell: bounding half-spaces (the Brillonin-zone faces)
/// plus representative vertices where ≥ 3 planes meet.
#[derive(Clone, Debug)]
pub struct WignerSeitzCell {
    /// Bounding half-spaces (face planes).
    pub faces: Vec<HalfSpace>,
    /// Vertices of the cell (in 1/Å⁻¹-reciprocal coordinates with the
    /// crystallographic 2π convention — same units as `Lattice::reciprocal`).
    pub vertices: Vec<Vec3>,
}

impl WignerSeitzCell {
    /// Wigner-Seitz cell of the reciprocal lattice of `lat` (the first
    /// Brillouin zone). Considers all reciprocal lattice points in a
    /// ±2-supercell neighborhood of the origin.
    pub fn brillouin_zone(lat: &Lattice) -> Self {
        let recip = lat.reciprocal();
        let mut points = Vec::new();
        for i in -2i32..=2 {
            for j in -2i32..=2 {
                for k in -2i32..=2 {
                    if i == 0 && j == 0 && k == 0 {
                        continue;
                    }
                    points.push(recip.a * i as f64 + recip.b * j as f64 + recip.c * k as f64);
                }
            }
        }
        Self::of_lattice_points(&points)
    }

    /// Wigner-Seitz cell for an explicit list of lattice points (the
    /// origin is implicit and excluded from the input).
    ///
    /// Two-pass construction: enumerate plane-triple intersections against
    /// all candidate half-spaces, keep faces with ≥3 vertices, then
    /// recompute the vertices from the surviving faces only.
    pub fn of_lattice_points(lattice_points: &[Vec3]) -> Self {
        let mut candidate_faces: Vec<HalfSpace> = Vec::new();
        for &g in lattice_points {
            if g.norm_sq() < 1e-18 {
                continue;
            }
            candidate_faces.push(HalfSpace::bisector(g));
        }

        let mut faces = candidate_faces;
        let mut vertices: Vec<Vec3> = Vec::new();
        for _pass in 0..2 {
            let nf = faces.len();
            vertices.clear();
            for i in 0..nf {
                for j in (i + 1)..nf {
                    for k in (j + 1)..nf {
                        if let Some(v) = intersect_three(&faces[i], &faces[j], &faces[k]) {
                            let inside = faces
                                .iter()
                                .all(|f| f.outside_by(v) < 1e-7 * (1.0 + f.offset.abs().sqrt()));
                            let inside_dbg = std::env::var("TPT_DBG").is_ok();
                            if inside_dbg && v.x.abs() > 0.7 && v.y.abs() > 0.7 && v.z.abs() > 0.7 {
                                eprintln!(
                                    "GEN ({:.6},{:.6},{:.6}) inside={}",
                                    v.x, v.y, v.z, inside
                                );
                            }
                            if inside && !vertices.iter().any(|u| (*u - v).norm() < 1e-7) {
                                vertices.push(v);
                            }
                        }
                    }
                }
            }
            // Keep faces with ≥3 vertices on their plane.
            faces.retain(|f| {
                vertices
                    .iter()
                    .filter(|v| (f.normal.dot(**v) - f.offset).abs() < 1e-6)
                    .count()
                    >= 3
            });
        }
        // A true vertex has ≥3 ON-plane normals spanning 3D (rank 3).
        // Points where a tangent (redundant) plane pierces a face have
        // rank-≤2 normal sets and are dropped.
        let on_normals = |v: &Vec3| -> Vec<Vec3> {
            faces
                .iter()
                .filter(|f| (f.normal.dot(*v) - f.offset).abs() < 1e-6)
                .map(|f| f.normal)
                .collect()
        };
        let rank3 = |normals: &[Vec3]| -> bool {
            normals.iter().any(|a| {
                normals.iter().any(|b| {
                    normals.iter().any(|c| {
                        a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                            + a.z * (b.x * c.y - b.y * c.x)
                            != 0.0
                    })
                })
            })
        };
        vertices.retain(|v| {
            let n = on_normals(v);
            n.len() >= 3 && rank3(&n)
        });
        // Drop faces whose plane no longer carries ≥3 of the true
        // vertices (tangent planes pierce faces but don't bound them).
        faces.retain(|f| {
            vertices
                .iter()
                .filter(|v| (f.normal.dot(**v) - f.offset).abs() < 1e-6)
                .count()
                >= 3
        });
        WignerSeitzCell { faces, vertices }
    }

    /// Volume via per-face triangle fans around the face centroid.
    pub fn volume(&self) -> f64 {
        let mut vol = 0.0;
        for f in &self.faces {
            let on_face: Vec<Vec3> = self
                .vertices
                .iter()
                .copied()
                .filter(|v| (f.normal.dot(*v) - f.offset).abs() < 1e-6)
                .collect();
            if on_face.len() < 3 {
                continue;
            }
            let normal = f.normal.normalize();
            let center = on_face.iter().fold(Vec3::ZERO, |a, v| a + *v) / on_face.len() as f64;
            let u0 = (on_face[0] - center).normalize();
            let u1 = normal.cross(u0);
            let mut sorted = on_face.clone();
            sorted.sort_by(|a, b| {
                let aa = (*a - center).normalize();
                let bb = (*b - center).normalize();
                let ang_a = num::atan2(aa.dot(u1), aa.dot(u0));
                let ang_b = num::atan2(bb.dot(u1), bb.dot(u0));
                ang_a.total_cmp(&ang_b)
            });
            for i in 0..sorted.len() {
                let a = sorted[i];
                let b = sorted[(i + 1) % sorted.len()];
                vol += a.dot(b.cross(center)) / 6.0;
            }
        }
        vol.abs()
    }
}

/// Intersection point of three non-parallel planes (Cramer over the
/// matrix whose ROWS are the normals).
fn intersect_three(a: &HalfSpace, b: &HalfSpace, c: &HalfSpace) -> Option<Vec3> {
    let det = a.normal.x * (b.normal.y * c.normal.z - b.normal.z * c.normal.y)
        - a.normal.y * (b.normal.x * c.normal.z - b.normal.z * c.normal.x)
        + a.normal.z * (b.normal.x * c.normal.y - b.normal.y * c.normal.x);
    if det.abs() < 1e-12 {
        return None;
    }
    let (dx, dy, dz) = (a.offset, b.offset, c.offset);
    // Cofactor expansion of Cramer's rule along the d-column.
    let x = dx * (b.normal.y * c.normal.z - b.normal.z * c.normal.y)
        - dy * (b.normal.x * c.normal.z - b.normal.z * c.normal.x)
        + dz * (b.normal.x * c.normal.y - b.normal.y * c.normal.x);
    let y = -dx * (b.normal.x * c.normal.z - b.normal.z * c.normal.x)
        + dy * (a.normal.x * c.normal.z - a.normal.z * c.normal.x)
        - dz * (a.normal.x * c.normal.z - a.normal.z * c.normal.x);
    let z = dx * (b.normal.x * c.normal.y - b.normal.y * c.normal.x)
        - dy * (a.normal.x * c.normal.y - a.normal.y * c.normal.x)
        + dz * (a.normal.x * b.normal.y - a.normal.y * b.normal.x);
    Some(Vec3::new(x / det, y / det, z / det))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cubic_bzone_is_a_cube() {
        // Simple cubic reciprocal lattice: BZ is a cube of side 2π/a.
        let lat = Lattice::cubic(4.0);
        let cell = WignerSeitzCell::brillouin_zone(&lat);
        // 6 faces, 8 vertices.
        assert_eq!(cell.faces.len(), 6, "faces: {:?}", cell.faces.len());
        assert_eq!(cell.vertices.len(), 8);
        // Vertices at (±π/a, ±π/a, ±π/a).
        let want = core::f64::consts::PI / 4.0;
        for v in &cell.vertices {
            for c in [v.x, v.y, v.z] {
                assert!((c.abs() - want).abs() < 1e-9, "vertex {v:?}");
            }
        }
        // Volume = (2π/a)³.
        assert!((cell.volume() - (2.0 * core::f64::consts::PI / 4.0f64).powi(3)).abs() < 1e-9);
    }

    #[test]
    fn contains_origin() {
        let lat = Lattice::cubic(4.0);
        let cell = WignerSeitzCell::brillouin_zone(&lat);
        for f in &cell.faces {
            assert!(f.outside_by(Vec3::ZERO) <= 1e-12);
        }
    }

    #[test]
    fn volume_matches_lattice_cell_for_cubic() {
        // The first BZ volume equals the reciprocal-cell volume 8π³/V.
        let lat = Lattice::cubic(3.0);
        let cell = WignerSeitzCell::brillouin_zone(&lat);
        let recip_vol = (2.0 * core::f64::consts::PI).powi(3) / lat.volume();
        assert!((cell.volume() - recip_vol).abs() < 1e-9);
    }
}
