//! Verlet neighbor lists with a skin region.
//!
//! The list stores all pairs within `cutoff + skin`; forces only consider
//! the true cutoff. Between rebuilds, atom displacement is tracked and the
//! list is refreshed once any atom has moved more than `skin/2` since the
//! last build — the standard half-skin criterion guaranteeing that no pair
//! inside the cutoff can be missed.

use tpt_chem_core::vec3::Vec3;

use crate::box3::Box3;
use crate::cell::CellList;

/// A neighbor pair list with displacement-triggered rebuilds.
#[derive(Clone, Debug)]
pub struct VerletList {
    /// All pairs `(a, b)`, `a < b`, within `cutoff + skin`.
    pub pairs: Vec<(usize, usize)>,
    /// Interaction cutoff (Å).
    pub cutoff: f64,
    /// Skin width (Å).
    pub skin: f64,
    /// Positions at the time of the last build.
    reference: Vec<Vec3>,
    /// Whether the last build used periodic wrapping.
    periodic: bool,
}

impl VerletList {
    /// Build a neighbor list for `positions`.
    ///
    /// * `cutoff` — interaction distance (Å).
    /// * `skin` — safety margin; `0.0` means rebuild every step.
    /// * `box_` — `Some` for periodic minimum-image distances, `None` for
    ///   open boundaries.
    pub fn build(positions: &[Vec3], cutoff: f64, skin: f64, box_: Option<Box3>) -> Self {
        let reach = cutoff + skin;
        let pairs = if let Some(b) = box_ {
            let dims = CellList::dimensions(b, reach.max(1e-6));
            let cells = CellList::build(b, dims, positions);
            collect_pairs(positions, &cells, reach, true, b)
        } else {
            // Open boundaries: an O(N²) scan is exact and simple; the
            // cell-list path is the hot one for periodic systems.
            let mut pairs = Vec::new();
            let reach2 = reach * reach;
            for a in 0..positions.len() {
                for b in (a + 1)..positions.len() {
                    if positions[a].dist_sq(positions[b]) < reach2 {
                        pairs.push((a, b));
                    }
                }
            }
            pairs
        };
        VerletList {
            pairs,
            cutoff,
            skin,
            reference: positions.to_vec(),
            periodic: box_.is_some(),
        }
    }

    /// Whether any atom has moved more than `skin/2` since the last build.
    ///
    /// With periodic boundaries, displacement is measured with the
    /// minimum-image convention so slow drift across a boundary does not
    /// trigger spurious rebuilds.
    pub fn needs_update(&self, positions: &[Vec3], box_: Option<Box3>) -> bool {
        let threshold = 0.25 * self.skin * self.skin; // (skin/2)²
        if self.reference.len() != positions.len() {
            return true;
        }
        for (r, p) in self.reference.iter().zip(positions.iter()) {
            let d = match (self.periodic, box_) {
                (true, Some(b)) => b.min_image(*p - *r),
                _ => *p - *r,
            };
            if d.norm_sq() > threshold {
                return true;
            }
        }
        false
    }

    /// Rebuild in place (same parameters as the constructor).
    pub fn update(&mut self, positions: &[Vec3], box_: Option<Box3>) {
        *self = VerletList::build(positions, self.cutoff, self.skin, box_);
    }
}

fn collect_pairs(
    positions: &[Vec3],
    cells: &CellList,
    reach: f64,
    periodic: bool,
    box_: Box3,
) -> Vec<(usize, usize)> {
    let _ = periodic;
    let reach2 = reach * reach;
    let mut pairs = Vec::new();
    cells.for_each_pair(|a, b| {
        let d = if positions[a].x.is_finite() && positions[b].x.is_finite() {
            box_.min_image(positions[b] - positions[a])
        } else {
            positions[b] - positions[a]
        };
        if d.norm_sq() < reach2 {
            pairs.push((a.min(b), a.max(b)));
        }
    });
    pairs
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn open_system_pair_count() {
        // Two atoms 3 Å apart, cutoff 4, skin 1 → one pair.
        let positions = [Vec3::ZERO, Vec3::new(3.0, 0.0, 0.0)];
        let vl = VerletList::build(&positions, 4.0, 1.0, None);
        assert_eq!(vl.pairs, vec![(0, 1)]);

        // Move the second atom to 6 Å — inside reach (5), list stale but
        // `needs_update` must flag it.
        let moved = [Vec3::ZERO, Vec3::new(6.0, 0.0, 0.0)];
        assert!(vl.needs_update(&moved, None));
        let updated = VerletList::build(&moved, 4.0, 1.0, None);
        assert!(updated.pairs.is_empty());
    }

    #[test]
    fn periodic_pair_count_matches_bruteforce() {
        let b = Box3::cubic(12.0);
        let mut positions = Vec::new();
        // 4×4×4 lattice.
        for i in 0..4 {
            for j in 0..4 {
                for k in 0..4 {
                    positions.push(Vec3::new(
                        i as f64 * 3.0 + 1.0,
                        j as f64 * 3.0 + 1.0,
                        k as f64 * 3.0 + 1.0,
                    ));
                }
            }
        }
        let vl = VerletList::build(&positions, 3.3, 0.5, Some(b));
        // Brute-force reference.
        let mut expected = 0;
        let reach2 = 3.8 * 3.8;
        for a in 0..positions.len() {
            for c in (a + 1)..positions.len() {
                if b.min_image(positions[c] - positions[a]).norm_sq() < reach2 {
                    expected += 1;
                }
            }
        }
        assert_eq!(vl.pairs.len(), expected);
    }

    #[test]
    fn no_rebuild_for_small_motion() {
        let positions = [Vec3::ZERO, Vec3::new(3.0, 0.0, 0.0)];
        let vl = VerletList::build(&positions, 4.0, 1.0, None);
        let bumped = [Vec3::new(0.1, 0.0, 0.0), Vec3::new(3.05, 0.0, 0.0)];
        assert!(!vl.needs_update(&bumped, None));
    }
}
