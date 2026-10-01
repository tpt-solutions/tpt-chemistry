//! Linked-cell spatial partitioning.
//!
//! Atoms are bucketed into a uniform grid of `nx × ny × nz` cells; pair
//! iteration touches each cell and its forward half-neighborhood
//! (13 neighbors + the cell itself), giving each unordered pair exactly
//! once — O(N) work for a fixed cutoff density.

use tpt_chem_core::vec3::Vec3;

use crate::box3::Box3;

/// Linked-cell grid over a periodic orthorhombic box.
#[derive(Clone, Debug)]
pub struct CellList {
    /// Cells along each axis (all ≥ 1).
    pub n_cells: [usize; 3],
    /// Cell edge lengths (Å).
    cell_size: Vec3,
    box_: Box3,
    /// Head atom index per cell, `-1` for empty. Flattened index
    /// `k = (i·ny + j)·nz + k` — see [`CellList::cell_index`].
    head: Vec<i32>,
    /// Next atom in the same cell, `-1` terminated.
    next: Vec<i32>,
}

impl CellList {
    /// Grid dimensions covering a `cutoff` in a box of the given size.
    ///
    /// The grid is chosen as large as possible while keeping every cell at
    /// least `cutoff` wide, capped at 1 cell per axis below that.
    pub fn dimensions(box_: Box3, cutoff: f64) -> [usize; 3] {
        let n = |l: f64| {
            if cutoff > 0.0 && l > 0.0 {
                ((l / cutoff).floor() as usize).max(1)
            } else {
                1
            }
        };
        [n(box_.length.x), n(box_.length.y), n(box_.length.z)]
    }

    /// Build the list with an explicit grid over the given (wrapped)
    /// positions. Used by the neighbor list, which sizes cells to the
    /// interaction cutoff.
    pub fn build(box_: Box3, n_cells: [usize; 3], positions: &[Vec3]) -> Self {
        let mut list = CellList {
            n_cells: [n_cells[0].max(1), n_cells[1].max(1), n_cells[2].max(1)],
            cell_size: Vec3::new(
                box_.length.x / n_cells[0].max(1) as f64,
                box_.length.y / n_cells[1].max(1) as f64,
                box_.length.z / n_cells[2].max(1) as f64,
            ),
            box_,
            head: Vec::new(),
            next: Vec::new(),
        };
        let total = list.total_cells();
        list.head = vec![-1; total];
        list.next = vec![-1; positions.len()];
        for (idx, p) in positions.iter().enumerate() {
            let (i, j, k) = list.cell_of(p);
            let c = list.cell_index(i, j, k);
            list.next[idx] = list.head[c];
            list.head[c] = idx as i32;
        }
        list
    }

    /// Build the list with a heuristic grid (~one atom per cell) for the
    /// given wrapped positions.
    pub fn new(box_: Box3, positions: &[Vec3]) -> Self {
        let target = (positions.len() as f64).cbrt().max(1.0);
        let want = |l: f64| {
            if l > 0.0 {
                ((l / target).floor() as usize).clamp(1, 256)
            } else {
                1
            }
        };
        let n = [
            want(box_.length.x),
            want(box_.length.y),
            want(box_.length.z),
        ];
        Self::build(box_, n, positions)
    }

    /// The box this list partitions.
    pub fn box_(&self) -> Box3 {
        self.box_
    }

    /// Flat cell index for grid coordinates. The caller must keep
    /// `i < nx`, `j < ny`, `k < nz`; the flattening itself is
    /// `k_total = (i·ny + j)·nz + k < nx·ny·nz` — this invariant is
    /// formally proven in `tpt-chem-verify` (Kani harness
    /// `cell_index_in_bounds`).
    #[inline]
    pub fn cell_index(&self, i: usize, j: usize, k: usize) -> usize {
        (i * self.n_cells[1] + j) * self.n_cells[2] + k
    }

    /// Total number of cells.
    pub fn total_cells(&self) -> usize {
        self.n_cells[0] * self.n_cells[1] * self.n_cells[2]
    }

    /// Grid coordinates of a position (already wrapped).
    #[inline]
    pub fn cell_of(&self, p: &Vec3) -> (usize, usize, usize) {
        let clamp = |x: f64, n_cells: usize, size: f64| -> usize {
            let c = (x / size).floor();
            // Saturate into [0, n_cells-1] for any finite input.
            if c < 0.0 {
                0
            } else if c >= n_cells as f64 {
                n_cells - 1
            } else {
                c as usize
            }
        };
        (
            clamp(p.x, self.n_cells[0], self.cell_size.x),
            clamp(p.y, self.n_cells[1], self.cell_size.y),
            clamp(p.z, self.n_cells[2], self.cell_size.z),
        )
    }

    /// Invoke `f(a, b)` for every unordered atom pair `(a, b)` whose cells
    /// are within one cell of each other (the caller applies the actual
    /// cutoff). Neighbor cells wrap periodically.
    ///
    /// Each unordered cell pair is visited exactly once for ANY grid
    /// dimension (including degenerate 1- and 2-cell axes): the 26
    /// neighbor offsets are enumerated and a cross-cell visit happens only
    /// when the wrapped target cell has a strictly larger flat index.
    pub fn for_each_pair<F: FnMut(usize, usize)>(&self, mut f: F) {
        let [nx, ny, nz] = self.n_cells;
        for ci in 0..nx {
            for cj in 0..ny {
                for ck in 0..nz {
                    let c0 = self.cell_index(ci, cj, ck);
                    // Same cell: all internal pairs.
                    let mut a = self.head[c0];
                    while a >= 0 {
                        let mut b = self.next[a as usize];
                        while b >= 0 {
                            f(a as usize, b as usize);
                            b = self.next[b as usize];
                        }
                        a = self.next[a as usize];
                    }
                    // 26 neighbors. Distinct offsets can wrap onto the
                    // same cell on small axes (e.g. +z and −z with nz = 2),
                    // so deduplicate the wrapped targets first; the flat
                    // index ordering then visits each unordered cell pair
                    // exactly once.
                    let mut targets = [usize::MAX; 26];
                    let mut n_targets = 0;
                    for di in [-1i32, 0, 1] {
                        for dj in [-1i32, 0, 1] {
                            for dk in [-1i32, 0, 1] {
                                if di == 0 && dj == 0 && dk == 0 {
                                    continue;
                                }
                                let cn = self.cell_index(
                                    wrap_step(ci, di, nx),
                                    wrap_step(cj, dj, ny),
                                    wrap_step(ck, dk, nz),
                                );
                                if cn == c0 || targets[..n_targets].contains(&cn) {
                                    continue;
                                }
                                targets[n_targets] = cn;
                                n_targets += 1;
                            }
                        }
                    }
                    for &cn in &targets[..n_targets] {
                        if cn <= c0 {
                            continue; // the partner cell visits this pair
                        }
                        let mut a = self.head[c0];
                        while a >= 0 {
                            let mut b = self.head[cn];
                            while b >= 0 {
                                f(a as usize, b as usize);
                                b = self.next[b as usize];
                            }
                            a = self.next[a as usize];
                        }
                    }
                }
            }
        }
    }
}

/// `(c + d) mod n` with wrapping (signed arithmetic, always in range).
#[inline]
pub(crate) fn wrap_step(c: usize, d: i32, n: usize) -> usize {
    (c as i64 + d as i64).rem_euclid(n as i64) as usize
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Number of unique pairs visited by the list must match the exact
    /// count C(N,2).
    #[test]
    fn visits_every_pair_exactly_once() {
        let box_ = Box3::cubic(10.0);
        let mut positions = Vec::new();
        // 3×3×3 = 27 atoms on a regular lattice — enough cells to exercise
        // wrapping in every direction.
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    positions.push(Vec3::new(
                        (i as f64 + 0.5) * 10.0 / 3.0,
                        (j as f64 + 0.5) * 10.0 / 3.0,
                        (k as f64 + 0.5) * 10.0 / 3.0,
                    ));
                }
            }
        }
        let list = CellList::new(box_, &positions);
        let n = positions.len();
        let mut count = 0usize;
        let mut seen = vec![false; n * n];
        list.for_each_pair(|a, b| {
            count += 1;
            seen[a * n + b] = true;
            assert!(!seen[b * n + a], "pair ({a},{b}) visited twice");
        });
        assert_eq!(count, n * (n - 1) / 2);
    }

    #[test]
    fn dimensions_respect_cutoff() {
        let d = CellList::dimensions(Box3::cubic(10.0), 2.5);
        assert_eq!(d, [4, 4, 4]); // floor(10/2.5) = 4 cells of 2.5 Å
        let d = CellList::dimensions(Box3::cubic(10.0), 20.0);
        assert_eq!(d, [1, 1, 1]); // cutoff larger than box
    }

    #[test]
    fn cell_of_matches_bucket() {
        let box_ = Box3::cubic(9.0);
        let mut list = CellList {
            n_cells: [3, 3, 3],
            cell_size: Vec3::new(3.0, 3.0, 3.0),
            box_,
            head: vec![-1; 27],
            next: Vec::new(),
        };
        list.head = vec![-1; list.total_cells()];
        let p = Vec3::new(8.999, 0.0, 2.999);
        let (i, j, k) = list.cell_of(&p);
        assert_eq!((i, j, k), (2, 0, 0));
        assert!(list.cell_index(i, j, k) < list.total_cells());
    }
}
