//! Nonbonded exclusion lists with 1-4 scaling.
//!
//! Atoms joined by one or two bonds (1-2, 1-3) do not interact through
//! the pair potential; atoms three bonds apart (1-4) interact with
//! scaled Lennard-Jones and Coulomb terms (the AMBER/OPLS fudge
//! factors). The pair kernels in [`crate::forces`] consult
//! [`crate::system::System::exclusions`].

use std::collections::{HashMap, VecDeque};

use crate::bonded::Bonded;

/// Scale factors applied to one pair's LJ and Coulomb interaction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PairScale {
    /// Lennard-Jones scale (0 = fully excluded).
    pub lj: f64,
    /// Coulomb scale (0 = fully excluded).
    pub coulomb: f64,
}

/// Per-pair scale factors; pairs not listed interact fully.
#[derive(Clone, Debug, Default)]
pub struct Exclusions {
    map: HashMap<(usize, usize), PairScale>,
}

fn key(a: usize, b: usize) -> (usize, usize) {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}

impl Exclusions {
    /// No exclusions.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the scale of one pair (replacing any previous entry).
    pub fn set(&mut self, a: usize, b: usize, scale: PairScale) {
        self.map.insert(key(a, b), scale);
    }

    /// Scale for a pair, if it is excluded or scaled.
    pub fn get(&self, a: usize, b: usize) -> Option<PairScale> {
        self.map.get(&key(a, b)).copied()
    }

    /// All listed pairs `(a, b, scale)` with `a < b`.
    pub fn iter(&self) -> impl Iterator<Item = (usize, usize, PairScale)> + '_ {
        self.map.iter().map(|(&(a, b), &s)| (a, b, s))
    }

    /// Number of listed pairs.
    pub fn len(&self) -> usize {
        self.map.len()
    }

    /// True when no pair is listed.
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }

    /// Derive exclusions from the bond graph of `bonded` over `n_atoms`
    /// atoms: 1-2 and 1-3 pairs are fully excluded, 1-4 pairs (shortest
    /// bond path of exactly three) get `fudge_lj` / `fudge_qq`. Pairs that
    /// are closer through a ring keep the shorter-path treatment.
    pub fn from_bonded(bonded: &Bonded, n_atoms: usize, fudge_lj: f64, fudge_qq: f64) -> Self {
        let mut adj = vec![Vec::new(); n_atoms];
        for ([a, b], _) in &bonded.bonds {
            adj[*a].push(*b);
            adj[*b].push(*a);
        }
        let mut ex = Exclusions::new();
        for start in 0..n_atoms {
            let mut dist = vec![usize::MAX; n_atoms];
            dist[start] = 0;
            let mut queue = VecDeque::from([start]);
            while let Some(u) = queue.pop_front() {
                if dist[u] == 3 {
                    continue;
                }
                for &v in &adj[u] {
                    if dist[v] == usize::MAX {
                        dist[v] = dist[u] + 1;
                        queue.push_back(v);
                    }
                }
            }
            for (other, &d) in dist.iter().enumerate() {
                if other > start && (1..=3).contains(&d) {
                    let scale = if d == 3 {
                        PairScale {
                            lj: fudge_lj,
                            coulomb: fudge_qq,
                        }
                    } else {
                        PairScale {
                            lj: 0.0,
                            coulomb: 0.0,
                        }
                    };
                    ex.set(start, other, scale);
                }
            }
        }
        ex
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_chem_core::forcefield::HarmonicBond;

    fn chain(n: usize) -> Bonded {
        let mut b = Bonded::new();
        for i in 0..n - 1 {
            b.bond(i, i + 1, HarmonicBond { k: 1.0, r0: 1.0 });
        }
        b
    }

    #[test]
    fn chain_classification() {
        let ex = Exclusions::from_bonded(&chain(6), 6, 0.5, 0.8333);
        assert_eq!(ex.get(0, 1).unwrap().lj, 0.0);
        assert_eq!(ex.get(2, 0).unwrap().coulomb, 0.0);
        assert_eq!(ex.get(0, 3).unwrap().lj, 0.5);
        assert_eq!(ex.get(3, 0).unwrap().coulomb, 0.8333);
        assert!(ex.get(0, 4).is_none());
        assert_eq!(ex.len(), 5 + 4 + 3);
    }

    #[test]
    fn four_ring_uses_shortest_path() {
        // In a 4-ring, atoms 0 and 2 are 1-3 (two bonds), never 1-4.
        let mut b = chain(4);
        b.bond(3, 0, HarmonicBond { k: 1.0, r0: 1.0 });
        let ex = Exclusions::from_bonded(&b, 4, 0.5, 0.5);
        assert_eq!(ex.get(0, 2).unwrap().lj, 0.0);
        assert_eq!(ex.get(1, 3).unwrap().lj, 0.0);
    }
}
