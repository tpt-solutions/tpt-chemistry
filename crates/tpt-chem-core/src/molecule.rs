//! Graph-based molecular representation.
//!
//! A [`Molecule`] is a vertex-labeled graph: atoms carry element/position/
//! partial-charge data, bonds are edges with an order. Atom insertion with
//! compile-time-known elements uses the const-generic API
//! `add_atom::<Z>()` (spec.txt §4); parsers that discover elements at
//! runtime use [`Molecule::add_atom_raw`].
//!
//! ```
//! use tpt_chem_core::molecule::{BondOrder, Molecule};
//!
//! let mut h2o = Molecule::new("water");
//! let o = h2o.add_atom::<8>([0.0, 0.0, 0.0].into());
//! let h1 = h2o.add_atom::<1>([0.0, 0.9238, -0.3402].into());
//! let h2 = h2o.add_atom::<1>([0.0, -0.9238, -0.3402].into());
//! h2o.add_bond(o, h1, BondOrder::Single).unwrap();
//! h2o.add_bond(o, h2, BondOrder::Single).unwrap();
//!
//! assert_eq!(h2o.degree(o), 2);
//! assert!((h2o.mass() - 18.015).abs() < 1e-3); // amu
//! ```

use alloc::string::{String, ToString};
use alloc::vec::Vec;
use core::fmt;

use crate::element::{self, MAX_Z};
use crate::vec3::Vec3;

/// Handle to an atom inside a [`Molecule`]. Indices are only handed out by
/// `add_atom*`, so a valid `AtomId` always addresses an existing atom.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AtomId(u16);

impl AtomId {
    /// Position of this atom in the molecule's atom list.
    pub fn index(self) -> usize {
        self.0 as usize
    }
}

impl From<u16> for AtomId {
    fn from(v: u16) -> Self {
        AtomId(v)
    }
}

impl From<usize> for AtomId {
    fn from(v: usize) -> Self {
        AtomId(v as u16)
    }
}

/// Chemical bond order.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BondOrder {
    /// Single bond (order 1).
    Single,
    /// Double bond (order 2).
    Double,
    /// Triple bond (order 3).
    Triple,
    /// Aromatic bond (delocalized, conventionally order 1.5).
    Aromatic,
}

impl BondOrder {
    /// Integral bond order; aromatic counts as 1.5.
    pub fn value(self) -> f64 {
        match self {
            BondOrder::Single => 1.0,
            BondOrder::Double => 2.0,
            BondOrder::Triple => 3.0,
            BondOrder::Aromatic => 1.5,
        }
    }
}

/// An atom record: element, Cartesian position (Å), partial charge (e), and
/// an optional label (e.g. PDB atom name).
#[derive(Clone, Debug, PartialEq)]
pub struct Atom {
    /// Atomic number.
    pub z: u8,
    /// Position in Å.
    pub pos: Vec3,
    /// Partial charge in units of the elementary charge.
    pub charge: f64,
    /// Optional atom label (e.g. `"CA"` in PDB files).
    pub label: Option<String>,
}

impl Atom {
    /// Standard atomic mass in amu.
    pub fn mass(&self) -> f64 {
        element::mass(self.z).unwrap_or(0.0)
    }
}

/// A bond between two atoms.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Bond {
    /// First endpoint (lower `AtomId`).
    pub a: AtomId,
    /// Second endpoint (higher `AtomId`).
    pub b: AtomId,
    /// Bond order.
    pub order: BondOrder,
}

/// Errors from molecule graph operations.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum MolError {
    /// Tried to bond an atom to itself.
    SelfBond,
    /// The two atoms are already bonded.
    DuplicateBond,
    /// An `AtomId` is out of range.
    BadAtomId,
}

impl fmt::Display for MolError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MolError::SelfBond => write!(f, "cannot bond an atom to itself"),
            MolError::DuplicateBond => write!(f, "bond already exists"),
            MolError::BadAtomId => write!(f, "atom id out of range"),
        }
    }
}

#[cfg(feature = "std")]
impl std::error::Error for MolError {}

/// A molecular graph: atoms + bonds, with adjacency lookups.
#[derive(Clone, Debug)]
pub struct Molecule {
    name: String,
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    adjacency: Vec<Vec<usize>>,
    counts: [u16; MAX_Z as usize],
}

impl Default for Molecule {
    fn default() -> Self {
        Molecule {
            name: String::new(),
            atoms: Vec::new(),
            bonds: Vec::new(),
            adjacency: Vec::new(),
            counts: [0; MAX_Z as usize],
        }
    }
}

impl Molecule {
    /// An empty molecule with the given name.
    pub fn new(name: &str) -> Self {
        Molecule {
            name: name.to_string(),
            ..Molecule::default()
        }
    }

    /// Add an atom whose element is known at compile time.
    ///
    /// ```
    /// use tpt_chem_core::element::Element;
    /// use tpt_chem_core::molecule::Molecule;
    ///
    /// let mut m = Molecule::new("methane");
    /// let c = m.add_atom::<6>([0.0; 3].into());
    /// assert_eq!(c.index(), 0);
    /// assert_eq!(Element::<6>::symbol(), "C");
    /// ```
    pub fn add_atom<const Z: u8>(&mut self, pos: Vec3) -> AtomId {
        self.add_atom_raw(Z, pos)
    }

    /// Add an atom with a runtime atomic number (parser path).
    pub fn add_atom_raw(&mut self, z: u8, pos: Vec3) -> AtomId {
        assert!(
            z >= 1 && (z as usize) <= MAX_Z as usize,
            "atomic number out of range"
        );
        let id = AtomId(self.atoms.len() as u16);
        self.atoms.push(Atom {
            z,
            pos,
            charge: 0.0,
            label: None,
        });
        self.adjacency.push(Vec::new());
        self.counts[(z - 1) as usize] += 1;
        id
    }

    /// Add a labeled atom (runtime element).
    pub fn add_atom_labeled(&mut self, z: u8, pos: Vec3, label: &str) -> AtomId {
        let id = self.add_atom_raw(z, pos);
        self.atoms[id.index()].label = Some(label.to_string());
        id
    }

    /// Add a bond between two atoms. The endpoints are stored in
    /// normalized order (lower id first).
    pub fn add_bond(&mut self, a: AtomId, b: AtomId, order: BondOrder) -> Result<(), MolError> {
        if a.index() >= self.atoms.len() || b.index() >= self.atoms.len() {
            return Err(MolError::BadAtomId);
        }
        if a == b {
            return Err(MolError::SelfBond);
        }
        let (lo, hi) = if a.index() < b.index() {
            (a, b)
        } else {
            (b, a)
        };
        if self.bond_between(lo, hi).is_some() {
            return Err(MolError::DuplicateBond);
        }
        let bi = self.bonds.len();
        self.bonds.push(Bond {
            a: lo,
            b: hi,
            order,
        });
        self.adjacency[lo.index()].push(bi);
        self.adjacency[hi.index()].push(bi);
        Ok(())
    }

    /// Number of atoms.
    pub fn len(&self) -> usize {
        self.atoms.len()
    }

    /// True if there are no atoms.
    pub fn is_empty(&self) -> bool {
        self.atoms.is_empty()
    }

    /// Molecule name.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Rename the molecule.
    pub fn set_name(&mut self, name: &str) {
        self.name = name.to_string();
    }

    /// Borrow atom `id`.
    pub fn atom(&self, id: AtomId) -> &Atom {
        &self.atoms[id.index()]
    }

    /// Borrow atom by list position, or `None` out of range.
    pub fn get_atom(&self, i: usize) -> Option<&Atom> {
        self.atoms.get(i)
    }

    /// Mutable borrow of atom `id`.
    pub fn atom_mut(&mut self, id: AtomId) -> &mut Atom {
        &mut self.atoms[id.index()]
    }

    /// All atoms, in insertion order.
    pub fn atoms(&self) -> impl Iterator<Item = &Atom> {
        self.atoms.iter()
    }

    /// All bonds, in insertion order.
    pub fn bonds(&self) -> impl Iterator<Item = &Bond> {
        self.bonds.iter()
    }

    /// Number of bonds.
    pub fn n_bonds(&self) -> usize {
        self.bonds.len()
    }

    /// The bond between two atoms, if any.
    pub fn bond_between(&self, a: AtomId, b: AtomId) -> Option<&Bond> {
        if a.index() >= self.atoms.len() || b.index() >= self.atoms.len() {
            return None;
        }
        let (lo, hi) = if a.index() < b.index() {
            (a, b)
        } else {
            (b, a)
        };
        self.adjacency[a.index()]
            .iter()
            .map(|&bi| &self.bonds[bi])
            .find(|bond| bond.a == lo && bond.b == hi)
    }

    /// Atoms directly bonded to `a`.
    pub fn neighbors(&self, a: AtomId) -> impl Iterator<Item = AtomId> + '_ {
        let a_index = a.index();
        self.adjacency[a_index].iter().map(move |&bi| {
            let bond = &self.bonds[bi];
            if bond.a == a {
                bond.b
            } else {
                bond.a
            }
        })
    }

    /// Number of bonds at atom `a`.
    pub fn degree(&self, a: AtomId) -> usize {
        self.adjacency[a.index()].len()
    }

    /// Distance between two atoms in Å.
    ///
    /// # Panics
    /// If either id is out of range.
    pub fn distance(&self, a: AtomId, b: AtomId) -> f64 {
        self.atoms[a.index()].pos.dist(self.atoms[b.index()].pos)
    }

    /// Bond angle at the middle atom `b` (a–b–c), in radians.
    ///
    /// # Panics
    /// If any id is out of range or the atoms are collinear with `a == c`.
    pub fn angle(&self, a: AtomId, b: AtomId, c: AtomId) -> f64 {
        let ab = self.atom(a).pos - self.atom(b).pos;
        let cb = self.atom(c).pos - self.atom(b).pos;
        let cos_theta = ab.dot(cb) / (ab.norm() * cb.norm());
        crate::num::acos(cos_theta)
    }

    /// Total standard mass in amu.
    pub fn mass(&self) -> f64 {
        self.atoms.iter().map(Atom::mass).sum()
    }

    /// Center of mass in Å.
    ///
    /// # Panics
    /// If the molecule is empty (total mass zero).
    pub fn center_of_mass(&self) -> Vec3 {
        let mut c = Vec3::ZERO;
        for a in &self.atoms {
            c += a.pos * a.mass();
        }
        c / self.mass()
    }

    /// Unweighted mean of the atomic positions in Å.
    ///
    /// # Panics
    /// If the molecule is empty.
    pub fn center(&self) -> Vec3 {
        let mut c = Vec3::ZERO;
        for a in &self.atoms {
            c += a.pos;
        }
        c / self.len() as f64
    }

    /// Sum of partial charges in units of e.
    pub fn total_charge(&self) -> f64 {
        self.atoms.iter().map(|a| a.charge).sum()
    }

    /// Hill-order molecular formula (C, H first, rest alphabetical).
    pub fn formula(&self) -> String {
        element::hill_formula(&self.counts)
    }

    /// Translate every atom by `d` (Å).
    pub fn translate(&mut self, d: Vec3) {
        for a in &mut self.atoms {
            a.pos += d;
        }
    }

    /// Move the center of mass to the origin.
    pub fn center_at_origin(&mut self) {
        let c = self.center_of_mass();
        self.translate(-c);
    }

    /// Infer bonds from geometry: two atoms are bonded if their distance is
    /// below `tolerance × (r_cov(i) + r_cov(j))` and above `0.4 Å`. Existing
    /// bonds are kept; returns how many bonds were added.
    pub fn infer_bonds(&mut self, tolerance: f64) -> usize {
        let n = self.atoms.len();
        let radii: Vec<f64> = (0..n)
            .map(|i| element::covalent_radius(self.atoms[i].z).unwrap_or(0.77))
            .collect();
        let mut added = 0;
        for i in 0..n {
            for j in (i + 1)..n {
                let d = self.atoms[i].pos.dist(self.atoms[j].pos);
                if d > 0.4
                    && d < tolerance * (radii[i] + radii[j])
                    && self
                        .add_bond(AtomId(i as u16), AtomId(j as u16), BondOrder::Single)
                        .is_ok()
                {
                    added += 1;
                }
            }
        }
        added
    }

    /// Remove all bonds.
    pub fn clear_bonds(&mut self) {
        self.bonds.clear();
        for adj in &mut self.adjacency {
            adj.clear();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::element::Element;

    fn c2h4() -> Molecule {
        let mut m = Molecule::new("C2H4 skeleton");
        let c1 = m.add_atom::<6>([0.0, 0.0, 0.77].into());
        let c2 = m.add_atom::<6>([0.0, 0.0, -0.77].into());
        m.add_bond(c1, c2, BondOrder::Single).unwrap();
        // Two H on each carbon, mirrored about z = 0.
        for (hx, hy, z, target) in [
            (-0.51, 0.89, 1.35, c1),
            (0.51, -0.89, 1.35, c1),
            (-0.51, 0.89, -1.35, c2),
            (0.51, -0.89, -1.35, c2),
        ] {
            let h = m.add_atom::<1>(Vec3::new(hx, hy, z));
            m.add_bond(h, target, BondOrder::Single).unwrap();
        }
        m
    }

    #[test]
    fn add_and_query() {
        let m = c2h4();
        assert_eq!(m.len(), 6);
        assert_eq!(m.n_bonds(), 5);
        assert_eq!(m.degree(AtomId(0)), 3);
        assert_eq!(m.neighbors(AtomId(0)).count(), 3);
        assert!(m.bond_between(AtomId(0), AtomId(1)).is_some());
        assert!(m.bond_between(AtomId(1), AtomId(0)).is_some());
        assert!(m.bond_between(AtomId(0), AtomId(2)).is_some());
        assert!(m.bond_between(AtomId(2), AtomId(4)).is_none());
    }

    #[test]
    fn bond_errors() {
        let mut m = Molecule::new("x");
        let a = m.add_atom::<1>(Vec3::ZERO);
        let b = m.add_atom::<1>(Vec3::new(1.0, 0.0, 0.0));
        assert_eq!(m.add_bond(a, a, BondOrder::Single), Err(MolError::SelfBond));
        m.add_bond(a, b, BondOrder::Single).unwrap();
        assert_eq!(
            m.add_bond(a, b, BondOrder::Single),
            Err(MolError::DuplicateBond)
        );
        assert_eq!(
            m.add_bond(b, a, BondOrder::Single),
            Err(MolError::DuplicateBond)
        );
        assert_eq!(
            m.add_bond(AtomId(99), b, BondOrder::Single),
            Err(MolError::BadAtomId)
        );
    }

    #[test]
    fn mass_formula_charge() {
        let mut m = c2h4();
        let expect = 2.0 * Element::<6>::mass() + 4.0 * Element::<1>::mass();
        assert!((m.mass() - expect).abs() < 1e-9);
        assert_eq!(m.formula(), "C2H4");
        assert!(m.total_charge().abs() < 1e-12);
        m.atom_mut(AtomId(0)).charge = 0.5;
        assert!((m.total_charge() - 0.5).abs() < 1e-12);
    }

    #[test]
    fn center_of_mass_midpoint() {
        let m = c2h4();
        let com = m.center_of_mass();
        assert!(com.x.abs() < 1e-12 && com.y.abs() < 1e-12 && com.z.abs() < 1e-12);
    }

    #[test]
    fn infer_bonds_water() {
        let mut m = Molecule::new("water");
        m.add_atom::<8>([0.0, 0.0, 0.117].into());
        m.add_atom::<1>([0.0, 0.757, -0.469].into());
        m.add_atom::<1>([0.0, -0.757, -0.469].into());
        assert_eq!(m.infer_bonds(1.3), 2);
        assert_eq!(m.degree(AtomId(0)), 2);
        // Re-running adds nothing (duplicates rejected).
        assert_eq!(m.infer_bonds(1.3), 0);
    }

    #[test]
    fn angle() {
        let mut m = Molecule::new("linear");
        m.add_atom::<1>([0.0, 0.0, -1.0].into());
        m.add_atom::<1>([0.0, 0.0, 0.0].into());
        m.add_atom::<1>([0.0, 0.0, 1.0].into());
        let th = m.angle(AtomId(0), AtomId(1), AtomId(2));
        assert!((th - core::f64::consts::PI).abs() < 1e-12);
    }
}
