//! SMILES reader: graph and composition from a line-notation string.
//!
//! Supports the organic subset (`B C N O P S F Cl Br I`, with implicit
//! hydrogens from default valences), bracket atoms (`[NH4+]`, `[Na+]`,
//! isotopes and charges), aromatic lowercase atoms, ring closures (`1`–`9`
//! and `%nn`), branches, bond symbols `- = # :`, and `.` for disconnected
//! components.
//!
//! The result has no 3D coordinates: every atom is placed at the origin.
//! Use it for formulas, connectivity and graph analysis, or embed the
//! geometry in a later step.
//!
//! Not handled: stereochemistry (`@`, `/`, `\` are accepted and ignored),
//! kekulization (aromatic bonds stay [`BondOrder::Aromatic`]), and atoms
//! outside the organic subset unless written in brackets.
//!
//! # Examples
//!
//! ```
//! use tpt_chem_io::parse_smiles;
//!
//! let benzene = parse_smiles("c1ccccc1").unwrap();
//! assert_eq!(benzene.formula(), "C6H6");
//! ```

use std::collections::BTreeMap;

use tpt_chem_core::element;
use tpt_chem_core::molecule::{AtomId, BondOrder, Molecule};
use tpt_chem_core::vec3::Vec3;

use crate::error::{IoError, Result};

const FORMAT: &str = "SMILES";

/// Parse a SMILES string into a molecule graph (no coordinates).
///
/// Only the first whitespace-separated token is read, so a trailing name or
/// comment after the SMILES is ignored, as is conventional.
///
/// # Errors
/// [`IoError`] on empty input, unbalanced branches or ring closures,
/// unknown elements, or malformed bracket atoms.
pub fn parse_smiles(input: &str) -> Result<Molecule> {
    let text = input
        .split_whitespace()
        .next()
        .ok_or(IoError::UnexpectedEof { format: FORMAT })?;
    let mut parser = Parser {
        s: text.as_bytes(),
        pos: 0,
        mol: Molecule::new(text),
        info: Vec::new(),
        current: None,
        branches: Vec::new(),
        pending_bond: None,
        rings: BTreeMap::new(),
    };
    parser.run()?;
    parser.finish()
}

/// Per-atom bookkeeping needed after the heavy-atom graph is built.
struct AtomInfo {
    id: AtomId,
    z: u8,
    aromatic: bool,
    /// Hydrogens written in brackets (`[NH4+]`); `None` for organic-subset atoms.
    bracket_h: Option<u32>,
    charge: i32,
}

struct Parser<'a> {
    s: &'a [u8],
    pos: usize,
    mol: Molecule,
    info: Vec<AtomInfo>,
    /// Atom that the next atom will bond to.
    current: Option<AtomId>,
    /// Branch points saved at each `(`.
    branches: Vec<AtomId>,
    /// Bond symbol written immediately before the next atom or ring digit.
    pending_bond: Option<BondOrder>,
    /// Open ring closures: digit → (opening atom, bond symbol at the opener).
    rings: BTreeMap<u32, (AtomId, Option<BondOrder>)>,
}

impl Parser<'_> {
    fn err(&self, message: &str) -> IoError {
        IoError::Parse {
            format: FORMAT,
            line: None,
            message: format!("{message} at column {}", self.pos + 1),
        }
    }

    fn peek(&self) -> Option<u8> {
        self.s.get(self.pos).copied()
    }

    fn run(&mut self) -> Result<()> {
        while let Some(c) = self.peek() {
            match c {
                b'[' => self.bracket_atom()?,
                b'B' | b'C' | b'N' | b'O' | b'P' | b'S' | b'F' | b'I' | b'b' | b'c' | b'n'
                | b'o' | b'p' | b's' => self.organic_atom()?,
                b'(' => {
                    let cur = self
                        .current
                        .ok_or_else(|| self.err("'(' with no atom before"))?;
                    self.branches.push(cur);
                    self.pos += 1;
                }
                b')' => {
                    let back = self
                        .branches
                        .pop()
                        .ok_or_else(|| self.err("')' without matching '('"))?;
                    self.current = Some(back);
                    self.pos += 1;
                }
                b'.' => {
                    self.current = None;
                    self.pending_bond = None;
                    self.pos += 1;
                }
                b'-' => self.bond_symbol(BondOrder::Single),
                b'=' => self.bond_symbol(BondOrder::Double),
                b'#' => self.bond_symbol(BondOrder::Triple),
                b':' => self.bond_symbol(BondOrder::Aromatic),
                b'/' | b'\\' => {
                    // Directional bonds only matter for stereochemistry, which is ignored.
                    self.pending_bond = Some(BondOrder::Single);
                    self.pos += 1;
                }
                b'0'..=b'9' | b'%' => self.ring_bond()?,
                _ => return Err(self.err(&format!("unexpected character '{}'", c as char))),
            }
        }
        Ok(())
    }

    fn bond_symbol(&mut self, order: BondOrder) {
        self.pending_bond = Some(order);
        self.pos += 1;
    }

    /// Default bond between two atoms when no symbol is written.
    fn default_order(&self, a: AtomId, b: AtomId) -> BondOrder {
        let aromatic = |id: AtomId| self.info.iter().any(|i| i.id == id && i.aromatic);
        if aromatic(a) && aromatic(b) {
            BondOrder::Aromatic
        } else {
            BondOrder::Single
        }
    }

    fn connect(&mut self, a: AtomId, b: AtomId, order: BondOrder) -> Result<()> {
        self.mol
            .add_bond(a, b, order)
            .map_err(|e| self.err(&format!("{e:?}")))
    }

    /// Add an atom bonded to the current one (if any) and make it current.
    fn push_atom(
        &mut self,
        z: u8,
        aromatic: bool,
        bracket_h: Option<u32>,
        charge: i32,
    ) -> Result<()> {
        let id = self.mol.add_atom_raw(z, Vec3::ZERO);
        if let Some(prev) = self.current {
            let order = match self.pending_bond.take() {
                Some(order) => order,
                None => {
                    let prev_aromatic = self.info.iter().any(|i| i.id == prev && i.aromatic);
                    if aromatic && prev_aromatic {
                        BondOrder::Aromatic
                    } else {
                        BondOrder::Single
                    }
                }
            };
            self.connect(prev, id, order)?;
        }
        self.pending_bond = None;
        self.info.push(AtomInfo {
            id,
            z,
            aromatic,
            bracket_h,
            charge,
        });
        self.current = Some(id);
        Ok(())
    }

    fn organic_atom(&mut self) -> Result<()> {
        let c = self.s[self.pos];
        let (z, len, aromatic) = if c.is_ascii_lowercase() {
            let z = match c {
                b'b' => 5,
                b'c' => 6,
                b'n' => 7,
                b'o' => 8,
                b'p' => 15,
                b's' => 16,
                _ => unreachable!("run() only dispatches aromatic organic letters"),
            };
            (z, 1, true)
        } else if self.s[self.pos..].starts_with(b"Cl") {
            (17, 2, false)
        } else if self.s[self.pos..].starts_with(b"Br") {
            (35, 2, false)
        } else {
            let z = match c {
                b'B' => 5,
                b'C' => 6,
                b'N' => 7,
                b'O' => 8,
                b'P' => 15,
                b'S' => 16,
                b'F' => 9,
                b'I' => 53,
                _ => unreachable!("run() only dispatches organic letters"),
            };
            (z, 1, false)
        };
        self.pos += len;
        self.push_atom(z, aromatic, None, 0)
    }

    fn bracket_atom(&mut self) -> Result<()> {
        self.pos += 1;
        // Isotope mass numbers are accepted and ignored.
        while self.peek().is_some_and(|c| c.is_ascii_digit()) {
            self.pos += 1;
        }

        let c = self
            .peek()
            .filter(u8::is_ascii_alphabetic)
            .ok_or_else(|| self.err("expected an element symbol in brackets"))?;
        let aromatic = c.is_ascii_lowercase();
        let start = self.pos;
        self.pos += 1;
        // Prefer a two-letter symbol when it names a real element (Cl, Se, Na).
        let two = match self.s.get(start..start + 2) {
            Some([first, second]) if second.is_ascii_lowercase() => {
                let mut symbol = [*first, *second];
                symbol[0] = symbol[0].to_ascii_uppercase();
                element_from_bytes(&symbol)
            }
            _ => None,
        };
        let z = match two {
            Some(z) => {
                self.pos += 1;
                z
            }
            None => {
                let upper = [c.to_ascii_uppercase()];
                if aromatic && !b"bcnops".contains(&c) {
                    return Err(self.err("aromatic symbol outside b c n o p s"));
                }
                element_from_bytes(&upper)
                    .ok_or_else(|| self.err("unknown element symbol in brackets"))?
            }
        };

        // Chirality marks are accepted and ignored.
        while self.peek() == Some(b'@') {
            self.pos += 1;
        }

        let mut bracket_h = 0u32;
        if self.peek() == Some(b'H') {
            self.pos += 1;
            bracket_h = 1;
            if let Some(n) = self.read_number() {
                bracket_h = n;
            }
        }

        let mut charge = 0i32;
        if let Some(sign) = self.peek().filter(|c| *c == b'+' || *c == b'-') {
            let unit = if sign == b'+' { 1 } else { -1 };
            self.pos += 1;
            charge = match self.read_number() {
                Some(n) => unit * n as i32,
                None => {
                    // Repeated signs ("++") each add one unit.
                    let mut total = unit;
                    while self.peek() == Some(sign) {
                        self.pos += 1;
                        total += unit;
                    }
                    total
                }
            };
        }

        // Atom-class number after ':' (used in reaction SMILES) is ignored.
        if self.peek() == Some(b':') {
            self.pos += 1;
            self.read_number();
        }

        if self.peek() != Some(b']') {
            return Err(self.err("expected ']' closing bracket atom"));
        }
        self.pos += 1;
        self.push_atom(z, aromatic, Some(bracket_h), charge)
    }

    fn read_number(&mut self) -> Option<u32> {
        let start = self.pos;
        let mut value = 0u32;
        while let Some(d) = self.peek().filter(u8::is_ascii_digit) {
            value = value * 10 + u32::from(d - b'0');
            self.pos += 1;
        }
        (self.pos > start).then_some(value)
    }

    fn ring_bond(&mut self) -> Result<()> {
        let num = if self.peek() == Some(b'%') {
            self.pos += 1;
            let two_digits = self
                .s
                .get(self.pos..self.pos + 2)
                .filter(|d| d.iter().all(u8::is_ascii_digit))
                .ok_or_else(|| self.err("'%' must be followed by two digits"))?;
            self.pos += 2;
            u32::from(two_digits[0] - b'0') * 10 + u32::from(two_digits[1] - b'0')
        } else {
            let d = self.s[self.pos];
            self.pos += 1;
            u32::from(d - b'0')
        };
        let cur = self
            .current
            .ok_or_else(|| self.err("ring bond with no atom before it"))?;
        let explicit = self.pending_bond.take();
        if let Some((open, open_bond)) = self.rings.remove(&num) {
            let order = explicit
                .or(open_bond)
                .unwrap_or_else(|| self.default_order(open, cur));
            self.connect(open, cur, order)
        } else {
            self.rings.insert(num, (cur, explicit));
            Ok(())
        }
    }

    fn finish(mut self) -> Result<Molecule> {
        if !self.branches.is_empty() {
            return Err(self.err("unclosed '('"));
        }
        if !self.rings.is_empty() {
            return Err(self.err("unclosed ring bond"));
        }

        // Sum of bond orders per heavy atom, for implicit-hydrogen counts.
        let mut sums = vec![0.0f64; self.info.len()];
        let index_of = |id: AtomId| id.index();
        for bond in self.mol.bonds() {
            sums[index_of(bond.a)] += bond.order.value();
            sums[index_of(bond.b)] += bond.order.value();
        }

        let mut total_charge = 0i32;
        let mut hydrogens = Vec::with_capacity(self.info.len());
        for (i, info) in self.info.iter().enumerate() {
            total_charge += info.charge;
            let h = match info.bracket_h {
                Some(h) => h,
                None => {
                    let valence = default_valence(info.z)
                        .ok_or_else(|| self.err("organic-subset atom without a valence"))?;
                    // Aromatic bonds count 1.5, so floor() gives the integer
                    // bond count (benzene carbon: 3.0 → one H).
                    (valence - sums[i].floor() as i32).max(0) as u32
                }
            };
            hydrogens.push((info.id, h));
        }

        for (heavy, h) in hydrogens {
            for _ in 0..h {
                let hydrogen = self.mol.add_atom_raw(1, Vec3::ZERO);
                self.connect(heavy, hydrogen, BondOrder::Single)?;
            }
        }
        self.mol.set_formal_charge(total_charge);
        Ok(self.mol)
    }
}

/// Default valence for implicit hydrogens on organic-subset atoms.
fn default_valence(z: u8) -> Option<i32> {
    match z {
        5 => Some(3),  // B
        6 => Some(4),  // C
        7 => Some(3),  // N
        8 => Some(2),  // O
        9 => Some(1),  // F
        15 => Some(3), // P
        16 => Some(2), // S
        17 => Some(1), // Cl
        35 => Some(1), // Br
        53 => Some(1), // I
        _ => None,
    }
}

fn element_from_bytes(symbol: &[u8]) -> Option<u8> {
    let text = core::str::from_utf8(symbol).ok()?;
    element::from_symbol(text)
}

/// Write a molecule as canonical SMILES.
///
/// Atoms are ranked from graph invariants (element, degree, hydrogen count,
/// aromaticity), refined by neighbourhood until the partition is stable. Any
/// tie that remains between equivalent atoms is broken by trying each
/// candidate and keeping the lexicographically smallest string. The result
/// depends only on the molecular graph, so the same molecule gives the same
/// string however its atoms are numbered, and it reads back through
/// [`parse_smiles`] to the same composition and bond orders.
///
/// Hydrogens bonded to a heavy atom are folded into that atom's count. Atoms
/// use the organic-subset form when their hydrogen count matches the default
/// valence, and brackets otherwise. Aromatic atoms are written in lowercase.
///
/// # Errors
/// [`IoError::Unsupported`] when the molecule has a net formal charge
/// (per-atom charges are not stored in [`Molecule`]), an element has no
/// symbol, or more than 99 ring closures are open at once.
pub fn write_smiles(mol: &Molecule) -> Result<String> {
    if mol.formal_charge() != 0 {
        return Err(IoError::Unsupported {
            format: FORMAT,
            message: "net formal charge: per-atom charges are not stored in Molecule".into(),
        });
    }
    let graph = Graph::from_molecule(mol)?;
    let ranks = refine(&graph, graph.initial_ranks());
    canonical_string(&graph, ranks)
}

/// Heavy-atom graph with hydrogens folded into counts, plus per-atom text.
struct Graph {
    z: Vec<u8>,
    folded: Vec<bool>,
    h: Vec<u32>,
    aromatic: Vec<bool>,
    text: Vec<String>,
    adj: Vec<Vec<(usize, BondOrder)>>,
    heavy: Vec<usize>,
}

impl Graph {
    fn from_molecule(mol: &Molecule) -> Result<Graph> {
        let z: Vec<u8> = mol.atoms().map(|a| a.z).collect();
        let n = z.len();

        // A hydrogen bonded to exactly one non-hydrogen is folded into its neighbour.
        let folded: Vec<bool> = (0..n)
            .map(|i| {
                let id = AtomId::from(i);
                z[i] == 1
                    && mol.degree(id) == 1
                    && mol
                        .neighbors(id)
                        .next()
                        .is_some_and(|nb| z[nb.index()] != 1)
            })
            .collect();
        let mut h = vec![0u32; n];
        for i in (0..n).filter(|&i| folded[i]) {
            let heavy = mol
                .neighbors(AtomId::from(i))
                .next()
                .map_or(0, |nb| nb.index());
            h[heavy] += 1;
        }

        let mut adj: Vec<Vec<(usize, BondOrder)>> = vec![Vec::new(); n];
        for bond in mol.bonds() {
            let (a, b) = (bond.a.index(), bond.b.index());
            if folded[a] || folded[b] {
                continue;
            }
            adj[a].push((b, bond.order));
            adj[b].push((a, bond.order));
        }
        for list in &mut adj {
            list.sort_by_key(|&(j, _)| j);
        }

        let aromatic: Vec<bool> = adj
            .iter()
            .map(|list| list.iter().any(|&(_, o)| o == BondOrder::Aromatic))
            .collect();
        let mut text = Vec::with_capacity(n);
        for i in 0..n {
            text.push(if folded[i] {
                String::new()
            } else {
                atom_text(z[i], h[i], aromatic[i], &adj[i])?
            });
        }
        let heavy = (0..n).filter(|&i| !folded[i]).collect();
        Ok(Graph {
            z,
            folded,
            h,
            aromatic,
            text,
            adj,
            heavy,
        })
    }

    /// Dense ranks from element, degree, hydrogen count and aromaticity.
    fn initial_ranks(&self) -> Vec<usize> {
        let keys: Vec<(u8, usize, u32, bool)> = self
            .heavy
            .iter()
            .map(|&i| (self.z[i], self.adj[i].len(), self.h[i], self.aromatic[i]))
            .collect();
        self.dense(&keys)
    }

    /// Dense ranks for the heavy atoms, from one key per heavy atom (in
    /// `self.heavy` order). Non-heavy slots get `usize::MAX`.
    fn dense<K: Ord + Clone>(&self, keys: &[K]) -> Vec<usize> {
        let mut unique: Vec<K> = keys.to_vec();
        unique.sort();
        unique.dedup();
        let mut ranks = vec![usize::MAX; self.z.len()];
        for (k, &i) in self.heavy.iter().enumerate() {
            ranks[i] = unique
                .binary_search(&keys[k])
                .expect("every key appears in the unique list");
        }
        ranks
    }
}

/// Numeric code for a bond order, used in refinement keys.
fn bond_code(order: BondOrder) -> u8 {
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 4,
    }
}

/// Refine a partition until no class splits further.
fn refine(graph: &Graph, mut ranks: Vec<usize>) -> Vec<usize> {
    let classes = |r: &[usize]| {
        let mut seen: Vec<usize> = graph.heavy.iter().map(|&i| r[i]).collect();
        seen.sort();
        seen.dedup();
        seen.len()
    };
    loop {
        let keys: Vec<(usize, Vec<(u8, usize)>)> = graph
            .heavy
            .iter()
            .map(|&i| {
                let mut neighbours: Vec<(u8, usize)> = graph.adj[i]
                    .iter()
                    .map(|&(j, o)| (bond_code(o), ranks[j]))
                    .collect();
                neighbours.sort();
                (ranks[i], neighbours)
            })
            .collect();
        let next = graph.dense(&keys);
        if classes(&next) == classes(&ranks) {
            return next;
        }
        ranks = next;
    }
}

/// Smallest string over all tie-breaks of a refined partition.
fn canonical_string(graph: &Graph, ranks: Vec<usize>) -> Result<String> {
    let mut class_size = vec![0usize; graph.z.len() + 1];
    for &i in &graph.heavy {
        class_size[ranks[i]] += 1;
    }
    let tied = (0..class_size.len()).find(|&r| class_size[r] > 1);
    let Some(tied) = tied else {
        return write_with_ranks(graph, &ranks);
    };

    let mut best: Option<String> = None;
    for &candidate in graph.heavy.iter().filter(|&&i| ranks[i] == tied) {
        // Individualise one atom: it sorts before the rest of its class.
        let keys: Vec<(usize, usize)> = graph
            .heavy
            .iter()
            .map(|&i| (ranks[i], usize::from(i != candidate)))
            .collect();
        let split = refine(graph, graph.dense(&keys));
        let text = canonical_string(graph, split)?;
        if best.as_ref().is_none_or(|b| text < *b) {
            best = Some(text);
        }
    }
    Ok(best.unwrap_or_default())
}

/// Emit SMILES with neighbours and roots visited in rank order.
fn write_with_ranks(graph: &Graph, ranks: &[usize]) -> Result<String> {
    let n = graph.z.len();
    let mut adj = graph.adj.clone();
    for list in &mut adj {
        list.sort_by_key(|&(j, _)| ranks[j]);
    }
    let mut roots = graph.heavy.clone();
    roots.sort_by_key(|&i| ranks[i]);

    let mut walk = Walk {
        children: vec![Vec::new(); n],
        closures: vec![Vec::new(); n],
        seen: graph.folded.clone(),
        recorded: Default::default(),
    };
    let mut starts = Vec::new();
    for &root in &roots {
        if !walk.seen[root] {
            walk.seen[root] = true;
            starts.push(root);
            walk.visit(&adj, root, None);
        }
    }

    let mut emit = Emit {
        text: &graph.text,
        adj: &adj,
        aromatic: &graph.aromatic,
        children: &walk.children,
        closures: &walk.closures,
        emitted: vec![false; n],
        digits: BTreeMap::new(),
        out: String::new(),
    };
    for (k, &root) in starts.iter().enumerate() {
        if k > 0 {
            emit.out.push('.');
        }
        emit.atom(root)?;
    }
    Ok(emit.out)
}

/// Text for one heavy atom: organic-subset symbol when the hydrogen count
/// fits the default valence, otherwise a bracket atom.
fn atom_text(z: u8, h: u32, aromatic: bool, bonds: &[(usize, BondOrder)]) -> Result<String> {
    let symbol = element::symbol(z).ok_or_else(|| IoError::Unsupported {
        format: FORMAT,
        message: format!("atomic number {z} has no element symbol"),
    })?;
    let written = |s: &str| {
        if aromatic {
            s.to_ascii_lowercase()
        } else {
            s.to_string()
        }
    };
    if let Some(valence) = default_valence(z) {
        // Aromatic bonds count 1.5; floor() gives the integer bond total.
        let sum: f64 = bonds.iter().map(|&(_, o)| o.value()).sum();
        let implicit = (valence - sum.floor() as i32).max(0) as u32;
        let aromatic_ok = !aromatic || matches!(z, 5 | 6 | 7 | 8 | 15 | 16);
        if implicit == h && aromatic_ok {
            return Ok(written(symbol));
        }
    }
    let mut out = String::from("[");
    out.push_str(&written(symbol));
    if h == 1 {
        out.push('H');
    } else if h > 1 {
        out.push_str(&format!("H{h}"));
    }
    out.push(']');
    Ok(out)
}

/// Bond a SMILES parser infers when no symbol is written.
fn default_bond(a_aromatic: bool, b_aromatic: bool) -> BondOrder {
    if a_aromatic && b_aromatic {
        BondOrder::Aromatic
    } else {
        BondOrder::Single
    }
}

/// Bond symbol to write, or nothing when the parser would infer the same order.
fn bond_symbol(order: BondOrder, default: BondOrder) -> &'static str {
    if order == default {
        return "";
    }
    match order {
        BondOrder::Single => "-",
        BondOrder::Double => "=",
        BondOrder::Triple => "#",
        BondOrder::Aromatic => ":",
    }
}

fn edge_key(a: usize, b: usize) -> (usize, usize) {
    (a.min(b), a.max(b))
}

/// Spanning-tree walk: tree children and ring-closure partners per atom.
struct Walk {
    children: Vec<Vec<usize>>,
    closures: Vec<Vec<usize>>,
    seen: Vec<bool>,
    recorded: std::collections::HashSet<(usize, usize)>,
}

impl Walk {
    fn visit(&mut self, adj: &[Vec<(usize, BondOrder)>], u: usize, parent: Option<usize>) {
        for &(v, _) in &adj[u] {
            if Some(v) == parent {
                continue;
            }
            if !self.seen[v] {
                self.seen[v] = true;
                self.children[u].push(v);
                self.visit(adj, v, Some(u));
            } else if self.recorded.insert(edge_key(u, v)) {
                // A non-tree edge: the first endpoint to see it closes a ring.
                self.closures[u].push(v);
                self.closures[v].push(u);
            }
        }
    }
}

/// Emission state for pass 2.
struct Emit<'a> {
    text: &'a [String],
    adj: &'a [Vec<(usize, BondOrder)>],
    aromatic: &'a [bool],
    children: &'a [Vec<usize>],
    closures: &'a [Vec<usize>],
    emitted: Vec<bool>,
    /// Ring digit assigned to each open ring closure.
    digits: BTreeMap<(usize, usize), u32>,
    out: String,
}

impl Emit<'_> {
    fn order(&self, a: usize, b: usize) -> BondOrder {
        self.adj[a]
            .iter()
            .find(|&&(j, _)| j == b)
            .map_or(BondOrder::Single, |&(_, o)| o)
    }

    fn default_between(&self, a: usize, b: usize) -> BondOrder {
        default_bond(self.aromatic[a], self.aromatic[b])
    }

    fn push_digit(&mut self, d: u32) -> Result<()> {
        if d > 99 {
            return Err(IoError::Unsupported {
                format: FORMAT,
                message: "more than 99 ring closures open at once".into(),
            });
        }
        if d < 10 {
            self.out.push(char::from(b'0' + d as u8));
        } else {
            self.out.push_str(&format!("%{d:02}"));
        }
        Ok(())
    }

    fn atom(&mut self, u: usize) -> Result<()> {
        self.out.push_str(&self.text[u]);
        self.emitted[u] = true;
        // Copy the shared slices out so the loops can call `&mut self` methods.
        let closures = self.closures;
        let children = self.children;

        for &p in &closures[u] {
            let key = edge_key(u, p);
            if self.emitted[p] {
                // Closing: the partner opened this ring earlier.
                let d = self.digits.remove(&key).unwrap_or(0);
                self.push_digit(d)?;
            } else {
                // Opening: smallest digit not already in use.
                let d = (1..)
                    .find(|d| !self.digits.values().any(|used| used == d))
                    .unwrap_or(1);
                self.digits.insert(key, d);
                let sym = bond_symbol(self.order(u, p), self.default_between(u, p));
                self.out.push_str(sym);
                self.push_digit(d)?;
            }
        }

        for (k, &c) in children[u].iter().enumerate() {
            let sym = bond_symbol(self.order(u, c), self.default_between(u, c));
            let last = k + 1 == children[u].len();
            if last {
                self.out.push_str(sym);
                self.atom(c)?;
            } else {
                self.out.push('(');
                self.out.push_str(sym);
                self.atom(c)?;
                self.out.push(')');
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn organic_subset_with_implicit_hydrogens() {
        assert_eq!(parse_smiles("C").unwrap().formula(), "CH4");
        assert_eq!(parse_smiles("CCO").unwrap().formula(), "C2H6O");
        assert_eq!(parse_smiles("O=C=O").unwrap().formula(), "CO2");
        assert_eq!(parse_smiles("C#N").unwrap().formula(), "CHN");
        assert_eq!(parse_smiles("CC(=O)O").unwrap().formula(), "C2H4O2");
        assert_eq!(parse_smiles("ClCBr").unwrap().formula(), "CBrClH2");
    }

    #[test]
    fn rings_and_aromaticity() {
        assert_eq!(parse_smiles("C1CCCCC1").unwrap().formula(), "C6H12");
        assert_eq!(parse_smiles("c1ccccc1").unwrap().formula(), "C6H6");
        assert_eq!(parse_smiles("c1ccncc1").unwrap().formula(), "C5H5N");
        // Thiophene: aromatic sulfur has no implicit hydrogen.
        assert_eq!(parse_smiles("c1ccsc1").unwrap().formula(), "C4H4S");
        // Percent-style ring closure: cyclopropane.
        assert_eq!(parse_smiles("C%12CC%12").unwrap().formula(), "C3H6");
    }

    #[test]
    fn aromatic_bonds_are_recorded() {
        // 6 aromatic C–C bonds plus 6 explicit C–H bonds.
        let benzene = parse_smiles("c1ccccc1").unwrap();
        assert_eq!(benzene.n_bonds(), 12);
        let aromatic = benzene
            .bonds()
            .filter(|b| b.order == BondOrder::Aromatic)
            .count();
        assert_eq!(aromatic, 6);
    }

    #[test]
    fn bracket_atoms_carry_charge_and_hydrogens() {
        let ammonium = parse_smiles("[NH4+]").unwrap();
        assert_eq!(ammonium.formula(), "H4N");
        assert_eq!(ammonium.formal_charge(), 1);

        let salt = parse_smiles("[Na+].[Cl-]").unwrap();
        assert_eq!(salt.formula(), "ClNa");
        assert_eq!(salt.formal_charge(), 0);

        let pyrrole = parse_smiles("c1cc[nH]c1").unwrap();
        assert_eq!(pyrrole.formula(), "C4H5N");
    }

    #[test]
    fn double_charge_and_stereo_marks_are_read() {
        let dication = parse_smiles("[Ca++]").unwrap();
        assert_eq!(dication.formal_charge(), 2);
        // Chirality and isotope are accepted and ignored.
        assert_eq!(parse_smiles("[13CH4]").unwrap().formula(), "CH4");
        assert_eq!(parse_smiles("F[C@H](Cl)Br").unwrap().formula(), "CBrClFH");
    }

    #[test]
    fn trailing_name_is_ignored() {
        assert_eq!(parse_smiles("CCO ethanol").unwrap().formula(), "C2H6O");
    }

    fn bond_order_counts(mol: &Molecule) -> [usize; 4] {
        let mut counts = [0usize; 4];
        for bond in mol.bonds() {
            let slot = match bond.order {
                BondOrder::Single => 0,
                BondOrder::Double => 1,
                BondOrder::Triple => 2,
                BondOrder::Aromatic => 3,
            };
            counts[slot] += 1;
        }
        counts
    }

    fn assert_round_trip(smiles: &str) {
        let first = parse_smiles(smiles).unwrap();
        let written = write_smiles(&first).unwrap();
        let second =
            parse_smiles(&written).unwrap_or_else(|e| panic!("{smiles} -> {written}: {e}"));
        assert_eq!(first.formula(), second.formula(), "{smiles} -> {written}");
        assert_eq!(first.n_bonds(), second.n_bonds(), "{smiles} -> {written}");
        assert_eq!(
            bond_order_counts(&first),
            bond_order_counts(&second),
            "{smiles} -> {written}"
        );
    }

    #[test]
    fn writes_simple_molecules_exactly() {
        assert_eq!(write_smiles(&parse_smiles("CCO").unwrap()).unwrap(), "CCO");
        assert_eq!(
            write_smiles(&parse_smiles("CC(=O)O").unwrap()).unwrap(),
            "CC(=O)O"
        );
        assert_eq!(
            write_smiles(&parse_smiles("c1ccccc1").unwrap()).unwrap(),
            "c1ccccc1"
        );
        assert_eq!(write_smiles(&parse_smiles("C#N").unwrap()).unwrap(), "C#N");
        assert_eq!(
            write_smiles(&parse_smiles("[H][H]").unwrap()).unwrap(),
            "[H][H]"
        );
    }

    #[test]
    fn round_trips_rings_branches_and_fused_systems() {
        for smiles in [
            "C1CCCCC1",
            "c1ccncc1",
            "c1ccsc1",
            "c1cc[nH]c1",
            "C%12CC%12",
            "c1ccc2ccccc2c1",
            "C1CC2CCC1C2",
            "CC(C)(C)c1ccccc1",
            "Cc1ccccc1-c1ccccc1",
            "O=C=O",
            "ClCBr",
            "C(=O)(O)N",
        ] {
            assert_round_trip(smiles);
        }
    }

    /// Same molecule with atoms added in the order `order` (new atom k is old atom `order[k]`).
    fn relabelled(mol: &Molecule, order: &[usize]) -> Molecule {
        let mut new_index = vec![0usize; order.len()];
        for (k, &old) in order.iter().enumerate() {
            new_index[old] = k;
        }
        let mut out = Molecule::new("");
        for &old in order {
            let atom = mol.atom(old.into());
            out.add_atom_raw(atom.z, atom.pos);
        }
        for bond in mol.bonds() {
            let a = new_index[bond.a.index()];
            let b = new_index[bond.b.index()];
            out.add_bond(a.into(), b.into(), bond.order).unwrap();
        }
        out
    }

    #[test]
    fn canonical_output_ignores_atom_numbering() {
        for smiles in [
            "CCO",
            "CC(=O)O",
            "OC(C)=O",
            "c1ccccc1",
            "c1ccncc1",
            "CC(C)(C)c1ccccc1",
            "C1CC2CCC1C2",
            "c1ccc2ccccc2c1",
            "CC(C)CC(C)C",
            "[H][H]",
            "CC.O",
        ] {
            let mol = parse_smiles(smiles).unwrap();
            let n = mol.len();
            let reference = write_smiles(&mol).unwrap();
            let reversed: Vec<usize> = (0..n).rev().collect();
            let rotated: Vec<usize> = (0..n).map(|k| (k + 3) % n).collect();
            for order in [reversed, rotated] {
                let other = relabelled(&mol, &order);
                assert_eq!(
                    write_smiles(&other).unwrap(),
                    reference,
                    "{smiles} relabelled as {order:?}"
                );
            }
        }
    }

    #[test]
    fn equivalent_spellings_share_one_canonical_form() {
        let a = write_smiles(&parse_smiles("CC(=O)O").unwrap()).unwrap();
        let b = write_smiles(&parse_smiles("OC(C)=O").unwrap()).unwrap();
        let c = write_smiles(&parse_smiles("O=C(O)C").unwrap()).unwrap();
        assert_eq!(a, b);
        assert_eq!(a, c);
    }

    #[test]
    fn net_charge_is_refused() {
        let ammonium = parse_smiles("[NH4+]").unwrap();
        assert!(matches!(
            write_smiles(&ammonium),
            Err(IoError::Unsupported { .. })
        ));
    }

    #[test]
    fn malformed_input_is_rejected() {
        assert!(parse_smiles("").is_err());
        assert!(parse_smiles("C(C").is_err());
        assert!(parse_smiles("C)").is_err());
        assert!(parse_smiles("C1CC").is_err());
        assert!(parse_smiles("[Xx]").is_err());
        assert!(parse_smiles("[NH4+").is_err());
        assert!(parse_smiles("CX").is_err());
    }
}
