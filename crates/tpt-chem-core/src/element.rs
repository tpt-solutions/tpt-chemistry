//! The periodic table with const-generic atom types.
//!
//! [`Element`] is a zero-sized type parameterized by the atomic number at
//! compile time: `Element<6>` is carbon, and `Element::<6>::symbol()` is a
//! `const fn` reading the static table. An invalid `Z` fails to compile —
//! spec.txt §4's "const-generic atom types".
//!
//! ```
//! use tpt_chem_core::element::Element;
//!
//! let z_carbon: u8 = Element::<6>::Z;
//! assert_eq!(Element::<6>::symbol(), "C");
//! assert_eq!(z_carbon, 6);
//! assert!((Element::<1>::mass() - 1.008).abs() < 1e-12);
//! ```
//!
//! Dynamic code (parsers, MD systems) uses the runtime accessors
//! [`symbol`], [`mass`], [`from_symbol`], … over the same table.

use alloc::format;
use alloc::string::String;
use alloc::vec::Vec;

/// Maximum supported atomic number (oganesson).
pub const MAX_Z: u8 = 118;

/// One row of the static periodic table:
/// `(symbol, name, standard atomic mass [amu], covalent radius [Å])`.
///
/// Masses are IUPAC standard atomic weights (conventional values for
/// interval elements, exact monoisotopic masses otherwise). Covalent radii
/// follow Cordero et al. (2008); values for the superheavy tail (Z ≥ 97)
/// are convention estimates.
const TABLE: [(&str, &str, f64, f64); MAX_Z as usize] = [
    ("H", "hydrogen", 1.008, 0.31),
    ("He", "helium", 4.002602, 0.28),
    ("Li", "lithium", 6.94, 1.28),
    ("Be", "beryllium", 9.0121831, 0.96),
    ("B", "boron", 10.81, 0.84),
    ("C", "carbon", 12.011, 0.76),
    ("N", "nitrogen", 14.007, 0.71),
    ("O", "oxygen", 15.999, 0.66),
    ("F", "fluorine", 18.998403163, 0.57),
    ("Ne", "neon", 20.1797, 0.58),
    ("Na", "sodium", 22.98976928, 1.66),
    ("Mg", "magnesium", 24.305, 1.41),
    ("Al", "aluminium", 26.9815385, 1.21),
    ("Si", "silicon", 28.085, 1.11),
    ("P", "phosphorus", 30.973761998, 1.07),
    ("S", "sulfur", 32.06, 1.05),
    ("Cl", "chlorine", 35.45, 1.02),
    ("Ar", "argon", 39.948, 1.06),
    ("K", "potassium", 39.0983, 2.03),
    ("Ca", "calcium", 40.078, 1.76),
    ("Sc", "scandium", 44.955908, 1.70),
    ("Ti", "titanium", 47.867, 1.60),
    ("V", "vanadium", 50.9415, 1.53),
    ("Cr", "chromium", 51.9961, 1.39),
    ("Mn", "manganese", 54.938044, 1.39),
    ("Fe", "iron", 55.845, 1.32),
    ("Co", "cobalt", 58.933194, 1.26),
    ("Ni", "nickel", 58.6934, 1.24),
    ("Cu", "copper", 63.546, 1.32),
    ("Zn", "zinc", 65.38, 1.22),
    ("Ga", "gallium", 69.723, 1.22),
    ("Ge", "germanium", 72.630, 1.20),
    ("As", "arsenic", 74.921595, 1.19),
    ("Se", "selenium", 78.971, 1.20),
    ("Br", "bromine", 79.904, 1.20),
    ("Kr", "krypton", 83.798, 1.16),
    ("Rb", "rubidium", 85.4678, 2.20),
    ("Sr", "strontium", 87.62, 1.95),
    ("Y", "yttrium", 88.90584, 1.90),
    ("Zr", "zirconium", 91.224, 1.75),
    ("Nb", "niobium", 92.90637, 1.64),
    ("Mo", "molybdenum", 95.95, 1.54),
    ("Tc", "technetium", 98.0, 1.47),
    ("Ru", "ruthenium", 101.07, 1.46),
    ("Rh", "rhodium", 102.90550, 1.42),
    ("Pd", "palladium", 106.42, 1.39),
    ("Ag", "silver", 107.8682, 1.45),
    ("Cd", "cadmium", 112.414, 1.44),
    ("In", "indium", 114.818, 1.42),
    ("Sn", "tin", 118.710, 1.39),
    ("Sb", "antimony", 121.760, 1.39),
    ("Te", "tellurium", 127.60, 1.38),
    ("I", "iodine", 126.90447, 1.39),
    ("Xe", "xenon", 131.293, 1.40),
    ("Cs", "caesium", 132.90545196, 2.44),
    ("Ba", "barium", 137.327, 2.15),
    ("La", "lanthanum", 138.90547, 2.07),
    ("Ce", "cerium", 140.116, 2.04),
    ("Pr", "praseodymium", 140.90766, 2.03),
    ("Nd", "neodymium", 144.242, 2.01),
    ("Pm", "promethium", 145.0, 1.99),
    ("Sm", "samarium", 150.36, 1.98),
    ("Eu", "europium", 151.964, 1.98),
    ("Gd", "gadolinium", 157.25, 1.96),
    ("Tb", "terbium", 158.92535, 1.94),
    ("Dy", "dysprosium", 162.500, 1.92),
    ("Ho", "holmium", 164.93033, 1.92),
    ("Er", "erbium", 167.259, 1.89),
    ("Tm", "thulium", 168.93422, 1.90),
    ("Yb", "ytterbium", 173.045, 1.87),
    ("Lu", "lutetium", 174.9668, 1.87),
    ("Hf", "hafnium", 178.49, 1.75),
    ("Ta", "tantalum", 180.94788, 1.70),
    ("W", "tungsten", 183.84, 1.62),
    ("Re", "rhenium", 186.207, 1.51),
    ("Os", "osmium", 190.23, 1.44),
    ("Ir", "iridium", 192.217, 1.41),
    ("Pt", "platinum", 195.084, 1.36),
    ("Au", "gold", 196.966569, 1.36),
    ("Hg", "mercury", 200.592, 1.32),
    ("Tl", "thallium", 204.38, 1.45),
    ("Pb", "lead", 207.2, 1.46),
    ("Bi", "bismuth", 208.98040, 1.48),
    ("Po", "polonium", 209.0, 1.40),
    ("At", "astatine", 210.0, 1.50),
    ("Rn", "radon", 222.0, 1.50),
    ("Fr", "francium", 223.0, 2.60),
    ("Ra", "radium", 226.0, 2.21),
    ("Ac", "actinium", 227.0, 2.15),
    ("Th", "thorium", 232.0377, 2.06),
    ("Pa", "protactinium", 231.03588, 2.00),
    ("U", "uranium", 238.02891, 1.96),
    ("Np", "neptunium", 237.0, 1.90),
    ("Pu", "plutonium", 244.0, 1.87),
    ("Am", "americium", 243.0, 1.80),
    ("Cm", "curium", 247.0, 1.69),
    ("Bk", "berkelium", 247.0, 1.66),
    ("Cf", "californium", 251.0, 1.68),
    ("Es", "einsteinium", 252.0, 1.65),
    ("Fm", "fermium", 257.0, 1.67),
    ("Md", "mendelevium", 258.0, 1.73),
    ("No", "nobelium", 259.0, 1.76),
    ("Lr", "lawrencium", 266.0, 1.61),
    ("Rf", "rutherfordium", 267.0, 1.60),
    ("Db", "dubnium", 268.0, 1.61),
    ("Sg", "seaborgium", 269.0, 1.62),
    ("Bh", "bohrium", 270.0, 1.63),
    ("Hs", "hassium", 269.0, 1.65),
    ("Mt", "meitnerium", 278.0, 1.50),
    ("Ds", "darmstadtium", 281.0, 1.50),
    ("Rg", "roentgenium", 282.0, 1.50),
    ("Cn", "copernicium", 285.0, 1.50),
    ("Nh", "nihonium", 286.0, 1.50),
    ("Fl", "flerovium", 289.0, 1.50),
    ("Mc", "moscovium", 290.0, 1.50),
    ("Lv", "livermorium", 293.0, 1.50),
    ("Ts", "tennessine", 294.0, 1.50),
    ("Og", "oganesson", 294.0, 1.50),
];

/// Compile-time-checked chemical element `Z`.
///
/// The type is zero-sized; `Z` is carried entirely in the type parameter and
/// validated at compile time — `Element<0>` and `Element<119>` do not
/// compile.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Element<const Z: u8>;

impl<const Z: u8> Element<Z> {
    /// Compile-time validity gate. Any use in a `const` context of an
    /// out-of-range `Z` (e.g. `Element::<0>::symbol()`) fails to compile;
    /// runtime misuse panics via the same check.
    const VALID: () = assert!(
        (Z as usize) >= 1 && (Z as usize) <= MAX_Z as usize,
        "atomic number out of range (1..=118)"
    );

    /// The atomic number as a runtime value.
    pub const Z: u8 = Z;

    /// Construct the element marker (exists so callers can annotate intent).
    pub const fn new() -> Self {
        let () = Self::VALID;
        Element::<Z>
    }

    /// Element symbol, e.g. `"C"`.
    pub const fn symbol() -> &'static str {
        let () = Self::VALID;
        TABLE[(Z as usize) - 1].0
    }

    /// Element name, e.g. `"carbon"`.
    pub const fn name() -> &'static str {
        let () = Self::VALID;
        TABLE[(Z as usize) - 1].1
    }

    /// Standard atomic mass in amu.
    pub const fn mass() -> f64 {
        let () = Self::VALID;
        TABLE[(Z as usize) - 1].2
    }

    /// Covalent radius in Å (Cordero et al. 2008).
    pub const fn covalent_radius() -> f64 {
        let () = Self::VALID;
        TABLE[(Z as usize) - 1].3
    }
}

/// Symbol of element `z`, or `None` out of range.
pub fn symbol(z: u8) -> Option<&'static str> {
    if (1..=MAX_Z).contains(&z) {
        Some(TABLE[(z - 1) as usize].0)
    } else {
        None
    }
}

/// Name of element `z`, or `None` out of range.
pub fn name(z: u8) -> Option<&'static str> {
    if (1..=MAX_Z).contains(&z) {
        Some(TABLE[(z - 1) as usize].1)
    } else {
        None
    }
}

/// Standard atomic mass [amu] of element `z`, or `None` out of range.
pub fn mass(z: u8) -> Option<f64> {
    if (1..=MAX_Z).contains(&z) {
        Some(TABLE[(z - 1) as usize].2)
    } else {
        None
    }
}

/// Covalent radius [Å] of element `z`, or `None` out of range.
pub fn covalent_radius(z: u8) -> Option<f64> {
    if (1..=MAX_Z).contains(&z) {
        Some(TABLE[(z - 1) as usize].3)
    } else {
        None
    }
}

/// Atomic number for a chemical symbol (case-insensitive), or `None`.
pub fn from_symbol(s: &str) -> Option<u8> {
    let needle = s.trim();
    if needle.is_empty() || needle.len() > 3 {
        return None;
    }
    TABLE
        .iter()
        .position(|(sym, _, _, _)| sym.eq_ignore_ascii_case(needle))
        .map(|i| (i + 1) as u8)
}

/// Hill-order molecular formula: C and H first (if present), then everything
/// else alphabetically. E.g. `"H2O"`, `"C2H6O"`.
pub fn hill_formula(counts: &[u8; MAX_Z as usize]) -> String {
    let mut out = Vec::new();
    let carbon = counts[5];
    let hydrogen = counts[0];
    if carbon > 0 {
        out.push(("C", carbon));
    }
    if hydrogen > 0 {
        out.push(("H", hydrogen));
    }
    for (i, &n) in counts.iter().enumerate() {
        let z = (i + 1) as u8;
        if n > 0 && z != 6 && z != 1 {
            out.push((symbol(z).unwrap_or("?"), n));
        }
    }
    // Hill order: C and H first, the rest alphabetical by symbol.
    out[if carbon > 0 { 1 } else { 0 }..].sort_by(|a, b| a.0.cmp(b.0));
    let mut s = String::new();
    for (sym, n) in out {
        s.push_str(sym);
        if n > 1 {
            s.push_str(&format!("{n}"));
        }
    }
    if s.is_empty() {
        s.push('∅');
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compile_time_lookup() {
        assert_eq!(Element::<6>::symbol(), "C");
        assert_eq!(Element::<6>::name(), "carbon");
        assert_eq!(Element::<79>::symbol(), "Au");
        assert_eq!(Element::<118>::symbol(), "Og");
        assert!((Element::<8>::mass() - 15.999).abs() < 1e-12);
        assert_eq!(Element::<6>::Z, 6);
    }

    #[test]
    fn runtime_lookup() {
        assert_eq!(symbol(1), Some("H"));
        assert_eq!(symbol(0), None);
        assert_eq!(symbol(119), None);
        assert_eq!(from_symbol("cl"), Some(17));
        assert_eq!(from_symbol("Fe"), Some(26));
        assert_eq!(from_symbol("Xx"), None);
        assert!((mass(26).unwrap() - 55.845).abs() < 1e-12);
        assert!(covalent_radius(1).unwrap() < covalent_radius(6).unwrap());
    }

    #[test]
    fn hill_formula_counts() {
        let mut counts = [0u8; MAX_Z as usize];
        counts[5] = 2; // C
        counts[0] = 6; // H
        counts[6] = 1; // N
        assert_eq!(hill_formula(&counts), "C2H6N");
        let mut water = [0u8; MAX_Z as usize];
        water[0] = 2;
        water[7] = 1; // O
        assert_eq!(hill_formula(&water), "H2O");
    }
}
