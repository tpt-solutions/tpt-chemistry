//! Space-group tables: all 230 standard settings against group axioms and
//! known operation counts, and the `by_symbol` lookups against the table.

use tpt_chem_crystal::symmetry::{compose, SpaceGroup, SymmetryOp};

/// Translations are equal modulo one lattice vector.
fn same_mod_one(a: [f64; 3], b: [f64; 3]) -> bool {
    a.iter().zip(&b).all(|(x, y)| {
        let d = (x - y).rem_euclid(1.0);
        d < 1e-9 || (1.0 - d) < 1e-9
    })
}

fn same_op(a: &SymmetryOp, b: &SymmetryOp) -> bool {
    a.r == b.r && same_mod_one(a.t, b.t)
}

fn contains(ops: &[SymmetryOp], op: &SymmetryOp) -> bool {
    ops.iter().any(|o| same_op(o, op))
}

/// Orbit of a generic position as a sorted list of canonical points.
fn generic_orbit(group: &SpaceGroup) -> Vec<[i64; 3]> {
    let seed = [0.1137, 0.2291, 0.3733];
    let mut points: Vec<[i64; 3]> = group
        .operations
        .iter()
        .map(|op| {
            let p = op.apply(seed);
            // Canonical integer key: 1e6 resolution after reducing mod 1.
            let key = |x: f64| ((x.rem_euclid(1.0) * 1e6).round() as i64) % 1_000_000;
            [key(p[0]), key(p[1]), key(p[2])]
        })
        .collect();
    points.sort();
    points.dedup();
    points
}

#[test]
fn every_number_from_1_to_230_has_a_group() {
    for n in 1..=230u16 {
        let group = SpaceGroup::from_number(n).unwrap_or_else(|| panic!("no group {n}"));
        assert_eq!(group.number, Some(n));
        assert!(!group.symbol.is_empty(), "group {n} has no symbol");
        assert!(
            contains(&group.operations, &SymmetryOp::identity()),
            "group {n} lacks the identity"
        );
    }
}

#[test]
fn out_of_range_numbers_are_none() {
    assert!(SpaceGroup::from_number(0).is_none());
    assert!(SpaceGroup::from_number(231).is_none());
}

#[test]
fn every_group_is_closed_under_composition() {
    for n in 1..=230u16 {
        let group = SpaceGroup::from_number(n).unwrap();
        for a in &group.operations {
            for b in &group.operations {
                let product = compose(a, b);
                assert!(
                    contains(&group.operations, &product),
                    "group {n} ({}) not closed: {a:?} ∘ {b:?}",
                    group.symbol
                );
            }
        }
    }
}

#[test]
fn known_operation_counts() {
    let expected = [
        (1, 1),
        (2, 2),
        (4, 2),
        (14, 4),
        (19, 4),
        (33, 4),
        (75, 4),
        (147, 6),
        (191, 24),
        (206, 48),
        (225, 192),
        (230, 96),
    ];
    for (n, count) in expected {
        let group = SpaceGroup::from_number(n).unwrap();
        assert_eq!(
            group.operations.len(),
            count,
            "group {n} ({}) operation count",
            group.symbol
        );
    }
}

#[test]
fn by_symbol_groups_match_the_table() {
    for symbol in [
        "P1", "P-1", "P21", "P21/c", "P212121", "Pna21", "P4", "P-3", "P6/mmm", "Fm-3m",
    ] {
        let coded = SpaceGroup::by_symbol(symbol).unwrap();
        let number = coded.number.expect("by_symbol groups carry a number");
        let table = SpaceGroup::from_number(number).unwrap();
        assert_eq!(
            coded.operations.len(),
            table.operations.len(),
            "{symbol}: operation count differs from table group {number} ({})",
            table.symbol
        );
        assert_eq!(
            generic_orbit(&coded),
            generic_orbit(&table),
            "{symbol}: generic orbit differs from table group {number} ({})",
            table.symbol
        );
    }
}
