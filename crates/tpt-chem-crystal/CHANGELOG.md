# Changelog

All notable changes to `tpt-chem-crystal` are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - not yet released

### Added

- Bravais lattices, Miller indices, and Wigner–Seitz cell construction.
- Space-group symmetry from P1 through Fm-3m, with Wyckoff orbit generation.
- Reciprocal-lattice k-point grids.
- Powder XRD simulation with Cromer–Mann form factors, the Lorentz–polarization factor, and extinction rules. Silicon and NaCl are validated against reference patterns.
- Optional `serde` feature for result types.

### Fixed

- Unwrap audit: floating-point sorts use `total_cmp`, and the reciprocal-lattice panics are documented.
