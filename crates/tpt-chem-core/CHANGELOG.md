# Changelog

All notable changes to `tpt-chem-core` are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - not yet released

### Added

- Foundation modules: elements, unit-safe types, molecular graph (`molecule`), force-field parameter types (`forcefield`), Cartesian and internal coordinates (`coords`), and a seedable RNG (`rng`).
- `no_std + alloc` support: disabling the default `std` feature builds the crate for bare-metal targets, with transcendental math via `libm`.
- Optional `serde` feature for serializing the core types.

### Changed

- Workspace-wide platform review: CI sibling-checkout setup, and correctness fixes across the crates.
- Analytic gradients, MD exclusions and constraints, and UHF support landed across the workspace, including the core types they use.

### Fixed

- Rustdoc link warnings; the workspace documents with zero warnings.
