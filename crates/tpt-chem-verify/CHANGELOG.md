# Changelog

All notable changes to `tpt-chem-verify` are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - not yet released

### Added

- `proptest` strategies for molecular geometries and force-field parameters.
- Energy and symplecticity conservation-law checks for MD.
- Mass-conservation properties for closed reaction networks (RK4, stiff solver, SSA, randomized states).
- Kani bounded-model-checking harnesses behind `cfg(kani)`: cell-index bounds, Ewald index range, and Verlet skin monotonicity. CI runs them on Linux.
- Property tests for the aqueous models.
