# Changelog

All notable changes to `tpt-chem-quantum` are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Gaussian `.cube` output for electron density and molecular-orbital grids (`cube` module).

## [0.1.0] - not yet released

### Added

- STO-3G basis sets for H–Ne, Obara–Saika s-shell integrals, McMurchie–Davidson nuclear attraction, a cyclic Jacobi eigensolver, and restricted Hartree–Fock SCF.
- McMurchie–Davidson electron-repulsion integrals based on the Hermite expansion.
- 6-31G basis sets.
- Unrestricted Hartree–Fock (UHF).
- Analytic nuclear gradients and geometry optimisation.
- Charged-species RHF.
- Molden orbital output (`write_molden`).
- Optional `serde` feature for result types, and optional `rayon` feature for the Fock build.

### Changed

- Refined analytic gradients.
- Electron-repulsion integrals rewritten on the McMurchie–Davidson Hermite expansion.

### Fixed

- Hartree–Fock and McMurchie–Davidson correctness fixes from the platform-review pass.
