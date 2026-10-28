# Changelog

All notable changes to `tpt-chem-io` are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- SDF (V2000) reader and writer: `parse_sdf`, `parse_sdf_all` and `write_sdf`.
- Idealized water, benzene and NaCl sample structures, used by the `data_samples` integration tests.

## [0.1.0] - not yet released

### Added

- XYZ, PDB, MOL2 and CIF parsers and writers, all hand-written with no regex dependency.
- Streaming XYZ trajectory I/O: `XyzTrajectoryReader` and `XyzTrajectoryWriter`.
- GROMACS `.gro` reader and writer.

### Changed

- Workspace-wide rustdoc link fixes; the workspace documents with zero warnings.
