# Changelog

All notable changes to `tpt-chem-md` are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `SystemBuilder` with index validation and seeded velocity initialisation.

## [0.1.0] - not yet released

### Added

- Serde feature for MD state and topology types.
- Velocity Verlet integrator, cell lists, Verlet neighbor lists, and Lennard-Jones and Coulomb forces.
- Ewald summation (real-space, reciprocal and self terms), with a V-normalized reciprocal term.
- Particle Mesh Ewald (PME) on an in-house radix-2 FFT.
- Leapfrog integrator; Berendsen thermostat and barostat; Langevin and Nosé–Hoover steps.
- Observables (energy, temperature, virial and pressure), bonded forces, and MD exclusions and constraints.
- Ensemble run loops (NVE, NVT, NPT), an energy minimiser, and a rigid-water model.
- Optional `rayon` feature for parallel force loops.

### Fixed

- Sign of the virial in the pressure calculation.
- Correctness fixes from the platform-review pass.
