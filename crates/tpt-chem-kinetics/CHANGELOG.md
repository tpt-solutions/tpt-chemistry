# Changelog

All notable changes to `tpt-chem-kinetics` are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - not yet released

### Added

- Mass-action reaction networks with a builder API.
- RK4 integrator and an implicit-Euler stiff solver.
- Gillespie stochastic simulation algorithm (SSA) and tau-leaping.
- Arrhenius and Eyring temperature-dependent rate laws.
- Steady-state iteration.
- Rate-constant fitting (`fit_rates`) and parameter sensitivity analysis (`sensitivity`).
- Optional `serde` feature for result and network types.
