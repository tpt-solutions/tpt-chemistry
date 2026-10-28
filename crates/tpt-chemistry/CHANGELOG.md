# Changelog

All notable changes to `tpt-chemistry` are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this crate adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `pipeline` module: CIF → structure → XRD helpers (requires `core`, `io` and `crystal`).

## [0.1.0] - not yet released

### Added

- Umbrella crate with feature-gated re-exports of every `tpt-chem-*` crate, plus short module aliases such as `tpt_chemistry::md`.
- `tpt-chem` command-line tool with `hf`, `xrd`, `md` and `kinetics` subcommands.
- `manifest` module: run manifests (provenance JSON) written alongside outputs.
- Serde umbrella feature.
- Examples: `h2_hf`, `water_hf`, `water_cube`, `xyz_energy`, `lj_argon_nve`, `nacl_pme`, `water_pme`, `lotka_volterra_ssa`, `michaelis_menten`, `cpp_calcium_phosphate`, `silicon_xrd`, `nacl_xrd` and `cif_to_xrd`.
- Optional `rayon` feature forwarding to parallel MD and Fock builds.
