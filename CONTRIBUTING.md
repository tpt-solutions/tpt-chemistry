# Contributing to tpt-chemistry

Thanks for your interest in contributing! This document covers the rules that
keep the workspace healthy.

## Licensing

- All code is dual-licensed `MIT OR Apache-2.0` (see `LICENSE-MIT` and
  `LICENSE-APACHE`). By opening a PR you agree to license your contribution
  under both.
- **No Apache-2.0-only dependencies.** Per `spec.txt` §2, every crate that
  enters the dependency graph must be MIT, dual `MIT OR Apache-2.0`, or
  otherwise MIT-compatible. `cargo deny check` enforces this — bare
  `Apache-2.0` is deliberately absent from the allow-list in `deny.toml`.

## Getting started

```sh
git clone <repo>
cd tpt-chemistry
just check     # fmt + clippy + test + cargo deny
```

A local `just` runner (`justfile`) wraps the common commands: `just fmt`,
`just clippy`, `just test`, `just deny`, `just no-std`, `just kani`.

## Ground rules

1. **Zero external solvers.** No LAMMPS/Gaussian/ORCA FFI, no BLAS/LAPACK
   system libraries. Algorithms are implemented from scratch in pure Rust or
   pulled from the in-house `tpt-math` crates.
2. **No `unsafe`.** The workspace lints `unsafe_code = "forbid"`.
3. **Unit-safe at the boundaries.** Physical quantities cross public APIs with
   explicit units (via `tpt-math-units`) or in documented simulation units
   (atomic units, Å/fs/kJ·mol⁻¹).
4. **`no_std + alloc` for the core.** `tpt-chem-core` must build for
   `thumbv6m-none-eabi` with `--no-default-features`. Only use `std` behind
   the default `std` feature.
5. **Physics is verified, not assumed.** New integrators need an energy- or
   phase-space-conservation property test (`tpt-chem-verify`); new numerical
   kernels need Kani harnesses or bounded-input tests where they can panic.
6. **Documentation.** Crate-level docs plus docs on every public item
   (`missing_docs` is a workspace-level warn and denied in CI via `-D
   warnings`).

## Style

- `cargo fmt` (100-column, edition 2021 defaults — see `rustfmt.toml`).
- `cargo clippy --all-targets --all-features -- -D warnings` must be clean.
- Numerics default to `f64`.

## Kani note

The [Kani Rust Verifier](https://github.com/model-checking/kani) does not
support Windows hosts. The proof harnesses in `crates/tpt-chem-verify` are
compiled with `cargo kani` on Linux (see `.github/workflows/ci.yml`), or
locally under WSL.
