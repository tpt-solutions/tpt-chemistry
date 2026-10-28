# tpt-chem-core

The foundation layer of [`tpt-chemistry`](../../README.md): molecular graphs,
chemical elements, force-field parameter types, physical constants and
unit-safe coordinate types. Every other `tpt-chem-*` crate depends on this one.

The crate is `no_std + alloc` compatible. With the default `std` feature
disabled it builds for bare-metal targets, with transcendental floating-point
math routed through [`libm`](https://crates.io/crates/libm).

## Contents

| Module | What it provides |
|---|---|
| `molecule` | `Molecule` — an atom/bond graph with `AtomId` handles, `BondOrder`, distances, molar mass and formula strings (`"H2O"`). |
| `element` | `Element` — the periodic-table data used throughout the workspace. |
| `coords` | `Cartesian`, `Internal` (Z-matrix) and `Position` — phantom-typed coordinate systems, so Cartesian and internal coordinates cannot be mixed by accident. |
| `forcefield` | Parameter types `LennardJones`, `Coulomb`, `HarmonicBond`, `HarmonicAngle` and `RyckaertBellemans`. |
| `potential` | Potential-energy abstractions. |
| `units` | Physical constants (e.g. `BOHR_ANGSTROM`) and the unit-conversion conventions below. |
| `vec3` | `Vec3`, a small 3-vector type used for every position and force. |
| `num` | Numeric helpers shared by the crates. |
| `special` | Special functions needed by the physics code. |
| `rng` | `Rng` — a small, seedable, deterministic generator, so simulations are reproducible. |
| `prelude` | Re-exports of the most-used types above. |

## Units convention

Every quantity crossing a public API is either explicitly unit-wrapped or
documented in one of two unit systems:

- **MD units** — length Å, time fs, mass amu, energy kJ·mol⁻¹.
- **Atomic units** — length Bohr, time ħ/Eₕ, mass mₑ, energy Eₕ (Hartree).

`tpt-chem-quantum` works in atomic units and converts from Å at its API
boundary; `tpt-chem-md` works in MD units.

## Installation

```sh
cargo add tpt-chem-core
```

Or in `Cargo.toml`:

```toml
[dependencies]
tpt-chem-core = "0.1"

# Serialize/deserialize the core types (no_std-compatible via `alloc`):
# tpt-chem-core = { version = "0.1", features = ["serde"] }

# no_std + alloc (bare metal):
# tpt-chem-core = { version = "0.1", default-features = false }
```

## Feature flags

| Feature | Default | Effect |
|---|---|---|
| `std` | yes | Enables the standard library (and `tpt-math-units/std`). Disable for `no_std + alloc`. |
| `serde` | no | Derives `Serialize`/`Deserialize` on the core types. |

## Quick start

Build a water molecule and measure its O–H distance:

```rust
use tpt_chem_core::molecule::{BondOrder, Molecule};

let mut water = Molecule::new("water");
let o = water.add_atom::<8>([0.000, 0.000, 0.117].into());
let h1 = water.add_atom::<1>([0.000, 0.757, -0.469].into());
let h2 = water.add_atom::<1>([0.000, -0.757, -0.469].into());
water.add_bond(o, h1, BondOrder::Single).unwrap();
water.add_bond(o, h2, BondOrder::Single).unwrap();

assert!((water.distance(o, h1) - 0.9577).abs() < 1e-3);
assert_eq!(water.formula(), "H2O");
```

Atoms are added with their atomic number as the const generic parameter
(`8` for oxygen, `1` for hydrogen), and positions are given in Ångström.

## Verifying the no_std build

The workspace `justfile` checks that the core builds for a bare-metal target:

```sh
just no-std    # cargo build for thumbv6m-none-eabi
```

## Status

Pre-release (`0.1.x`): expect additions and occasional breaking changes
until `1.0`. See
[CHANGELOG.md](CHANGELOG.md) for the history.

## Related crates

- [`tpt-chem-io`](../tpt-chem-io) — read and write these molecules as XYZ, PDB, MOL2, SDF and CIF.
- [`tpt-chem-md`](../tpt-chem-md) — molecular dynamics using these force-field and vector types.
- [`tpt-chem-quantum`](../tpt-chem-quantum) — Hartree–Fock on top of `Molecule`.
- [`tpt-chemistry`](../tpt-chemistry) — umbrella crate with feature-gated re-exports.

## License

Dual-licensed under [MIT](../../LICENSE-MIT) or [Apache-2.0](../../LICENSE-APACHE), at your option.
