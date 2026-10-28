# tpt-chem-crystal

Solid-state physics and crystallography for [`tpt-chemistry`](../../README.md):
Bravais lattices, Miller indices, Wigner–Seitz cells, space-group symmetry with
Wyckoff orbit generation, reciprocal-lattice k-point grids, and powder X-ray
diffraction (XRD) simulation.

Atomic positions are fractional in the direct lattice. Lattice vectors are the
rows of the conventional cell.

## Capabilities

| Module | What it provides |
|---|---|
| `bravais` | `Bravais` centering types and `Lattice` construction from cell parameters. |
| `miller` | `Miller` indices and reciprocal-lattice vectors. |
| `wigner` | Wigner–Seitz cell construction (half-space intersection with a rank-3 vertex filter). |
| `symmetry` | Space-group operations for P1 through Fm-3m, with Wyckoff orbit generation. |
| `kpoints` | `KGrid` and `KPoint` — Monkhorst–Pack style reciprocal-space grids. |
| `build` | Assembling crystal structures from lattice, basis and symmetry. |
| `xrd` | `XrdSimulator` and `Reflection` — powder diffraction patterns from structure factors. |

The XRD simulator uses Cromer–Mann atomic form factors, the Lorentz–polarization
factor, and the systematic-extinction rules of the space group. Silicon (diamond)
and sodium chloride (rock salt) are validated against reference patterns.

## Installation

```toml
[dependencies]
tpt-chem-crystal = "0.1"

# Serialize result and network types:
# tpt-chem-crystal = { version = "0.1", features = ["serde"] }
```

## Feature flags

| Feature | Effect |
|---|---|
| `serde` | Serialize result types (forwards to `tpt-chem-core/serde`). |

## Quick start

The examples show the typical flow. Silicon's diamond-structure powder pattern:

```sh
cargo run -p tpt-chemistry --features core,crystal --example silicon_xrd   # prints CSV
cargo run -p tpt-chemistry --features core,crystal --example nacl_xrd      # rock-salt extinctions
```

Parsing a CIF file and simulating its diffraction pattern runs end-to-end
through [`tpt-chem-io`](../tpt-chem-io):

```sh
cargo run -p tpt-chemistry --features core,io,crystal --example cif_to_xrd
```

## Testing

`tests/serde_roundtrip.rs` checks serialization (`--features serde`). Symmetry
and XRD results are checked against the reference diamond and rock-salt
patterns.

## Status

Pre-release (`0.1.x`). Space-group symmetry is implemented from P1 through
Fm-3m; the full set of 230 space groups is not yet covered. See
[CHANGELOG.md](CHANGELOG.md).

## Related crates

- [`tpt-chem-io`](../tpt-chem-io) — reads CIF structures into this crate's inputs.
- [`tpt-chem-core`](../tpt-chem-core) — shared types and units.
- [`tpt-chemistry`](../tpt-chemistry) — umbrella crate; enable with the `crystal` feature.

## License

Dual-licensed under [MIT](../../LICENSE-MIT) or [Apache-2.0](../../LICENSE-APACHE), at your option.
