# tpt-chem-io

Chemical file-format support for [`tpt-chemistry`](../../README.md): readers and
writers for XYZ, PDB, MOL2, SDF (V2000), CIF and GROMACS `.gro`, plus streaming
XYZ trajectory I/O for large MD outputs.

All parsers are hand-written. There is no regex engine and no external parsing
dependency; the only workspace dependency is [`tpt-chem-core`](../tpt-chem-core).

## Supported formats

| Format | Module | Read | Write | Notes |
|---|---|:-:|:-:|---|
| XYZ | `xyz` | ✔ | ✔ | Single-frame and multi-frame (`parse_xyz_all`). |
| XYZ trajectory | `trajectory` | ✔ | ✔ | Streaming `XyzTrajectoryReader` / `XyzTrajectoryWriter` for large files. |
| PDB | `pdb` | ✔ | ✔ | Atom records for biomolecular structures. |
| MOL2 | `mol2` | ✔ | ✔ | Tripos MOL2 atoms and bonds. |
| SDF (V2000) | `sdf` | ✔ | ✔ | Single (`parse_sdf`) and multi-record (`parse_sdf_all`) files. |
| CIF | `cif` | ✔ | ✔ | Crystallographic Information Files; feeds `tpt-chem-crystal`. |
| GROMACS `.gro` | `gro` | ✔ | ✔ | GROMACS coordinate files (positions in nm). |

Errors are reported through `IoError`, and every fallible function returns
`tpt_chem_io::Result<T>`.

## Installation

```toml
[dependencies]
tpt-chem-io = "0.1"
```

This crate has no feature flags.

## Quick start

```rust
use tpt_chem_io::{parse_xyz, write_xyz};

let source = "3\nwater (comment line)\nO  0.000  0.000  0.117\nH  0.000  0.757 -0.469\nH  0.000 -0.757 -0.469\n";
let mol = parse_xyz(source).unwrap();
assert_eq!(mol.formula(), "H2O");

let text = write_xyz(&mol, "round trip");
let again = parse_xyz(&text).unwrap();
assert_eq!(again.atom(0usize.into()).pos, mol.atom(0usize.into()).pos);
```

### Streaming a trajectory

`XyzTrajectoryReader` yields one frame at a time, so a multi-gigabyte
trajectory never has to fit in memory. Use `XyzTrajectoryWriter` to emit
frames as a simulation runs:

```rust
use tpt_chem_io::{XyzTrajectoryReader, XyzTrajectoryWriter};
```

See the type documentation on [docs.rs](https://docs.rs/tpt-chem-io) for the
full frame-iteration API.

## Test data

The integration tests in `tests/data_samples.rs` round-trip the idealized
water, benzene and NaCl samples in the workspace [`data/`](../../data) directory,
so the parsers are exercised on realistic input, not just hand-written
snippets.

## Status

Pre-release (`0.1.x`). SDF support covers V2000 only; V3000 is not yet
supported. See [CHANGELOG.md](CHANGELOG.md).

## Related crates

- [`tpt-chem-core`](../tpt-chem-core) — the `Molecule` type these formats produce and consume.
- [`tpt-chem-crystal`](../tpt-chem-crystal) — consumes CIF structures for XRD simulation.
- [`tpt-chemistry`](../tpt-chemistry) — the `cif_to_xrd` example chains CIF parsing into diffraction.

## License

Dual-licensed under [MIT](../../LICENSE-MIT) or [Apache-2.0](../../LICENSE-APACHE), at your option.
