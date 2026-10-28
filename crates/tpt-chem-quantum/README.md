# tpt-chem-quantum

Electronic-structure primitives for [`tpt-chemistry`](../../README.md): Gaussian
basis sets, one- and two-electron integrals, and Hartree–Fock self-consistent
field (SCF) solvers, all implemented from scratch in Rust.

- Gaussian-type-orbital basis sets: STO-3G for H–Ne, plus 6-31G.
- Obara–Saika overlap and kinetic integrals, and McMurchie–Davidson nuclear-attraction and electron-repulsion integrals (ERIs), built on the Hermite expansion.
- Boys function evaluation for the Gaussian integrals.
- Restricted closed-shell Hartree–Fock (RHF) and unrestricted Hartree–Fock (UHF).
- Charged-species RHF.
- Analytic nuclear gradients and geometry optimisation.
- The Roothaan–Hall equations are solved with a from-scratch cyclic Jacobi eigensolver, so the crate has no LAPACK or `nalgebra` dependency.
- Output: Molden orbital files (`write_molden`) and Gaussian `.cube` grids of electron density and molecular orbitals (`cube`).

Everything runs in **atomic units** (Bohr, Hartree). Molecular geometries from
[`tpt-chem-core`](../tpt-chem-core) are in Ångström and are converted at the API
boundary.

## Installation

```toml
[dependencies]
tpt-chem-quantum = "0.1"

# Optional features:
# tpt-chem-quantum = { version = "0.1", features = ["rayon", "serde"] }
```

## Feature flags

| Feature | Effect |
|---|---|
| `rayon` | Parallel Fock-matrix build over shell quartets (uses `rayon`). |
| `serde` | Serialize result types (forwards to `tpt-chem-core/serde`). |

## Quick start

Restricted Hartree–Fock for H₂ at its experimental bond length:

```rust
use tpt_chem_core::molecule::Molecule;
use tpt_chem_core::units::BOHR_ANGSTROM;
use tpt_chem_core::vec3::Vec3;
use tpt_chem_quantum::hf::rhf_energy;

let mut h2 = Molecule::new("H2");
let r = 1.4 * BOHR_ANGSTROM; // 1.4 Bohr in Angstrom
h2.add_atom::<1>(Vec3::new(0.0, 0.0, -r / 2.0));
h2.add_atom::<1>(Vec3::new(0.0, 0.0, r / 2.0));

// HF/STO-3G literature value: −1.11675931 Eₕ (Schwartz & Schaad 1967).
let e = rhf_energy(&h2).unwrap();
assert!((e - (-1.11675931)).abs() < 5e-4);
```

The `prelude` module collects the main entry points:

```rust
use tpt_chem_quantum::prelude::*; // BasisSet, rhf, rhf_with_basis, HfResult,
                                  // uhf_with_basis, optimize_geometry, write_molden, …
```

## Examples

```sh
cargo run -p tpt-chemistry --features core,quantum --example h2_hf      # bond-length scan
cargo run -p tpt-chemistry --features core,quantum --example water_hf   # SCF energies, orbitals, Molden
cargo run -p tpt-chemistry --features core,quantum --example water_cube # density and MO grids
```

## Accuracy

- H₂ / STO-3G at R = 1.4 Bohr agrees with the literature value to within 75 µEₕ, and is tested against it.
- Water with STO-3G reproduces the reference energies exactly after the McMurchie–Davidson ERI rewrite.

## Testing

`tests/serde_roundtrip.rs` checks result serialization (`--features serde`).
Physical-law property tests are in [`tpt-chem-verify`](../tpt-chem-verify).

## Status

Pre-release (`0.1.x`). Basis-set coverage is limited to STO-3G (H–Ne) and
6-31G; larger correlation-consistent bases and post-HF methods (MP2, CCSD) are
not implemented. See [CHANGELOG.md](CHANGELOG.md).

## Related crates

- [`tpt-chem-core`](../tpt-chem-core) — the `Molecule` type and unit conventions.
- [`tpt-chem-io`](../tpt-chem-io) — read molecular geometries from XYZ, MOL2 and other formats.
- [`tpt-chemistry`](../tpt-chemistry) — umbrella crate with feature-gated re-exports.

## License

Dual-licensed under [MIT](../../LICENSE-MIT) or [Apache-2.0](../../LICENSE-APACHE), at your option.
