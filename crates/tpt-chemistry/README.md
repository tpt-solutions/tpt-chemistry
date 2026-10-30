# tpt-chemistry

A pure-Rust, AI-native computational chemistry and molecular simulation
library: molecular dynamics, quantum chemistry (Hartree–Fock with Gaussian
basis sets), reaction kinetics, aqueous thermodynamics, and crystallography.

This is the **umbrella crate**. It re-exports the `tpt-chem-*` crates behind
Cargo features, so you can pull in only the physics you need. It also ships the
`tpt-chem` command-line tool and the runnable examples for the whole workspace.

No C/C++ FFI, no external solver. No LAMMPS, no Gaussian, no Quantum ESPRESSO:
SCF iteration, McMurchie–Davidson integrals, Ewald summation, Particle Mesh
Ewald, neighbor lists and the solvers are all implemented from scratch.

## Feature flags

| Feature | Crate | Contents |
|---|---|---|
| `core` | [`tpt-chem-core`](../tpt-chem-core) | Molecular graphs, elements, force fields, unit-safe types (`no_std`). |
| `io` | [`tpt-chem-io`](../tpt-chem-io) | XYZ, PDB, MOL2, SDF, CIF and GRO readers and writers; trajectory streaming. |
| `md` | [`tpt-chem-md`](../tpt-chem-md) | Integrators, thermostats and barostats, cell and neighbor lists, Ewald and PME. |
| `quantum` | [`tpt-chem-quantum`](../tpt-chem-quantum) | Gaussian basis sets, ERIs, RHF/UHF SCF, geometry optimisation, Molden and cube output. |
| `kinetics` | [`tpt-chem-kinetics`](../tpt-chem-kinetics) | Mass-action networks, stiff and stochastic solvers, Arrhenius/Eyring, rate fitting. |
| `aqueous` | [`tpt-chem-aqueous`](../tpt-chem-aqueous) | Activity coefficients (Debye–Hückel, Davies, Pitzer, eNRTL), Ksp, Ca-phosphate precipitation. |
| `crystal` | [`tpt-chem-crystal`](../tpt-chem-crystal) | Bravais lattices, space groups, k-point grids, XRD simulation. |
| `verify` | [`tpt-chem-verify`](../tpt-chem-verify) | proptest strategies, conservation checks, edge-case properties. |
| `rayon` | — | Parallel MD force loops and Fock build (forwarded to `tpt-chem-md` and `tpt-chem-quantum`). |
| `serde` | — | Serialization of result and state types across the enabled crates. |

Every feature is off by default, so an empty build compiles only the manifest
and run-provenance helpers.

## Installation

```toml
[dependencies]
tpt-chemistry = { version = "0.1", features = ["core", "quantum", "md"] }
```

Each module is available under two names, for example `tpt_chemistry::md` and
`tpt_chemistry::tpt_chem_md`. The short names are used in the examples below.

## Quick start

```rust
use tpt_chemistry::prelude::*;

let mut water = Molecule::new("water");
let o = water.add_atom::<8>([0.0, 0.0, 0.117].into());
let h1 = water.add_atom::<1>([0.0, 0.757, -0.469].into());
let h2 = water.add_atom::<1>([0.0, -0.757, -0.469].into());
water.add_bond(o, h1, BondOrder::Single).unwrap();
water.add_bond(o, h2, BondOrder::Single).unwrap();
assert_eq!(water.formula(), "H2O");
```

The `prelude` module requires the `core` feature.

## Command-line tool

The `tpt-chem` binary has four subcommands: `hf`, `xrd`, `md` and `kinetics`.
It needs every physics feature to build:

```sh
cargo run -p tpt-chemistry --features core,io,quantum,crystal,md,kinetics --bin tpt-chem -- --help
```

## Examples

Run any example with `cargo run -p tpt-chemistry --features <features> --example <name>`.

| Example | Features | Shows |
|---|---|---|
| `h2_hf` | `core`, `quantum` | RHF/STO-3G bond-length scan and dissociation. |
| `water_hf` | `core`, `quantum` | SCF energies, orbital spectrum, Molden output. |
| `water_cube` | `core`, `quantum` | RHF electron density and HOMO written as Gaussian `.cube` files. |
| `xyz_energy` | `core`, `io`, `quantum` | Read an XYZ geometry and print its RHF/STO-3G energy. |
| `lj_argon_nve` | `core`, `md` | NVT equilibration, then NVE production and energy drift. |
| `nacl_pme` | `core`, `md` | Direct Ewald vs Particle Mesh Ewald convergence. |
| `water_pme` | `core`, `md` | 64 rigid TIP3P waters with SHAKE/RATTLE constraints and PME electrostatics. |
| `lotka_volterra_ssa` | `core`, `kinetics` | Gillespie SSA vs mass-action ODE. |
| `michaelis_menten` | `core`, `kinetics` | Enzyme kinetics in both regimes. |
| `cpp_calcium_phosphate` | `core`, `aqueous` | Calcium-phosphate onset pH vs temperature, with and without a chelator. |
| `silicon_xrd` | `core`, `crystal` | Diamond-structure powder pattern as CSV. |
| `nacl_xrd` | `core`, `crystal` | Rock-salt extinctions and peak families. |
| `cif_to_xrd` | `core`, `io`, `crystal` | CIF text → structure → diffraction pipeline. |

## Modules

- **`manifest`** — run manifests: provenance JSON written alongside simulation outputs, so every result records the inputs and versions that produced it.
- **`pipeline`** (requires `core`, `io`, `crystal`) — CIF → structure → XRD helpers.

## Building the workspace

This workspace depends on the sibling repository
[`tpt-math`](https://github.com/tpt-solutions/tpt-math) through path
dependencies (`../tpt-math/...`). A fresh clone needs both checkouts side by side:

```sh
git clone https://github.com/tpt-solutions/tpt-math.git
git clone <this repo> && cd tpt-chemistry
cargo build && cargo test
```

## Status

Pre-release (`0.1.x`). The umbrella API mirrors the constituent crates, so
breaking changes there are reflected here. See [CHANGELOG.md](CHANGELOG.md).

## Related crates

The crates behind each feature are listed in the table above. Design rationale
lives in the workspace [`spec.txt`](../../spec.txt), and build progress in
[`todo.md`](../../todo.md).

## License

Dual-licensed under [MIT](../../LICENSE-MIT) or [Apache-2.0](../../LICENSE-APACHE), at your option.
