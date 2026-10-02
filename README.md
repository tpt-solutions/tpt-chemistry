# tpt-chemistry

A pure-Rust, AI-native computational chemistry and molecular simulation
library: molecular dynamics (MD), quantum chemistry primitives
(Hartree–Fock/DFT groundwork), reaction kinetics, and crystallography —
**without C/C++ FFI or external solver dependencies**.

No LAMMPS. No Gaussian. No Quantum ESPRESSO. SCF iteration, electron
repulsion integrals, Ewald summation, and neighbor lists are all implemented
from scratch, leaning on the sibling [`tpt-math`](../tpt-math) crates for
linear algebra and FFTs.

Part of the TPT science stack. Full design rationale lives in
[`spec.txt`](spec.txt); the build progress in [`todo.md`](todo.md).

## Design principles

- **Zero external solver dependencies** — everything from scratch in safe,
  pure Rust.
- **Strict licensing** — workspace is `MIT OR Apache-2.0`, and the dependency
  graph is held to a *stricter* bar: no Apache-2.0-only crates
  (`deny.toml` enforces it in CI).
- **Type-state & dimensional safety** — phantom types distinguish Cartesian
  from internal (Z-matrix) coordinates; physical quantities are wrapped with
  compile-time units (`tpt-math-units`).
- **Native formal verification** — Kani bounded model checking for critical
  loops, `proptest` for statistical verification of physical laws (energy
  conservation, mass conservation, phase-space volume preservation).
- **`no_std + alloc` core** — `tpt-chem-core` builds for bare-metal targets.

## Crates

| Crate | Purpose |
|---|---|
| [`tpt-chem-core`](crates/tpt-chem-core) | Molecular graphs, force fields, physical constants, unit-safe types |
| [`tpt-chem-md`](crates/tpt-chem-md) | Integrators (Velocity Verlet, Leapfrog), thermostats, cell/neighbor lists, Ewald/PME |
| [`tpt-chem-quantum`](crates/tpt-chem-quantum) | Gaussian basis sets, Obara–Saika ERIs, SCF / Hartree–Fock |
| [`tpt-chem-kinetics`](crates/tpt-chem-kinetics) | Mass-action ODEs, Gillespie SSA, Arrhenius/Eyring rates |
| [`tpt-chem-crystal`](crates/tpt-chem-crystal) | Bravais lattices, space groups, reciprocal space, XRD simulation |
| [`tpt-chem-io`](crates/tpt-chem-io) | XYZ, PDB, MOL2, CIF readers/writers, trajectory streaming |
| [`tpt-chem-verify`](crates/tpt-chem-verify) | proptest strategies, conservation-law checks, Kani harnesses |
| [`tpt-chemistry`](crates/tpt-chemistry) | Feature-gated umbrella crate re-exporting all of the above |

## Synergies (spec.txt §5)

- **`tpt-materials`** — micro-to-macro: MD/DFT results feed homogenized
  properties (elastic constants, thermal conductivity).
- **`tpt-math`** — dense linear algebra and ODE machinery.
- **`tpt-dsp` / `tpt-math-signal-fft`** — 3D FFTs for Particle Mesh Ewald and
  vibrational density of states.
- **`tpt-fem` / `tpt-physics`** — QM/MM coupling of the quantum region into
  continuum FEM/CFD environments.
- **AI-native agents (`tpt-eve`, `tpt-anima`)** — forces and energy surfaces
  expose autodiff hooks for training machine-learned interatomic potentials
  and differentiable molecular design.

## Quick start

```sh
just check      # fmt + clippy + test + cargo deny
just no-std     # verify the core builds for thumbv6m-none-eabi
just kani       # bounded model checking (Linux/WSL only)
```

```rust
use tpt_chem_core::molecule::{BondOrder, Molecule};

let mut mol = Molecule::new("H2");
let h1 = mol.add_atom::<2>([0.0, 0.0, 0.0].into());
let h2 = mol.add_atom::<2>([0.0, 0.0, 0.74].into());
mol.add_bond(h1, h2, BondOrder::Single).unwrap();
assert_eq!(mol.mass(), 2.0 * 1.008);
```

Hartree–Fock for H₂ at R = 1.4 Bohr (STO-3G, lit. −1.11676 Eₕ):

```rust
use tpt_chem_core::molecule::Molecule;
use tpt_chem_core::units::BOHR_ANGSTROM;
use tpt_chem_core::vec3::Vec3;
use tpt_chem_quantum::hf::rhf_energy;

let mut h2 = Molecule::new("H2");
let r = 1.4 * BOHR_ANGSTROM;
h2.add_atom::<1>(Vec3::new(0.0, 0.0, -r / 2.0));
h2.add_atom::<1>(Vec3::new(0.0, 0.0, r / 2.0));
assert!((rhf_energy(&h2).unwrap() - (-1.11675931)).abs() < 5e-4);
```

## License

Dual-licensed under [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE), at
your option.
