# tpt-chem-md

The classical molecular dynamics engine of [`tpt-chemistry`](../../README.md).

It provides symplectic Velocity Verlet integration, from-scratch cell lists and
Verlet neighbor lists for O(N) force evaluation, Lennard-Jones and Coulomb
interactions, Ewald summation, and Particle Mesh Ewald (PME) long-range
electrostatics on an in-house radix-2 FFT. There is no LAMMPS, GROMACS or
OpenMM underneath: every algorithm is implemented in safe Rust.

## Units

MD unit system throughout:

| Quantity | Unit |
|---|---|
| Position | Å |
| Time | fs |
| Mass | amu |
| Energy | kJ·mol⁻¹ |
| Force | kJ·mol⁻¹·Å⁻¹ |

Conversion constants are exported at the crate root: `KJ_MOL_PER_MDU` and
`ACC_PER_FORCE_OVER_MASS`.

## Capabilities

- **Integrators** (`integrator`) — `VelocityVerlet`; `Leapfrog` (re-exported at the root).
- **Thermostats and barostats** (`thermostat`) — Berendsen thermostat and barostat, Nosé–Hoover and Langevin steps.
- **Ensembles and run loops** (`ensemble`) — NVE, NVT and NPT run loops.
- **Neighbor search** (`cell`, `neighbors`) — cell lists and Verlet lists with a skin distance.
- **Non-bonded forces** (`forces`) — Lennard-Jones and real-space Coulomb, with `ForceModel` as the extension point.
- **Long-range electrostatics** (`ewald`, `pme`, `bspline`, `fft`) — direct Ewald summation (real, reciprocal and self terms) and Particle Mesh Ewald.
- **Bonded interactions and topology** (`bonded`, `topology`, `exclusions`, `water`) — bonds, angles, exclusions and scaled 1-4 pairs, plus a rigid-water model.
- **Constraints** (`constraints`) — holonomic bond constraints.
- **Periodic box** (`box3`) — `Box3`, including cubic boxes.
- **Minimisation** (`minimiser`) — energy minimisation before dynamics.
- **Observables** (`observables`) — kinetic and potential energy, temperature, virial and pressure.
- **Building systems** (`system`) — `System` for direct construction and `SystemBuilder` for validated, index-checked construction with seeded velocities.

## Installation

```toml
[dependencies]
tpt-chem-md = "0.1"

# Optional features:
# tpt-chem-md = { version = "0.1", features = ["rayon", "serde"] }
```

## Feature flags

| Feature | Effect |
|---|---|
| `rayon` | Parallel force loops over pairs and neighbor lists (uses `rayon`). |
| `serde` | Serialize MD state and topology types (forwards to `tpt-chem-core/serde`). |

## Quick start

NVE dynamics of a small Lennard-Jones cluster:

```rust
use tpt_chem_core::forcefield::LennardJones;
use tpt_chem_core::rng::Rng;
use tpt_chem_md::{integrator::VelocityVerlet, system::System};

let mut sys = System::new();
let lj = LennardJones { sigma: 3.4, epsilon: 0.997 };
for i in 0..13 {
    sys.add_atom(lj, 0.0, 39.95, [i as f64 * 4.0, 0.0, 0.0].into());
}
let mut rng = Rng::new(42);
sys.init_velocities(94.4, &mut rng); // ~94 K: Ar thermal energy
let vv = VelocityVerlet::new(1.0);   // 1 fs steps

let _ = vv.run(&mut sys, 100, None);
let e1 = vv.run(&mut sys, 100, None);
// Symplectic integration: total energy drift stays bounded.
let e_pot = sys.potential_energy;
let e_tot = e_pot + sys.kinetic_energy();
assert!((e_tot - e1).abs() / e1.abs() < 1e-2);
```

The `prelude` module imports the usual set of types:

```rust
use tpt_chem_md::prelude::*;
```

## Examples

Runnable demos live in the umbrella crate; these two exercise this crate directly:

```sh
cargo run -p tpt-chemistry --features core,md --example lj_argon_nve
cargo run -p tpt-chemistry --features core,md --example nacl_pme
```

`lj_argon_nve` reports energy drift over a trajectory. `nacl_pme` compares
direct Ewald against PME convergence. A third example, `water_pme`, runs 64
rigid TIP3P waters with SHAKE/RATTLE constraints and PME electrostatics
(equilibrated to 300 K, then NVE). Run it with `--release`.

## Testing

The crate ships integration tests for `SystemBuilder` (`tests/builder.rs`) and
serde round-tripping (`tests/serde_roundtrip.rs`, run with `--features serde`).
Physical-law property tests (energy and symplecticity conservation) live in
[`tpt-chem-verify`](../tpt-chem-verify).

## Status

Pre-release (`0.1.x`). Electrostatics are checked against direct Ewald sums;
the Ewald reciprocal term is V-normalized as documented in the module. See
[CHANGELOG.md](CHANGELOG.md).

## Related crates

- [`tpt-chem-core`](../tpt-chem-core) — atoms, force-field parameter types and `Vec3`.
- [`tpt-chem-io`](../tpt-chem-io) — trajectory writers for MD output.
- [`tpt-chem-verify`](../tpt-chem-verify) — conservation-law checks for this engine.
- [`tpt-chemistry`](../tpt-chemistry) — umbrella crate with feature-gated re-exports.

## License

Dual-licensed under [MIT](../../LICENSE-MIT) or [Apache-2.0](../../LICENSE-APACHE), at your option.
