# tpt-chem-verify

Verification harnesses for [`tpt-chemistry`](../../README.md). This crate holds
the checks that prove the physics stays correct as the rest of the workspace
changes: property-based tests, reusable conservation-law assertions, and
edge-case properties at boundaries and coincident centres.

It is a test-support crate. It is not needed to use the simulators; the
umbrella crate exposes it only behind the `verify` feature.

## What is here

| Module | What it provides |
|---|---|
| `strategies` | `proptest` strategies that generate valid molecular geometries and force-field parameters, plus `edge_f64` corner values. |
| `conservation` | Reusable conservation-law assertions, such as `assert_energy_conserved` for MD and symplecticity checks. |
| `mass_conservation` | Mass-conservation properties for closed reaction networks: RK4, the stiff solver, the SSA, and randomized states. |
| `edge_cases` | Property tests, compiled only under `cfg(test)`: cell-index bounds, the Verlet rebuild criterion (open and periodic), SCF density positive semi-definiteness, and integrals at coincident centres and near-zero separations. |

The aqueous-model property tests are a private `#[cfg(test)]` module.

## Installation

Add it as a dev-dependency so it only builds for tests:

```toml
[dev-dependencies]
tpt-chem-verify = "0.1"
```

## Quick start

Reuse the energy-conservation check on a custom system:

```rust
use tpt_chem_core::forcefield::LennardJones;
use tpt_chem_core::rng::Rng;
use tpt_chem_core::vec3::Vec3;
use tpt_chem_md::box3::Box3;
use tpt_chem_md::system::System;
use tpt_chem_verify::conservation::assert_energy_conserved;

let mut sys = System::new();
let lj = LennardJones { sigma: 3.4, epsilon: 0.997 };
let mut rng = Rng::new(1);
for i in 0..8 {
    sys.add_atom(lj, 0.0, 39.95, Vec3::new(i as f64 * 4.5 + 2.0, 2.0, 2.0));
}
sys.set_box(Box3::cubic(40.0));
sys.init_velocities(94.0, &mut rng);
// 1 fs Velocity Verlet on a clean geometry conserves energy tightly.
assert_energy_conserved(&mut sys, 200, 1.0, 5e-3);
```

## Running the checks

```sh
cargo test -p tpt-chem-verify            # proptest properties and conservation checks
```

## Adding a new property

1. Put a generator in `strategies` if you need a new kind of random input.
2. Put a reusable invariant in `conservation` or `mass_conservation`.
3. Call it from a `proptest!` block in your test file.

## Status

Pre-release (`0.1.x`). `proptest` regression seeds are stored in
`proptest-regressions/`; commit them so failures reproduce. See
[CHANGELOG.md](CHANGELOG.md).

## Related crates

- [`tpt-chem-md`](../tpt-chem-md), [`tpt-chem-kinetics`](../tpt-chem-kinetics), [`tpt-chem-quantum`](../tpt-chem-quantum), [`tpt-chem-aqueous`](../tpt-chem-aqueous) — the crates whose invariants are checked here.

## License

Dual-licensed under [MIT](../../LICENSE-MIT) or [Apache-2.0](../../LICENSE-APACHE), at your option.
