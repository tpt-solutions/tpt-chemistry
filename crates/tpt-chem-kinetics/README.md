# tpt-chem-kinetics

Reaction kinetics and population dynamics for [`tpt-chemistry`](../../README.md).

The crate covers three regimes:

- **Deterministic mass-action kinetics** — ODE systems integrated with classic Runge–Kutta (RK4), or with an implicit-Euler stiff solver for systems with widely separated time scales.
- **Stochastic kinetics** — the Gillespie stochastic simulation algorithm (SSA) and tau-leaping, for low-copy-number regimes where fluctuations matter.
- **Temperature dependence** — Arrhenius and Eyring rate laws.

It also provides equilibrium and steady-state analysis, plus rate-constant fitting
and sensitivity analysis against observed data.

## Capabilities

| Module | What it provides |
|---|---|
| `network` | `ReactionNetwork` built with a fluent builder; `Reaction` and `ReactionInput` describe stoichiometry and rate constants. |
| `solver` | `integrate` (RK4), `integrate_stiff` (implicit Euler) and `ssa` (Gillespie). |
| `arrhenius` | Arrhenius and Eyring temperature-dependent rate constants. |
| `fit` | `fit_rates` to estimate rate constants from `Observation`s, with `FitSettings` and `FitResult`; `sensitivity` for parameter sensitivities. |

## Installation

```toml
[dependencies]
tpt-chem-kinetics = "0.1"

# Serialize result and network types:
# tpt-chem-kinetics = { version = "0.1", features = ["serde"] }
```

## Feature flags

| Feature | Effect |
|---|---|
| `serde` | Serialize result and network types. |

## Quick start

First-order decay A → B, integrated deterministically:

```rust
use tpt_chem_kinetics::{network::ReactionNetwork, solver};

// A → B with rate constant k = 1.0 (time units are whatever the
// rate constants are expressed in).
let net = ReactionNetwork::builder()
    .species(&["A", "B"])
    .reaction(&[("A", 1.0)], &[], 1.0)
    .build();
let y0 = [10.0, 0.0];

let traj = solver::integrate(&net, &y0, 0.0, 5.0, 500);
let a_final = traj.last().unwrap()[0];
assert!((a_final - 10.0 * (-5.0_f64).exp()).abs() < 1e-6);
```

## Examples

```sh
cargo run -p tpt-chemistry --features core,kinetics --example lotka_volterra_ssa  # Gillespie SSA vs mass-action ODE
cargo run -p tpt-chemistry --features core,kinetics --example michaelis_menten    # enzyme kinetics, both regimes
```

## Testing

- Mass conservation for closed reaction networks is property-tested in [`tpt-chem-verify`](../tpt-chem-verify) (RK4, stiff solver, SSA, randomized states).
- `tests/serde_roundtrip.rs` checks serialization (`--features serde`).

## Status

Pre-release (`0.1.x`). Networks are mass-action; enzyme kinetics such as
Michaelis–Menten are expressed as mass-action networks (see the
`michaelis_menten` example). See [CHANGELOG.md](CHANGELOG.md).

## Related crates

- [`tpt-chem-core`](../tpt-chem-core) — shared chemistry types and constants.
- [`tpt-chem-aqueous`](../tpt-chem-aqueous) — thermodynamics whose equilibrium constants feed kinetic models.
- [`tpt-chemistry`](../tpt-chemistry) — umbrella crate with feature-gated re-exports.

## License

Dual-licensed under [MIT](../../LICENSE-MIT) or [Apache-2.0](../../LICENSE-APACHE), at your option.
