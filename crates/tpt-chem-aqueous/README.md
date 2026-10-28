# tpt-chem-aqueous

Aqueous-solution thermodynamics for [`tpt-chemistry`](../../README.md):
activity coefficients, solubility products and saturation indices, and
chelator-controlled calcium-phosphate precipitation.

## Capabilities

| Module | What it provides |
|---|---|
| `activity` | `ActivityModel` implementations for multi-component electrolytes: limiting-law and extended Debye–Hückel, Davies, Pitzer virial equations, and single-solvent electrolyte NRTL (eNRTL). |
| `equilibrium` | Van't Hoff temperature dependence of equilibrium constants. |
| `solubility` | Solid phases, ion-activity products and saturation indices. |
| `phosphate` | The H₃PO₄ / calcium-phosphate constant set and solid phases. |
| `chelator` | Multi-site Ca²⁺-binding ligands, such as phosphopeptides. |
| `precipitation` | pH × temperature scans and onset-pH search, predicting when a calcium-phosphate phase becomes supersaturated ("chalkiness"). Includes a kinetic nucleation onset. |
| `species`, `solvent`, `constants` | Ionic species descriptions, water properties and tabulated constants. |
| `solver` | Internal root-finding used by the equilibrium solvers. |

Errors are reported through `AqueousError`, and fallible functions return
`tpt_chem_aqueous::Result<T>`.

## Units

- Temperatures in kelvin.
- Concentrations as molalities (mol per kg of water).
- Enthalpies in J·mol⁻¹.
- pH on the activity scale: `a(H⁺) = 10^(−pH)`.
- Equilibrium constants as `log₁₀ K` at 25 °C.

## Installation

```toml
[dependencies]
tpt-chem-aqueous = "0.1"
```

This crate has no dependencies and no feature flags.

## Quick start

Mean activity coefficient of 1 m NaCl from the Pitzer model:

```rust
use tpt_chem_aqueous::activity::{ActivityModel, pitzer::{Pitzer, PitzerPair}};
use tpt_chem_aqueous::species::Species;

let mut pitzer = Pitzer::new();
pitzer.add_pair(PitzerPair::new("Na+", "Cl-", 0.0765, 0.2664, 0.0, 0.00127));
let sol = [Species::new("Na+", 1, 1.0), Species::new("Cl-", -1, 1.0)];
let ln_g = pitzer.ln_gamma(&sol, 298.15).unwrap();
let gamma_pm = ((ln_g[0] + ln_g[1]) / 2.0).exp();
assert!((gamma_pm - 0.657).abs() < 0.005);
```

## Examples

```sh
cargo run -p tpt-chemistry --features core,aqueous --example cpp_calcium_phosphate
```

This scans calcium-phosphate onset pH against temperature, with and without a
chelator.

## Important caveats

- **Thermodynamic, not kinetic.** A saturation index above zero means the
  solution is supersaturated with respect to a phase. It does *not* mean the
  phase has nucleated. The kinetic nucleation onset is a separate, explicitly
  labelled model.
- **Default constants are illustrative.** Tabulated constants shipped as
  defaults are literature-typical values. Replace them with measured data for
  quantitative work; every value is overridable.

## Testing

Proptest-based verification of the activity models and solvers lives in
[`tpt-chem-verify`](../tpt-chem-verify) (`aqueous` module).

## Status

Pre-release (`0.1.x`). Activity models are exercised by doc examples and by
the property tests in `tpt-chem-verify`. See [CHANGELOG.md](CHANGELOG.md).

## Related crates

- [`tpt-chem-kinetics`](../tpt-chem-kinetics) — kinetic models that consume equilibrium constants.
- [`tpt-chemistry`](../tpt-chemistry) — umbrella crate; enable with the `aqueous` feature.

## License

Dual-licensed under [MIT](../../LICENSE-MIT) or [Apache-2.0](../../LICENSE-APACHE), at your option.
