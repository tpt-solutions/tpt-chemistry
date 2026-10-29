# Quick start

tpt-chemistry depends on the sibling repository `tpt-math` through path
dependencies, so clone both side by side:

```sh
git clone https://github.com/tpt-solutions/tpt-math
git clone https://github.com/tpt-solutions/tpt-chemistry
cd tpt-chemistry
cargo test --workspace
```

## Run a Hartree–Fock energy

The starter example reads an XYZ file (coordinates in Å) and prints the
RHF/STO-3G energy:

```sh
cargo run -p tpt-chemistry --features core,io,quantum --example xyz_energy -- data/water.xyz
```

Output for `data/water.xyz`:

```text
molecule : water (H2O)
atoms    : 3
basis fns: 7
energy   : <value> Hartree
```

## Make your own project

The `templates/starter` directory is a [cargo-generate](https://github.com/cargo-generate/cargo-generate)
template with the same program:

```sh
cargo generate --path templates/starter --name my-calc
```

The generated project depends on the umbrella crate by local path. Its
placeholder defaults to `../tpt-chemistry/crates/tpt-chemistry`.

## Other examples

Each example is feature-gated, so pass the features it needs with
`--features` and the name with `--example`:

| Example | Shows |
| --- | --- |
| `h2_hf` | Bond-length scan of H₂ |
| `water_hf` | Orbitals, energies, Molden output |
| `lj_argon_nve` | NVT equilibration then NVE drift |
| `nacl_pme` | Ewald versus PME electrostatics |
| `water_pme` | Rigid SPC/TIP3P water with PME |
| `lotka_volterra_ssa` | Stochastic simulation |
| `michaelis_menten` | Stochastic versus deterministic kinetics |
| `silicon_xrd`, `nacl_xrd` | Powder XRD patterns |
| `cif_to_xrd` | End-to-end CIF to XRD |
| `cpp_calcium_phosphate` | Aqueous speciation scan |
