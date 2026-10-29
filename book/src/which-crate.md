# Which crate do I need?

Start with the umbrella crate `tpt-chemistry` and turn on the features you
use. Depend on a single `tpt-chem-*` crate only when you want the smallest
build.

| I want to… | Crate or feature | Notes |
| --- | --- | --- |
| Parse and write structures | `tpt-chem-io` (`io`) | XYZ, PDB, MOL2, SDF, CIF, GRO, SMILES, Gaussian |
| Build molecules and force-field parameters | `tpt-chem-core` (`core`) | `no_std + alloc` |
| Run classical MD | `tpt-chem-md` (`md`) | Ewald, PME, thermostats, SHAKE, ML potential |
| Run Hartree–Fock or optimise geometry | `tpt-chem-quantum` (`quantum`) | STO-3G and 6-31G (H–Ne) |
| Simulate reaction networks | `tpt-chem-kinetics` (`kinetics`) | ODE, SSA, tau-leaping |
| Work with crystals and XRD | `tpt-chem-crystal` (`crystal`) | Limited space-group set |
| Compute aqueous activities | `tpt-chem-aqueous` (`aqueous`) | Constants are illustrative |
| Check conservation laws | `tpt-chem-verify` (`verify`) | proptest harnesses |
| Use it from Python | `tpt-chem-py` | Built with maturin, see [Python and the browser](bindings.md) |
| Run HF in a browser | `tpt-chem-wasm` | Built with wasm-pack |
