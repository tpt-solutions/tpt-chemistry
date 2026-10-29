# Validation

These are the checks that the test suite runs. Each row names the test that
enforces it, so you can read the assertion yourself. Most rows compare against
values recorded in this repository. The section
[Independent cross-checks](#independent-cross-checks) compares against other
programs and measured data.

## Quantum: Hartree–Fock energies (STO-3G)

| System | Energy (Eₕ) | Tolerance | Test |
| --- | --- | --- | --- |
| H₂ at R = 0.7414 Å | −1.11668439 | 1e-6 | `h2_sto3g_literature_value` (`tpt-chem-quantum/src/hf.rs`) |
| He (BSE exponents) | −2.80778516 | 1e-6 | `he_atom_sto3g` |
| H₂O, Szabo–Ostlund geometry | −74.96590117 | 1e-6 | `water_sto3g_near_equilibrium` |

The H₂ test records a point about literature values. The Szabo–Ostlund value
−1.11675931 Eₕ is quoted at R = 1.4 Bohr. At 0.7414 Å the same basis gives
−1.11668439 Eₕ. The He difference from the quoted literature value is
12 µEₕ, and the source attributes it to rounding in the contraction
coefficients.

## Quantum: gradients and optimisation

Analytic gradients for H₂ and water agree with finite differences
(`tpt-chem-quantum/src/gradient.rs`). Geometry optimisation reproduces the
HF/STO-3G minima for H₂ and water. The water test is `#[ignore]` in debug
builds, so run it with `cargo test --release -- --ignored`.

## Molecular dynamics

| Check | Result | Source |
| --- | --- | --- |
| PME against direct Ewald, 32³ mesh | matching-cutoff energies agree to 5.4e-5 relative | `pme_matches_ewald_with_matching_wavenumber_cutoff` |
| PME forces | exact gradient of the mesh energy (finite differences) | `tpt-chem-md/src/pme.rs` |
| Velocity Verlet | energy conserved across proptest runs | `tpt-chem-verify` |
| Rigid water, PME, NVE | about 1 kJ·mol⁻¹ drift over 500 fs in a 64-molecule box | `examples/water_pme.rs` |
| Minimisation, L-BFGS | 82 steps versus 522 for steepest descent on a 27-atom cluster | `minimiser::tests` |
| Radial ML potential | analytic forces match finite differences; NVE drift below 1e-3 relative | `tpt-chem-md/src/mlp.rs` |

## Parallel speed-up

Benzene RHF/STO-3G (36 basis functions), release build, one run each, on the
development machine recorded in `todo.md`: serial 122.0 s, rayon 39.0 s (about
3.1×), same energy −227.89100277 Eₕ. Your timings will differ.

## Benchmarks

`cargo bench` runs criterion benchmarks in `tpt-chem-md` and `tpt-chem-quantum`.
A Windows development run recorded these times:

| Benchmark | Time |
| --- | --- |
| All-pairs LJ + Coulomb, 216 argon atoms, cutoff 8.5 Å | 0.37–0.57 ms per evaluation (two runs) |
| RHF/STO-3G water, end to end | about 197 ms |

CI runs the same benchmarks on each change and compares them with the stored
baseline on `master`.

## Machine-learned potential

The radial model in `tpt-chem-md::mlp` is trained on classical Lennard-Jones
energies of four-atom clusters. Full-batch gradient descent over 3000 epochs
cuts the training MSE from 28.2 to 0.096 kJ²·mol⁻². This shows that the
training loop works. It does not show transferability, because the training
data is synthetic and the model has one species.

## Quantum: cc-pVDZ with spherical d functions

| Check | Result | Test |
| --- | --- | --- |
| Basis size, water | 24 functions (O: 3s2p1d = 14; each H: 2s1p = 5) | `water_has_24_functions` |
| Rotational invariance, water | energy unchanged to 1e-8 Eh under a rotation | `water_energy_is_rotation_invariant` |
| Orthonormal spherical d shell | overlap matrix equals the identity to 1e-10 | `spherical_d_functions_are_orthonormal` |
| Addition theorem, d shell | Σ_m φ_m² depends only on the radius | `spherical_d_density_is_rotation_invariant` |

The water energy at the experimental geometry, −76.0267987 Eh, is recorded as
a regression value. It is not checked against a reference program.

## Independent cross-checks

These compare against other programs or measured data, not against values
recorded in this repository.

| Reference | What is compared | Result | Test |
| --- | --- | --- | --- |
| OpenMM 8.6 (no cutoff) | isolated six-atom cluster, LJ + Coulomb | relative agreement 1e-6 | `isolated_lj_and_coulomb_match_openmm_no_cutoff` |
| OpenMM 8.6 PME | periodic eight-ion box, Coulomb | Ewald within 1e-3 kJ/mol | `periodic_coulomb_matches_openmm_pme` |
| OpenMM 8.6, amber14 | N-methylacetamide: energy and forces | within 1e-3 kJ/mol and 1e-3 kJ/mol/Å | `tests/amber_openmm.rs` |
| CRC Handbook (97th ed.) | mean activity coefficients of NaCl, KCl, CaCl₂, MgCl₂ | worst deviations 0.003, 0.003, 0.0014, 0.007 | `activity/pitzer.rs` (`crc_validation`) |
| OpenMM 8.6 DCDFile | DCD files with unit cells: reading OpenMM's output, and our writer's frame records against OpenMM's bytes | exact coordinates and cells; identical frame bytes | `crates/tpt-chem-io/tests/dcd_openmm.rs` |
| spglib 2.8 database | general-position operations of 230 space groups | closure and operation counts for all 230 | `tests/space_groups.rs` |
| ASE calculator (optional, LGPL) | RHF/STO-3G forces against finite differences | agreement to 2e-4 eV/Å | `crates/tpt-chem-py/tests/test_ase_calculator.py` (skips without ASE) |
| OpenMM 8.6, OPLS-AA (foyer parameters) | n-butane: energy and forces | within 1e-3 kJ/mol and 1e-3 kJ/mol/Å | `crates/tpt-chem-md/tests/opls_openmm.rs` |
| Analytic all-pairs forces | Lennard-Jones forces from dual numbers | agreement to 1e-9 | `tpt-chem-md/tests/autodiff_forces.rs` |

A bug was found by the OpenMM comparison. The Ewald heuristic had
α·r_c ≈ 1.77, which left the energy about 0.6% low. The heuristic now targets
1e-8 convergence.

The Pitzer model reproduces the CRC series with a binary-only parameter set.
CaCl₂ is accurate only up to 1 mol/kg, where the deviation is 0.016.
