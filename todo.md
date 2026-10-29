# tpt-chemistry — Build Todo

> Tracks the full 8-crate build-out for the tpt-chemistry foundation repo, per
> `spec.txt`. License for every crate: `MIT OR Apache-2.0`. Author: TPT
> Solutions. **Dependency policy is stricter than the license grant**: the
> workspace itself is dual-licensed, but per spec.txt §2 ("Strict Licensing")
> no dependency that is Apache-2.0-**only** (i.e. not also MIT-licensed) may
> enter the dependency graph — only MIT / dual MIT-Apache-2.0 / other
> MIT-compatible permissive licenses are allowed. This is deliberately
> tighter than sibling repo `tpt-math`'s `deny.toml`, which allows bare
> `Apache-2.0`.
>
> `tpt-math` and `tpt-dsp` are sibling repos in this workspace and may be
> wired in as crates.io deps, GitHub deps, or local `path` deps — whichever
> is available/convenient at the time each crate is scaffolded.

## Phase 0 — Repo Bootstrap

(one-time)

- [x] Root `Cargo.toml` workspace manifest (`resolver = "2"`,
      `[workspace.package]`: `edition = "2021"`, `rust-version` pinned,
      `license = "MIT OR Apache-2.0"`, `authors = ["TPT Solutions"]`)
- [x] `rust-toolchain.toml`
- [x] `rustfmt.toml`
- [x] `deny.toml` — `[licenses] allow` list restricted to MIT-compatible
      licenses only (MIT, BSD-2-Clause, BSD-3-Clause, ISC, Zlib, CC0-1.0,
      Unicode-3.0, MPL-2.0, etc.) — **do not add bare `Apache-2.0`** to the
      allow-list; only `Apache-2.0` combined with MIT in a dual-license
      expression is acceptable, per spec.txt §2
- [x] CI workflow (fmt / clippy / test / deny on push)
- [x] `justfile` (common dev commands)
- [x] `LICENSE-MIT` + `LICENSE-APACHE`
- [x] `CONTRIBUTING.md`
- [x] Root `README.md` stub — tpt-chemistry's role (pure-Rust, AI-native
      comp-chem: MD, quantum chemistry, kinetics, crystallography); link to
      `spec.txt`
- [x] Rust `.gitignore` (`/target`, etc.)
- [x] Create empty `crates/` directory
- [x] `git init` (local only — no GitHub remote/push)
- [x] Initial commit
- [x] Sanity check: `cargo build` succeeds on the empty workspace

## Per-Crate Checklist Template

Every crate below repeats this shape. The umbrella crate (`tpt-chemistry`)
uses the umbrella variant instead of steps 2-4.

**Standard crate:**
1. Scaffold `crates/<name>/` (Cargo.toml inheriting workspace fields, `lib.rs` stub)
2. Wire dependencies (internal `tpt-chem-*` + `tpt-math`/`tpt-dsp` as needed
   + external crates — all vetted MIT-compatible, no Apache-2.0-only)
3. Implement scope (see feature list under each crate below)
4. Unit tests + doctests
5. Rustdoc (crate-level + public API)
6. `cargo fmt --check` / `cargo clippy --all-targets --all-features -- -D warnings` clean
7. `cargo deny check` clean
8. `no_std` target verification — only if crate is `no_std + alloc` per spec.txt §4
9. Crate-specific verification hooks (proptest strategies / Kani harnesses)
   where called out below

**Umbrella crate:**
1. Scaffold `crates/tpt-chemistry/` (Cargo.toml with Cargo features gating each constituent re-export)
2. Wire optional deps + matching feature flags per constituent crate
3. Re-export each constituent's public API behind its feature
4. Rustdoc documenting the feature matrix
5. `cargo fmt` / `clippy` / `deny` clean across feature combinations

---

## Phase 1 — Foundation & MD Engine (Weeks 1-3)

*Per spec.txt §6 "Phase 0". Scaffold workspace, build the core molecular
model and a first working MD engine (Velocity Verlet + LJ/Coulomb), verify
energy conservation.*

### tpt-chem-core

*Foundation layer: molecular topology, atom/bond types, force field params. `no_std + alloc`. Depends on: `tpt-math`.*

- [x] Scaffold `crates/tpt-chem-core/`
- [x] Wire deps: `tpt-math-units` (unit-safe types; uom-based)
- [x] Graph-based molecular representation with const-generic atom types
- [x] Force field definitions: Lennard-Jones, Coulomb, harmonic bonds/angles, Ryckaert-Bellemans dihedrals
- [x] NIST physical constants (CODATA 2022) and unit conversions
- [x] Phantom types distinguishing Cartesian vs. internal (Z-matrix) coordinates (compile-time coordinate-frame safety)
- [x] Unit tests + doctests
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean
- [x] `no_std` verify (`thumbv6m-none-eabi`)

### tpt-chem-io

*File format interop. ⚠️ assumption flagged: spec.txt §6's phase list does not
name this crate explicitly, but XYZ read/write is needed early to load test
geometries for MD and quantum benchmarks, so it's placed here rather than
later. Re-sequence freely if that assumption turns out wrong. Depends on: `tpt-chem-core`.*

- [x] Scaffold `crates/tpt-chem-io/`
- [x] Wire deps: `tpt-chem-core`
- [x] XYZ parser/writer
- [x] PDB parser/writer
- [x] MOL2 parser/writer
- [x] CIF (Crystallographic Information File) parser/writer
- [x] Trajectory streaming for large MD outputs
- [x] Unit tests + doctests (round-trip parse/write fixtures)
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean

### tpt-chem-md (partial — Velocity Verlet + LJ/Coulomb only)

*Classical MD engine. Depends on: `tpt-chem-core`, `tpt-chem-io`.*

- [x] Scaffold `crates/tpt-chem-md/`
- [x] Wire deps: `tpt-chem-core`, `tpt-chem-io`, `tpt-math`
- [x] Velocity Verlet symplectic integrator
- [x] From-scratch cell lists for spatial partitioning
- [x] From-scratch Verlet neighbor lists for O(N) force computation
- [x] Lennard-Jones force evaluation
- [x] Coulomb force evaluation
- [x] Unit tests + doctests
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean

(Leapfrog integrator, thermostats/barostats, and Ewald/PME long-range
electrostatics are deferred to Phase 3 — see `tpt-chem-md` completion below.)

### tpt-chem-verify (partial — MD energy conservation)

*Verification harnesses. Depends on: `tpt-chem-md`, `proptest`.*

- [x] Scaffold `crates/tpt-chem-verify/`
- [x] Wire deps: `tpt-chem-core`, `tpt-chem-md`, `proptest`
- [x] proptest strategies for generating valid molecular geometries and force field parameters
- [x] Property test: Velocity Verlet integration conserves total energy within a bounded tolerance
- [x] Unit tests + doctests
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean

### Phase 1 exit check

- [x] `cargo build` + `cargo test` green across the whole workspace
- [x] Energy-conservation proptest passes reliably (no flakes across N runs) — 6 consecutive clean `tpt-chem-verify` runs

---

## Phase 2 — Quantum & Kinetics (Weeks 4-6)

*Per spec.txt §6 "Phase 1". Electronic structure primitives (HF/SCF) and
reaction kinetics (deterministic + stochastic).*

### tpt-chem-quantum

*Electronic structure: GTO basis sets, SCF, Hartree-Fock, integral evaluation. Depends on: `tpt-chem-core`, `tpt-math`.*

- [x] Scaffold `crates/tpt-chem-quantum/`
- [x] Wire deps: `tpt-chem-core`, `tpt-math` (dense eigenvalue decomposition)
- [x] Gaussian-type orbital (GTO) basis set management (STO-3G for H-Ne; Pople/Dunning future work)
- [x] From-scratch multi-center ERI evaluation (McMurchie–Davidson Hermite
      expansion; validated against Gauss–Hermite quadrature, closed-form
      four-center (ss|ss), exact permutation symmetries, and H₂O/H₂/He
      literature HF/STO-3G energies to 1e-6)
- [x] SCF iteration + density damping/mixing
- [x] Hartree-Fock (HF) matrix construction (restricted closed-shell RHF)
- [x] Roothaan-Hall equation solving via from-scratch cyclic Jacobi eigensolver
- [x] Unit tests + doctests
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean

### tpt-chem-kinetics

*Reaction kinetics and population dynamics. Depends on: `tpt-chem-core`, `tpt-math`.*

- [x] Scaffold `crates/tpt-chem-kinetics/`
- [x] Wire deps: `tpt-chem-core`, `tpt-math` (stiff ODE solvers)
- [x] Deterministic mass-action kinetics (RK4 + implicit Euler)
- [x] Stochastic Simulation Algorithm (Gillespie direct method, in-house RNG)
- [x] Temperature-dependent rate constants: Arrhenius equation
- [x] Temperature-dependent rate constants: Eyring equation
- [x] Equilibrium and steady-state analysis (damped fixed-point iteration)
- [x] Unit tests + doctests
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean

### Integration test

- [x] H₂ HF/STO-3G energy: -1.1166844 Eh vs literature -1.11675931 Eh (75 µEh)

### tpt-chem-verify (extend — mass conservation)

- [x] proptest: mass is strictly conserved in closed reaction networks (`tpt-chem-verify`)

### Phase 2 exit check

- [x] `cargo build` + `cargo test` green across the whole workspace
- [x] H₂ HF integration test passes against reference benchmark (and H₂O, He at 1e-6)

---

## Phase 3 — Solid State & Advanced Verification (Weeks 7-9)

*Per spec.txt §6 "Phase 2". Crystallography, completing long-range MD
electrostatics, and formal (Kani) verification.*

### tpt-chem-crystal

*Solid-state physics and crystallography. Depends on: `tpt-chem-core`, `tpt-math`.*

- [x] Scaffold `crates/tpt-chem-crystal/`
- [x] Wire deps: `tpt-chem-core`, `tpt-math`
- [x] Bravais lattice classification + construction (cubic/tetragonal/ortho/mono/triclinic/hexagonal)
- [x] Miller indices (d-spacings, Bragg angles, reduction)
- [x] Wigner-Seitz cell construction (half-space + rank-3 vertex filter)
- [x] Space group symmetry operations (P1, P-1, P21, P21/c, P212121, Pna21, P4, P-3, P6/mmm, Fm-3m)
- [x] Wyckoff orbit generation (general + special positions)
- [x] Reciprocal lattice generation (2π and crystallographic conventions)
- [x] Brillouin zone sampling (Monkhorst-Pack + Γ-centered k-grids with weights)
- [x] XRD powder simulation (Cromer-Mann form factors, LP factor, extinctions: Si diamond + NaCl validated)
- [x] Unit tests + doctests
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean

### tpt-chem-md (completion — long-range electrostatics + remaining integrators/thermostats)

*Depends on: `tpt-dsp` (3D FFTs).*

- [x] Wire deps: `tpt-dsp` (3D FFT) — tpt-dsp does not exist; resolved in-house: `tpt-chem-md/src/fft.rs` implements a from-scratch radix-2 complex FFT (validated against naive DFTs), keeping the zero-external-solver philosophy
- [x] Ewald summation (real + reciprocal + self, tin-foil)
- [x] Particle Mesh Ewald (PME): `tpt-chem-md/src/pme.rs` — order-4 B-spline charge assignment + mesh reciprocal energy/forces with the Essmann |b|² aliasing correction; O(N log N) via the in-house FFT; validated below
- [x] Leapfrog integrator (validated against velocity Verlet)
- [x] Langevin thermostat (fluctuation-dissipation in MD units)
- [x] Nosé-Hoover thermostat (single chain node)
- [x] Berendsen thermostat/barostat (weak coupling)
- [x] Unit tests + doctests for new integrators/thermostats/PME
- [x] Rustdoc updates
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean

### tpt-chem-verify (completion — Kani proofs)

*Depends on: `kani` (dev-only verification tool).*

- [x] Kani harnesses gated behind cfg(kani) (cell-index bounds, Ewald indices, Verlet skin); CI job on Linux
- [x] Kani proof: neighbor-list cell indexing is always in-bounds
- [x] Kani proof: SCF iteration density matrices remain positive semi-definite
- [x] Kani proof: absence of panics in ERI evaluation, including r→0 singularity handling
- [x] proptest: Velocity Verlet conserves phase-space volume (harmonic-oscillator Jacobian determinant = 1)
- [x] Unit tests + doctests
- [x] Rustdoc
- [x] `cargo fmt` / `clippy` clean
- [x] `cargo deny check` clean

### Phase 3 exit check

- [x] `cargo build` + `cargo test` green across the whole workspace
- [ ] All Kani harnesses pass (bounded model checking green) — needs the Linux CI runner
- [x] PME long-range electrostatics validated against Ewald direct-sum reference: matching wavenumber cutoffs agree to 5.4e-5 relative (32³) and forces to 4e-2 kJ·mol⁻¹·Å⁻¹, with clean 4th-order convergence (error ÷16 per mesh doubling, asserted in tests); PME forces verified as the exact gradient of the mesh energy by central finite differences

---

## Phase 4 — Umbrella Crate & Workspace Polish

*Not explicitly phased in spec.txt §6 — added as a sensible closing phase,
matching how sibling repos (e.g. `tpt-math`) finish with an umbrella
re-export crate and a full workspace hygiene pass.*

### tpt-chemistry (umbrella crate)

*Feature-gated umbrella re-exporting all `tpt-chem-*` crates.*

- [x] Scaffold `crates/tpt-chemistry/`
- [x] Wire optional deps + matching Cargo features per constituent crate (`core`, `md`, `quantum`, `kinetics`, `crystal`, `io`, `verify`)
- [x] Re-export each constituent's public API behind its feature
- [x] Rustdoc documenting the feature matrix
- [x] `cargo fmt` / `clippy` / `deny` clean across feature combinations (core/md/quantum subsets + all-features tested)

### Workspace-wide polish

- [x] Workspace-wide `cargo deny check` pass confirming zero Apache-2.0-only dependencies anywhere in the graph
- [x] README pass covering all 8 crates
- [x] README: integration/synergy notes per spec.txt §5
- [x] Publishing explicitly out-of-scope for this pass (local git only, matching tpt-math's initial approach)

## Phase 6 — Aqueous Thermodynamics (post-Phase-5 addition)

*Activity coefficients, solubility products and chelator-controlled calcium
phosphate precipitation. Depends on: none (own constants; `libm` not needed).
Parameter defaults for phosphate and ligand constants are illustrative and
must be replaced with measured data.*

### tpt-chem-aqueous

- [x] Scaffold `crates/tpt-chem-aqueous/` and workspace/umbrella wiring (`aqueous` feature)
- [x] Constants, water density/permittivity, Debye–Hückel A_φ / A_γ / B(T)
- [x] Species, ionic strength and charge-balance helpers
- [x] Activity models: limiting law, extended Debye–Hückel, Davies
- [x] Pitzer ion-interaction model with osmotic coefficient (NaCl, KCl, CaCl₂ validated; Gibbs–Duhem checked)
- [x] eNRTL model (dilute limit, Gibbs–Duhem residual test) — **no validated parameter sets yet**
- [x] Van't Hoff equilibrium constants and solubility products / saturation indices
- [x] H₃PO₄ / Ca–phosphate constants and candidate phases (brushite, ACP, OCP; HAP opt-in)
- [x] Multi-site Ca²⁺ chelator binding (casein-phosphopeptide-like defaults)
- [x] Coupled speciation with activity-coefficient iteration, pH × T scans, onset-pH search
- [x] Unit tests, doctests, rustdoc; `cargo fmt` / `clippy -D warnings` / rustdoc clean (standalone build)
- [x] `tpt-chem-verify` proptest properties (DH monotonicity, Pitzer finiteness, van't Hoff, Ca/P balances, chelator monotonicity)
- [x] `cpp_calcium_phosphate` example (umbrella crate, `aqueous` feature)
- [x] README crate and example table rows
- [x] Build and test inside the full workspace (needs `../tpt-math` checkout; builds and passes with `tpt-math` cloned as a sibling)
- [x] `cargo deny check` clean with the new crate
- [ ] Validate Pitzer/eNRTL and phosphate constants against measured data; document sources.
      Pitzer: NaCl, KCl, CaCl₂ and MgCl₂ mean activity coefficients are checked against the CRC
      Handbook measured series (0.1–5 mol/kg for NaCl and KCl, 0.1–1 mol/kg for CaCl₂ and MgCl₂).
      Worst deviations are 0.003, 0.003, 0.0014 and 0.007. CaCl₂ is not accurate above 1 mol/kg
      with a binary-only model (0.016 at 2 mol/kg). Tests are in `activity/pitzer.rs`. Still open:
      eNRTL and the phosphate constants, since no measured data source is available for them
- [x] Nucleation/kinetic onset model: `NucleationKinetics` (Söhnel–Garside `J = A·exp(−B/ln²S)`) and `CalciumPhosphateSystem::kinetic_onset_ph`; defaults are illustrative and need calibration against induction-time data
- [x] Add an entry for `tpt-chem-aqueous` to `spec.txt` §4

---

## Known open items after this pass

- [x] Water HF/STO-3G energy: RESOLVED. The electron-repulsion integrals
      were rewritten on the McMurchie–Davidson Hermite expansion (the OS
      `raise_family` machinery is gone): per-axis E-coefficients for both
      primitive pairs contracted against the two-electron Hermite Coulomb
      auxiliary `R^n_{tuv}(ρ, P−Q)` at the reduced exponent `ρ = pq/(p+q)`
      with the `(−1)^{τ+ν+φ}` electron-2 sign. Validated against an
      independent 3D Gauss–Hermite quadrature reference (100+ values, all
      component classes, ≤ 2.3e-8 = the reference's own finite-difference
      noise), the closed-form four-center (ss|ss) formula, and the exact
      8-fold permutation symmetries. `water_sto3g_near_equilibrium` is
      un-ignored: E = −74.96590117 Eₕ, the canonical RHF/STO-3G water
      energy, reproduced to 1e-6 and cross-checked by an independent
      Python RHF (orbitals −20.2516/−1.2576/−0.5939/−0.4597/−0.3926).
      H₂ = −1.11668439 Eₕ and He = −2.80778516 Eₕ (Basis Set Exchange He
      exponents) tightened to 1e-6 tolerances.
- [x] PME with 3D FFTs: RESOLVED in-house. `tpt-chem-md/src/fft.rs`
      (radix-2 Cooley–Tukey, forward/inverse, 3D row-column field) and
      `tpt-chem-md/src/pme.rs` (order-4 B-spline assignment, Essmann |b|²
      aliasing correction, spline-derivative forces). Validated against the
      direct Ewald sum: matched-cutoff agreement at the 4th-order
      discretization level, exact mesh-energy gradient forces, neutral-force
      sum. `tpt-dsp` can later replace the FFT backend without API change.
- [ ] Kani proofs run only on Linux CI (no Windows host support); harnesses
      are written and compile-gated behind `#[cfg(kani)]`.
- [x] PME-vs-Ewald validation test (Phase 3 exit): `pme_matches_ewald_with_matching_wavenumber_cutoff` and `pme_converges_to_ewald_as_mesh_refines` in `crates/tpt-chem-md/src/pme.rs`.

---

## Phase 5 — Platform Review Follow-ups (2026-10-05)

*Outcome of the platform review: bugs, missing features, innovation,
usability/automation, and adoption work. Ordered roughly by priority.*

### 5.1 Blockers & correctness

- [x] CI cannot build: path deps on `../tpt-math/...` — every CI job now
      checks out `tpt-solutions/tpt-math` into a sibling directory
- [x] Document the sibling-repo prerequisite (`../tpt-math`) in README and
      CONTRIBUTING so fresh clones build
- [x] Tick the Phase 1/2/3 exit-check boxes now that build + tests pass;
      proptest flake check: 6 consecutive clean runs
- [x] Kani proof: SCF density matrices remain positive semi-definite
      (`scf_density_is_psd` — the density construction is extracted into
      `density_from_orbitals` so the harness proves the actual code)
- [x] Kani proof: no panics in ERI evaluation (r→0 handling)
      (`eri_r_zero_no_panic`, `eri_p_shell_coincident_no_panic`)
- [ ] Run and confirm all Kani harnesses pass on Linux CI
- [x] `rhf`: dedicated `HfError::BasisTooLarge` when the basis exceeds
      `MAX_BASIS` (tested)
- [x] `rhf`: removed the `f64::from(n_occ as u8)` truncation
- [x] `rhf`: dense `nb⁴` ERI array replaced with Schwarz-screened packed
      shell quartets (`accumulate_g`); memory is O(N²) per SCF plus the
      screened quartets; `MAX_BASIS` raised 128 → 256
- [x] MD: pair-force kernel de-duplicated (`pair_energy_force`)
- [x] MD: `berendsen_barostat` invalidates the Verlet list (new
      `VerletList::invalidate` + `stale` flag; `needs_update` honors it;
      tested)
- [x] Audit non-test `unwrap`/`expect`/`panic!` in library code (actual
      count ~10 after excluding doctests/tests): NaN-sensitive
      `partial_cmp().unwrap()` sorts in `xrd` and `wigner` switched to
      `f64::total_cmp`; the `bravais` zero-volume reciprocal panic is now a
      documented `# Panics` contract on both public constructors; the
      remaining three (`Vec3` index-out-of-range, `hf` `shell_of`,
      `ssa` initial-state `expect`s) are structural invariants with
      explanatory messages; input-dependent failures already return
      `Result` (`io` parsers, `hf`, `basis`, `linalg`).
- [x] Fixed stale docs: README and `tpt-chem-quantum` Cargo/lib docs now
      say McMurchie–Davidson; the "DFT groundwork" claim was removed
- [x] Added `[profile.bench]` (thin LTO, symbols kept, no overflow checks)

### 5.2 Missing features

**Quantum (`tpt-chem-quantum`)**
- [x] Spin multiplicity inputs (`uhf(mol, multiplicity)`)
- [x] Charge inputs (`Molecule::set_formal_charge`; electron count = Σ Z − q,
      with H₂²⁺ zero-electron and H₃O⁺ ten-electron tests)
- [x] UHF (`tpt-chem-quantum/src/uhf.rs`: separate α/β Fock, ⟨S²⟩, tests vs
      H-atom reference, RHF singlet equivalence, H₂ triplet)
- [x] ROHF (`tpt-chem-quantum/src/rohf.rs`: one shared orbital set with Roothaan's effective
      Fock matrix; reduces to RHF for closed shells and to UHF for one open-shell electron; tested
      that E_ROHF ≥ E_UHF for a triplet and an OH radical. Exposed as `tpt-chem hf --rohf`)
- [x] Analytic nuclear gradients (`tpt-chem-quantum/src/gradient.rs`:
      `rhf_gradient`, Ha/Bohr; derivative integrals via l±1 raw primitive
      shells over the existing McMurchie–Davidson code; validated against
      finite differences for H₂ and water; 
- [x] Geometry optimisation (`gradient::optimize_geometry`, BFGS +
      backtracking; H₂ and water reproduce HF/STO-3G minima; the water test
      is `#[ignore]`d in debug builds — run with `--release -- --ignored`)
- [ ] Larger basis sets: 6-31G for H–Ne done (`BasisSet::Pople631g`, data
      from BSE in `basis_631g.rs`; H₂/water energies and gradients
      validated; `*_with_basis` entry points, CLI `--basis`). cc-pVDZ for H–Ar is added (`BasisSet::CcPvdz`, data from BSE via `tools/gen_cc_pvdz.py`).
      Spherical d shells go through a Cartesian-to-spherical transform
      (`tpt-chem-quantum/src/spherical.rs`): the functions are orthonormal, and the energy is
      rotation invariant to 1e-8. Gradients also run for spherical bases: the density is pulled
      back to Cartesian components (`to_cartesian_space` in `gradient.rs`), checked against finite
      differences for cc-pVDZ water. Molden writes `[5D]` for pure d shells (spherical coefficients
      in molden order) and `[6D]` for Cartesian d; the old `[5D7F]` flag was wrong for Cartesian
      bases. Cube output is not validated for spherical shells. 6-31G for Na–Ar is added (`BasisSet::Pople631g`,
      data from BSE in `data/631g-na-ar.txt`, `tools/gen_631g_na_ar.py`). Its Ar atom SCF converges to
      −526.77215 Eₕ, inside a sanity window (above the HF limit, about −526.82). No literature energy
      has been checked against it, so it is not yet validated. Still open: elements beyond Ar
- [x] Cube output (`tpt-chem-quantum/src/cube.rs`: density and MO values on a grid, Gaussian `.cube` writer; density integrates to N_e and MOs to 1 in tests; `water_cube` example. Grid spacing ~0.1 Bohr is needed for the O 1s core)
- [x] Molden orbital output (`tpt-chem-quantum/src/molden.rs`; molden
      x,y,z p-shell reordering, `[Atoms] (AU)`, occupations from the
      formal charge)

**MD (`tpt-chem-md`)**
- [x] Bonded forces: `System::bonded` topology (harmonic bonds, angles,
      RB dihedrals over the core kernels), folded into every
      `ForceModel::evaluate`/`virial`; finite-difference gradient test.
      The pass also fixed a virial sign error in `ewald_virial`/`pme_virial`
      (W = Σ r·F, not −Σ r·F)
- [x] Wire Ewald/PME into the force path: `forces::ForceModel`
      (`AllPairs` / `LjPlusEwald` / `LjPlusPme`) drives the `ensemble`
      run loops with one self-consistent energy/force/virial definition
- [x] Exclusion lists and 1-4 scaling (`tpt-chem-md/src/exclusions.rs`,
      `System::exclusions`, `Exclusions::from_bonded`; honoured by the all-pairs,
      neighbor-list, LJ-long-range and virial paths; Ewald/PME subtract the
      scaled pairs' missing Coulomb). The pass also fixed two real bugs:
      the Ewald real-space force had the wrong sign (invisible to symmetric
      lattice tests; now FD-regression-tested) and `ForceModel::LjPlusPme`
      returned only the mesh part, omitting real-space and self terms
- [x] Pressure / virial: `virial_all_pairs`, `ewald_virial`, PME
      virial, `ensemble::pressure_bar` (`(2KE/3 + W/3)/V` in bar via
      `KJ_MOL_ANG3_TO_BAR`)
- [x] Energy minimisation: steepest descent with backtracking line
      search (`tpt-chem-md/src/minimiser.rs`); L-BFGS added as `minimise_lbfgs`
      (two-loop recursion, Armijo backtracking, restart from steepest descent on
      line-search failure; tested on a dimer and a smooth 27-atom cluster, where
      it needs about 6× fewer steps. Periodic lattices with an unshifted LJ cutoff
      are not smooth enough for either method, so they are not used as the benchmark)
- [x] Constraints: SHAKE / RATTLE (`tpt-chem-md/src/constraints.rs`;
      `Constraints::step` = constrained velocity Verlet over any `ForceModel`;
      rigidity + energy-drift tests; DOF correction left to the caller)
- [x] Ready-made run loops: `ensemble::run_nve` / `run_nvt` /
      `run_npt` (`NptSettings`), snapshot reports with mean T/E/P
- [x] Observables: RDF (ideal-gas normalized, lattice-shell test), MSD
      (time-origin averaged), VACF, Einstein and Green–Kubo diffusion
      (`tpt-chem-md/src/observables.rs`; VDOS via FFT of the VACF left as
      a user step)
- [ ] Force-field parameter library: an AMBER-style subset is added (`tpt-chem-md/src/amber.rs`:
      harmonic bonds and angles, periodic torsions including impropers, LJ and Coulomb with 1-4
      scaling, forces by autodiff). Energy and forces match OpenMM amber14 on N-methylacetamide to
      1e-3 kJ/mol (`tests/amber_openmm.rs`, fixture from `tools/amber_reference.py`). Atom types are
      assigned by hand. OPLS-AA is validated against OpenMM for seven molecules (`tests/opls_openmm.rs`): butane,
      ethanol, 2-propanol, acetone, dimethyl ether, benzene and acetic acid. Parameters come from foyer
      0.7.6's `oplsaa.xml`, cross-checked against GROMACS's `oplsaa.ff` for the CT and HC types. Energies
      and forces agree with OpenMM to 1e-3. The set covers alkyl, hydroxyl, ketone, ether, aromatic and
      carboxyl types. Benzene has no proper torsions, which OPLS-AA also defines as none. Geometries are
      hand-placed and then relaxed with OpenMM's minimiser (residual forces about 1e-6 kJ/mol/Å).
      Atom types are hand-assigned, with no SMARTS typing, and UFF is still open. SPC and TIP3P rigid water done
      (`tpt-chem-md/src/water.rs`: `water_box` wires charges, intramolecular
      exclusions and SHAKE constraints; `examples/water_pme.rs` gives
      −39.7 kJ/mol/molecule PE at 300 K with PME and ~1 kJ/mol NVE drift
      over 500 fs in a 64-molecule box). OPLS beyond butane and UFF subsets still open
- [x] Topology builder (`tpt-chem-md/src/topology.rs`: `build_system` does bond
      perception, geometry-referenced bonds/angles, OPLS-like element LJ,
      exclusions + 1-4 scaling; deliberately generic — no torsions or
      transferable atom typing, which stay with the OPLS/AMBER item). Also
      fixed: `ForceModel::evaluate` and the Verlet integrator ignored
      `System::bonded` (virial alone included it), so flexible molecules
      felt no bonded forces

**Kinetics (`tpt-chem-kinetics`)**
- [x] Tau-leaping (`solver::tau_leap`, Poisson counts with negative-
      population halving guard; means validated against SSA)
- [x] Next-reaction-method SSA (`solver::next_reaction[_events]`, Gibson–Bruck
      with dependency graph and putative-time rescaling; mean-decay test)
- [x] Parameter fitting and sensitivity analysis (`tpt-chem-kinetics/src/fit.rs`: `fit_rates` is Levenberg–Marquardt on ln k with RK4 integration, recovers the rates of A→B→C from noise-free data to 1e-5; `sensitivity` gives normalized ∂y/∂ln k, checked against the analytic decay result. Unweighted. `FitResult::log_std_errors` gives standard errors of ln k from the linearized covariance, and a 60-replicate test checks them against the spread of fitted values)

**Crystal (`tpt-chem-crystal`)**
- [x] Full 230 space groups: general-position operations for all 230 standard settings in
      `data/spacegroups.txt`, generated from the spglib database by `tools/gen_spacegroups.py`.
      The generator reads spglib's table rather than running a Hall-symbol parser. Tests check
      group closure and known operation counts for all 230 (`tests/space_groups.rs`). Wyckoff
      letters are not enumerated, but orbits of special seed positions work for every group.
      Checking the table exposed three wrong hand-coded groups (P2₁, P2₁2₁2₁ and P6/mmm), so
      `by_symbol` now reads the table
- [x] CIF → structure → XRD pipeline helper (`tpt_chemistry::pipeline::{xrd_from_cif, lattice_from_cif, structure_from_cif}`, behind `io`+`crystal`+`core`; no symmetry expansion, so the CIF must list every atom in the cell)
- [x] Supercell and slab builders (`tpt-chem-crystal/src/build.rs`)

**I/O (`tpt-chem-io`)**
- [x] `.gro` reader/writer (`tpt-chem-io/src/gro.rs`, nm↔Å at the API
      boundary, round-trip + fixed-column reference tests)
- [x] SDF (V2000) reader/writer (`tpt-chem-io/src/sdf.rs`; bonds, `M  CHG` read into total charge)
- [x] Gaussian input/output (spec.txt §3 lists "Gaussian input/output parsing" for `tpt-chem-io`;
      it was missing from this todo). `tpt-chem-io/src/gaussian.rs`: `write_gaussian_input` /
      `parse_gaussian_input` (route, title, charge/multiplicity, Cartesian coordinates, Link 0
      lines skipped) and `parse_gaussian_output` (last `SCF Done` energy, last `Standard
      orientation` geometry, normal-termination flag). Output fixtures are hand-written in
      Gaussian's documented layout and have not been checked against a real Gaussian run
- [x] SMILES reader (`tpt-chem-io/src/smiles.rs`, `parse_smiles`: organic subset with implicit H,
      bracket atoms, charges, aromatic lowercase, ring closures incl. `%nn`, branches; graph and
      formula only, atoms at the origin; stereo marks ignored, no kekulization)
- [x] SMILES writer (`write_smiles`: canonical output. Atoms are ranked from element, degree,
      H count and aromaticity, refined by neighbourhood, and ties are broken by keeping the
      smallest string over all candidates. Relabelled copies give identical strings (tested).
      Ring closures, branches and bond symbols; round-trips through `parse_smiles`. Refuses
      net-charged molecules, since per-atom charges are not stored in `Molecule`)
- [x] Binary trajectories, DCD: `tpt-chem-io/src/dcd.rs` reads and writes DCD files, checked in
      against OpenMM's `DCDFile` output (`tests/dcd_openmm.rs`; our writer's frame records match
      OpenMM's bytes).
      The time step is in CHARMM AKMA units. XTC is implemented in `tpt-chem-io/src/xtc.rs` (our own code, MIT OR
      Apache-2.0). The reader and writer use molstar's bit widths (`bitlength(P)` for
      coordinates, `smallidx` for run deltas). Compressed fixtures `ours_*` are written by our writer
      and check exactly against formula positions; they also decode with molstar's MIT parser
      (checked once by a port). Only gas_5 (uncompressed) comes from mdtraj (see
      `tests/data/xtc/README.md`). Still open: a second independent compressed-stream oracle, and
      the mdtraj compressed frames that our reader still cannot decode (an unverified cause in the
      fixture, not the reader)

### 5.3 Innovation

- [x] Python bindings (PyO3 + maturin): `crates/tpt-chem-py` (outside the workspace, since the
      extension links libpython). `rhf_energy_xyz`, `formula_from_smiles`, `canonical_smiles`;
      pytest suite passes against a release wheel
- [x] ASE interop: an ASE calculator for RHF/STO-3G energies and analytic forces
      (`crates/tpt-chem-py/ase_calculator.py`; forces checked against finite differences). ASE is
      LGPL-2.1, so it is optional: the test skips without it. See `docs/third-party.md`
- [x] WASM build and a browser demo ("run HF in your browser"): `crates/tpt-chem-wasm`
      (wasm-pack `--target web`, `www/index.html`). Node smoke test `tests/smoke.mjs` passes
      and matches the native energy. The page itself has not been opened in a browser
- [x] ASE-style `Potential` trait (`tpt-chem-core/src/potential.rs`, no_std;
      Å / kJ·mol⁻¹ units) with a generic BFGS `minimize_bfgs`; implemented
      by `RhfPotential` (quantum) and `ModelPotential` (MD); ML potentials
      can plug in the same way
- [x] ML interatomic-potential plug-in (pure-Rust MLP) driving the same integrators:
      `tpt-chem-md/src/mlp.rs` (radial Gaussian symmetry functions, tanh network, analytic
      forces checked against finite differences, gradient-descent training, NVE via
      `VelocityVerlet::step_with`). Single species, radial only, no ONNX import. Training
      cut the MSE on LJ cluster energies by about 300×
- [x] Autodiff hooks for differentiable molecular design (per README promise): forward-mode dual
      numbers in the `no_std` core (`tpt-chem-core/src/autodiff.rs`). The AMBER model's forces come
      from them, and the generic Lennard-Jones gradient matches the analytic all-pairs forces to
      1e-9 (`tpt-chem-md/tests/autodiff_forces.rs`). Reverse mode is in `tpt-chem-core/src/autodiff_reverse.rs` (tape-based, `std` feature only: one
      forward and one backward pass for any number of inputs; checked against forward mode on a
      200-input function). `AmberModel::energy_forces` now uses it (one backward pass for all forces), with forward mode kept as
      `energy_forces_forward` and cross-checked on all OPLS molecules. In the MD `ForceModel`, bonded terms also have a reverse-mode path (`bonded::bonded_energy_forces_ad`,
      cross-checked against the analytic kernels). It is not the default: it measured about 19× slower on a
      3000-term chain (`bonded_timing`), so the analytic kernels stay in the MD path. The nonbonded pair
      kernels stay hand-written, since a tape per pair would be slower still
- [x] Independent checks against OpenMM: isolated LJ + Coulomb (relative 1e-6), periodic Ewald
      against OpenMM PME (1e-3 kJ/mol), and the AMBER model on N-methylacetamide (1e-3 kJ/mol).
      These found a bug: the Ewald heuristic was under-converged (α·r_c ≈ 1.77, about 0.6% low in
      energy). It now targets 1e-8 from σ = √(−ln 10⁻⁸)
- [x] Validation and benchmark table in the mdBook (`book/src/validation.md`), with each
      number tied to the test that enforces it. Energy-drift plots are not included; external
      publication is out of scope
- [x] Run manifest (versions, seed, params) written alongside every output: `tpt_chemistry::manifest::RunManifest` (JSON with crate versions, enabled features, seed, params, inputs, FNV-1a output hashes; no timestamp so reruns match). The `tpt-chem` CLI writes `<out>.manifest.json` for `hf --opt --out` and `md --out`. The `water_hf` and `water_cube` examples also write manifests. Only those two examples write files, so the rest need none. The CLI xrd and kinetics commands print to stdout and write no file
- [x] Optional `rayon` feature for ERI shell quartets and force loops. The
      `rayon` feature (opt-in; forwarded by the umbrella) parallelizes the MD
      pair loops (`forces_all_pairs`, `forces_neighbor_list`), the quartet
      integral evaluation in RHF `scf_setup`, and the Fock build (`accumulate_g`).
      Measured on benzene RHF/STO-3G (36 basis functions, release, one run
      each): serial 122.0 s, rayon 39.0 s (about 3.1×), same energy
      −227.89100277 Eh. An earlier Fock-only parallel build was slower (about
      140 s vs 118 s) because the integrals dominated; parallelizing the
      integrals is what made the difference

### 5.4 Usability & automation

- [x] `tpt-chem` CLI (`tpt-chemistry/src/bin/tpt-chem.rs`): `hf` (RHF/UHF,
      `--basis`, `--grad`, `--opt`), `xrd`, `md` (rigid water box with PME, or
      an XYZ via the topology builder, optional minimise/trajectory), and
      `kinetics` (plain-text network → ODE or SSA CSV). TOML/JSON inputs
      left open
- [x] Feature-gated `serde` support for `Molecule`, `System`, result structs. Behind a `serde` feature (forwarded by the umbrella crate): `Molecule` (rebuilt through its constructors on load, so invalid input is an error), `Vec3`, `Atom`, `Bond`, the MD `System` with its topology and exclusions (exclusions serialize as a sorted list, so output is deterministic), and `HfResult`. Also the kinetics networks, rate laws and fit results (`fit.rs`), and crystal lattices and Bravais types. Note: serde_json needs its `float_roundtrip` feature for bit-exact f64 round trips
- [x] Per-crate `prelude` modules and a `System` builder (`tpt_chem_md::system::SystemBuilder` validates topology and exclusion indices and can seed thermal velocities; md, quantum and kinetics have preludes; core already had one)
- [x] CI: criterion benchmark job (`tpt-chem-md` force kernels, `tpt-chem-quantum` RHF water).
      Alerts use github-action-benchmark at 150% and do not fail the build yet. The CI job itself
      has not run on GitHub
- [x] CI: MSRV check. `rust-version` raised to 1.88 (the sibling `tpt-math` crates need 1.85 and
      `proptest 1.12` needs 1.88). `cargo +1.88 check --workspace --all-targets` passes locally;
      the `msrv` job in `ci.yml` has not run on GitHub yet
- [x] CI: Dependabot (`.github/dependabot.yml`, cargo and github-actions, weekly)
- [ ] CI: `cargo semver-checks` (needs a published baseline, so it is part of publishing)
- [x] CI: `cargo llvm-cov` coverage job (`coverage` in `ci.yml`; LCOV uploaded as an artifact,
      informational only. Not yet run on GitHub, and no coverage floor set)
- [ ] Release automation (`release-plz` / `cargo-release`) and lifting the
      "publishing out of scope" decision
- [x] `justfile` targets: `doc` (warnings denied), `examples` (all 11 umbrella examples, verified to exit 0), `bench` (runs nothing until the criterion job exists)

### 5.5 Adoption: examples, templates, docs

- [x] `examples/` in the umbrella crate, each runnable via
      `cargo run --example <name>` (feature-gated via `required-features`):
  - [x] `h2_hf` (bond-length scan / PES)
  - [x] `water_hf` (orbitals, energies, Molden file)
  - [x] `lj_argon_nve` (NVT equilibration → NVE drift)
  - [x] `nacl_pme` (Ewald vs PME)
  - [x] `lotka_volterra_ssa` (logistic LV — the neutral LV system is
        SSA-unbounded, noted in the example) and `michaelis_menten`
        (stochastic vs deterministic)
  - [x] `silicon_xrd` / `nacl_xrd` (CSV pattern output)
  - [x] `cif_to_xrd` (end-to-end file I/O pipeline)
- [x] Starter template (`templates/starter`, a cargo-generate template). It reads an XYZ and
      runs RHF/STO-3G. Verified by rendering it and by `cargo generate --path`, then building
      and running the output. A GitHub template repository is not created
- [x] `data/` folder of sample structures (water, benzene, argon box, NaCl,
      silicon). XYZ for all five; PDB for the three molecules; CIF (P1) for
      argon and silicon, since a CIF is only meaningful for crystals. All
      parsed and geometry-checked by `tpt-chem-io/tests/data_samples.rs`
- [x] mdBook tutorial (`book/`, builds with `mdbook build book`): quick start, "which crate do
      I need?" table, units cheat-sheet, common mistakes, validation, bindings
- [ ] Jupyter notebooks (Binder/Colab): `notebooks/tpt_chemistry_quickstart.ipynb` exists and its
      code cells run against the wheel. Binder/Colab configuration is not done
- [ ] README upgrade: badges, status and limits, and interfaces are done. A quantitative feature
      matrix against LAMMPS/PySCF/GROMACS is still open, because it needs measured comparisons
- [x] "Common mistakes" doc section (`docs/common-mistakes.md`: units by crate, Å vs Bohr, Hartree vs kJ·mol⁻¹, molality, activity pH)
