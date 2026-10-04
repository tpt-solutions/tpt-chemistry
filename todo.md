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
- [ ] Kani proof: SCF iteration density matrices remain positive semi-definite
- [ ] Kani proof: absence of panics in ERI evaluation, including r→0 singularity handling
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
- [ ] Spin multiplicity inputs
- [x] Charge inputs (`Molecule::set_formal_charge`; electron count = Σ Z − q,
      with H₂²⁺ zero-electron and H₃O⁺ ten-electron tests)
- [ ] UHF / ROHF (radicals, ions, open shells)
- [ ] Analytic nuclear gradients
- [ ] Geometry optimisation (uses gradients)
- [ ] Larger basis sets (6-31G, cc-pVDZ), heavier elements
- [ ] Cube output
- [x] Molden orbital output (`tpt-chem-quantum/src/molden.rs`; molden
      x,y,z p-shell reordering, `[Atoms] (AU)`, occupations from the
      formal charge)

**MD (`tpt-chem-md`)**
- [ ] Bonded forces (harmonic bonds, angles, RB dihedrals) in `System`
- [x] Wire Ewald/PME into the force path: `forces::ForceModel`
      (`AllPairs` / `LjPlusEwald` / `LjPlusPme`) drives the `ensemble`
      run loops with one self-consistent energy/force/virial definition
- [ ] Exclusion lists and 1-4 scaling
- [x] Pressure / virial: `virial_all_pairs`, `ewald_virial`, PME
      virial, `ensemble::pressure_bar` (`(2KE/3 + W/3)/V` in bar via
      `KJ_MOL_ANG3_TO_BAR`)
- [x] Energy minimisation: steepest descent with backtracking line
      search (`tpt-chem-md/src/minimiser.rs`); L-BFGS left open
- [ ] Constraints: SHAKE / RATTLE
- [x] Ready-made run loops: `ensemble::run_nve` / `run_nvt` /
      `run_npt` (`NptSettings`), snapshot reports with mean T/E/P
- [ ] Observables: RDF, MSD, diffusion coefficient, VACF/VDOS
- [ ] Force-field parameter library (OPLS/AMBER/UFF subset, SPC/TIP3P water)
- [ ] Topology builder: bond perception from geometry, atom-type assignment

**Kinetics (`tpt-chem-kinetics`)**
- [ ] Tau-leaping and next-reaction-method SSA
- [ ] Parameter fitting and sensitivity analysis

**Crystal (`tpt-chem-crystal`)**
- [ ] Full 230 space groups (generate from Hall symbols)
- [ ] CIF → structure → XRD pipeline helper
- [ ] Supercell and slab builders

**I/O (`tpt-chem-io`)**
- [ ] `.gro`, SDF, SMILES readers/writers
- [ ] Binary trajectories (DCD / XTC)

### 5.3 Innovation

- [ ] Python bindings (PyO3 + maturin), ASE/MDAnalysis interop
- [ ] WASM build and a browser demo ("run HF in your browser")
- [ ] ASE-style `Calculator` / `Potential` trait (energy + forces) shared by
      MD, HF and ML potentials
- [ ] ML interatomic-potential plug-in (pure-Rust MLP or ONNX) driving the
      same integrators
- [ ] Autodiff hooks for differentiable molecular design (per README promise)
- [ ] Published validation/benchmark table (reference energies, energy-drift
      plots)
- [ ] Run manifest (versions, seed, params) written alongside every output
- [ ] Optional `rayon` feature for ERI shell quartets and force loops

### 5.4 Usability & automation

- [ ] `tpt-chem` CLI: `hf`, `md`, `xrd`, `kinetics` subcommands with
      TOML/JSON inputs
- [ ] Feature-gated `serde` support for `Molecule`, `System`, result structs
- [ ] Per-crate `prelude` modules and a `System` builder
- [ ] CI: criterion benchmark job with regression alerts
- [ ] CI: `cargo doc -D warnings`, MSRV check, `cargo semver-checks`,
      Dependabot, `cargo llvm-cov` coverage
- [ ] Release automation (`release-plz` / `cargo-release`) and lifting the
      "publishing out of scope" decision
- [ ] `justfile` targets: `doc`, `examples`, `bench`

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
- [ ] Starter template (`cargo generate` / GitHub template repo) that loads
      an XYZ and runs a job
- [ ] `data/` folder of sample structures (water, benzene, argon box, NaCl,
      silicon in XYZ/PDB/CIF)
- [ ] mdBook tutorial: 5-minute quick start, "which crate do I need?" table,
      units cheat-sheet (Å/kJ·mol⁻¹/fs vs Bohr/Hartree), validation page
- [ ] Jupyter notebooks (Binder/Colab) once Python bindings exist
- [ ] README upgrade: badges, status/accuracy table, feature matrix vs
      LAMMPS/PySCF/GROMACS
- [ ] "Common mistakes" doc section (Å vs Bohr mix-ups, etc.)
