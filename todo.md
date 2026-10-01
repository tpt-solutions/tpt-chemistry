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

- [ ] Scaffold `crates/tpt-chem-core/`
- [ ] Wire deps: `tpt-math` (linear algebra), `tpt-math-units` (unit-safe types)
- [ ] Graph-based molecular representation with const-generic atom types
- [ ] Force field definitions: Lennard-Jones, Coulomb, harmonic bonds/angles, Ryckaert-Bellemans dihedrals
- [ ] NIST physical constants and unit conversions
- [ ] Phantom types distinguishing Cartesian vs. internal (Z-matrix) coordinates (compile-time coordinate-frame safety)
- [ ] Unit tests + doctests
- [ ] Rustdoc
- [ ] `cargo fmt` / `clippy` clean
- [ ] `cargo deny check` clean
- [ ] `no_std` verify (target TBD, e.g. `thumbv6m-none-eabi`)

### tpt-chem-io

*File format interop. ⚠️ assumption flagged: spec.txt §6's phase list does not
name this crate explicitly, but XYZ read/write is needed early to load test
geometries for MD and quantum benchmarks, so it's placed here rather than
later. Re-sequence freely if that assumption turns out wrong. Depends on: `tpt-chem-core`.*

- [ ] Scaffold `crates/tpt-chem-io/`
- [ ] Wire deps: `tpt-chem-core`
- [ ] XYZ parser/writer
- [ ] PDB parser/writer
- [ ] MOL2 parser/writer
- [ ] CIF (Crystallographic Information File) parser/writer
- [ ] Trajectory streaming for large MD outputs
- [ ] Unit tests + doctests (round-trip parse/write fixtures)
- [ ] Rustdoc
- [ ] `cargo fmt` / `clippy` clean
- [ ] `cargo deny check` clean

### tpt-chem-md (partial — Velocity Verlet + LJ/Coulomb only)

*Classical MD engine. Depends on: `tpt-chem-core`, `tpt-chem-io`.*

- [ ] Scaffold `crates/tpt-chem-md/`
- [ ] Wire deps: `tpt-chem-core`, `tpt-chem-io`, `tpt-math`
- [ ] Velocity Verlet symplectic integrator
- [ ] From-scratch cell lists for spatial partitioning
- [ ] From-scratch Verlet neighbor lists for O(N) force computation
- [ ] Lennard-Jones force evaluation
- [ ] Coulomb force evaluation
- [ ] Unit tests + doctests
- [ ] Rustdoc
- [ ] `cargo fmt` / `clippy` clean
- [ ] `cargo deny check` clean

(Leapfrog integrator, thermostats/barostats, and Ewald/PME long-range
electrostatics are deferred to Phase 3 — see `tpt-chem-md` completion below.)

### tpt-chem-verify (partial — MD energy conservation)

*Verification harnesses. Depends on: `tpt-chem-md`, `proptest`.*

- [ ] Scaffold `crates/tpt-chem-verify/`
- [ ] Wire deps: `tpt-chem-core`, `tpt-chem-md`, `proptest`
- [ ] proptest strategies for generating valid molecular geometries and force field parameters
- [ ] Property test: Velocity Verlet integration conserves total energy within a bounded tolerance
- [ ] Unit tests + doctests
- [ ] Rustdoc
- [ ] `cargo fmt` / `clippy` clean
- [ ] `cargo deny check` clean

### Phase 1 exit check

- [ ] `cargo build` + `cargo test` green across the whole workspace
- [ ] Energy-conservation proptest passes reliably (no flakes across N runs)

---

## Phase 2 — Quantum & Kinetics (Weeks 4-6)

*Per spec.txt §6 "Phase 1". Electronic structure primitives (HF/SCF) and
reaction kinetics (deterministic + stochastic).*

### tpt-chem-quantum

*Electronic structure: GTO basis sets, SCF, Hartree-Fock, integral evaluation. Depends on: `tpt-chem-core`, `tpt-math`.*

- [ ] Scaffold `crates/tpt-chem-quantum/`
- [ ] Wire deps: `tpt-chem-core`, `tpt-math` (dense eigenvalue decomposition)
- [ ] Gaussian-type orbital (GTO) basis set management (STO-nG, Pople, Dunning sets)
- [ ] From-scratch multi-center electron repulsion integral (ERI) evaluation (Obara-Saika or McMurchie-Davidson scheme)
- [ ] Self-Consistent Field (SCF) iteration + density mixing
- [ ] Hartree-Fock (HF) matrix construction
- [ ] Roothaan-Hall equation solving via `tpt-math` eigen-solvers
- [ ] Unit tests + doctests
- [ ] Rustdoc
- [ ] `cargo fmt` / `clippy` clean
- [ ] `cargo deny check` clean

### tpt-chem-kinetics

*Reaction kinetics and population dynamics. Depends on: `tpt-chem-core`, `tpt-math`.*

- [ ] Scaffold `crates/tpt-chem-kinetics/`
- [ ] Wire deps: `tpt-chem-core`, `tpt-math` (stiff ODE solvers)
- [ ] Deterministic mass-action kinetics (stiff ODE systems)
- [ ] Stochastic Simulation Algorithm (SSA / Gillespie algorithm) for low-copy-number regimes
- [ ] Temperature-dependent rate constants: Arrhenius equation
- [ ] Temperature-dependent rate constants: Eyring equation
- [ ] Equilibrium and steady-state analysis
- [ ] Unit tests + doctests
- [ ] Rustdoc
- [ ] `cargo fmt` / `clippy` clean
- [ ] `cargo deny check` clean

### Integration test

- [ ] Hartree-Fock energy of H₂ matches analytical/NIST benchmark (within tolerance)

### tpt-chem-verify (extend — mass conservation)

- [ ] proptest: mass is strictly conserved in closed reaction networks (`tpt-chem-kinetics`)

### Phase 2 exit check

- [ ] `cargo build` + `cargo test` green across the whole workspace
- [ ] H₂ HF integration test passes against reference benchmark

---

## Phase 3 — Solid State & Advanced Verification (Weeks 7-9)

*Per spec.txt §6 "Phase 2". Crystallography, completing long-range MD
electrostatics, and formal (Kani) verification.*

### tpt-chem-crystal

*Solid-state physics and crystallography. Depends on: `tpt-chem-core`, `tpt-math`.*

- [ ] Scaffold `crates/tpt-chem-crystal/`
- [ ] Wire deps: `tpt-chem-core`, `tpt-math`
- [ ] Bravais lattice generation
- [ ] Miller indices
- [ ] Wigner-Seitz cell construction
- [ ] Space group symmetry operations
- [ ] Wyckoff positions
- [ ] Reciprocal lattice generation
- [ ] Brillouin zone sampling (k-point grids)
- [ ] X-ray diffraction (XRD) pattern simulation from structure factors
- [ ] Unit tests + doctests
- [ ] Rustdoc
- [ ] `cargo fmt` / `clippy` clean
- [ ] `cargo deny check` clean

### tpt-chem-md (completion — long-range electrostatics + remaining integrators/thermostats)

*Depends on: `tpt-dsp` (3D FFTs).*

- [ ] Wire deps: `tpt-dsp` (3D FFT)
- [ ] Ewald summation for long-range electrostatics
- [ ] Particle Mesh Ewald (PME) using `tpt-dsp` 3D FFTs
- [ ] Leapfrog integrator
- [ ] Langevin thermostat
- [ ] Nosé-Hoover thermostat
- [ ] Berendsen thermostat/barostat
- [ ] Unit tests + doctests for new integrators/thermostats/PME
- [ ] Rustdoc updates
- [ ] `cargo fmt` / `clippy` clean
- [ ] `cargo deny check` clean

### tpt-chem-verify (completion — Kani proofs)

*Depends on: `kani` (dev-only verification tool).*

- [ ] Wire deps: `kani` (Kani Rust Verifier, dev-dependency / separate verification harness crate)
- [ ] Kani proof: neighbor-list cell indexing is always in-bounds
- [ ] Kani proof: SCF iteration density matrices remain positive semi-definite
- [ ] Kani proof: absence of panics in ERI evaluation, including r→0 singularity handling
- [ ] proptest: symplectic integrators (Velocity Verlet, Leapfrog) conserve phase-space volume
- [ ] Unit tests + doctests
- [ ] Rustdoc
- [ ] `cargo fmt` / `clippy` clean
- [ ] `cargo deny check` clean

### Phase 3 exit check

- [ ] `cargo build` + `cargo test` green across the whole workspace
- [ ] All Kani harnesses pass (bounded model checking green)
- [ ] PME long-range electrostatics validated against Ewald direct-sum reference on a small test system

---

## Phase 4 — Umbrella Crate & Workspace Polish

*Not explicitly phased in spec.txt §6 — added as a sensible closing phase,
matching how sibling repos (e.g. `tpt-math`) finish with an umbrella
re-export crate and a full workspace hygiene pass.*

### tpt-chemistry (umbrella crate)

*Feature-gated umbrella re-exporting all `tpt-chem-*` crates.*

- [ ] Scaffold `crates/tpt-chemistry/`
- [ ] Wire optional deps + matching Cargo features per constituent crate (`core`, `md`, `quantum`, `kinetics`, `crystal`, `io`, `verify`)
- [ ] Re-export each constituent's public API behind its feature
- [ ] Rustdoc documenting the feature matrix
- [ ] `cargo fmt` / `clippy` / `deny` clean across feature combinations (test each feature subset, not just `--all-features`)

### Workspace-wide polish

- [ ] Workspace-wide `cargo deny check` pass confirming zero Apache-2.0-only dependencies anywhere in the graph
- [ ] README pass covering all 8 crates
- [ ] README: integration/synergy notes per spec.txt §5 — `tpt-materials` (micro-to-macro property feed), `tpt-math` (linear algebra/ODE), `tpt-dsp` (3D FFTs for PME + vibrational DOS), `tpt-fem`/`tpt-physics` (QM/MM coupling), AI-native agents (`tpt-eve`/`tpt-anima` — autodiff hooks for MLIPs / differentiable molecular design)
- [ ] Decide on and document crates.io / GitHub publishing status (or explicitly mark out-of-scope for this pass, as `tpt-math` initially did)
