# tpt-chemistry

tpt-chemistry is a pure-Rust library for computational chemistry: classical
molecular dynamics, Hartree–Fock, reaction kinetics, crystallography, aqueous
thermodynamics, and file formats. It has no C/C++ FFI and no external solver.
Every algorithm (integrals, SCF, Ewald and PME, neighbour lists, ODE and
stochastic solvers) is implemented in this workspace.

## What works today

- **MD:** velocity Verlet, leapfrog, Langevin, Nosé–Hoover, Berendsen, Ewald and PME,
  SHAKE/RATTLE, minimisation (steepest descent and L-BFGS), observables, and a
  radial machine-learned potential that plugs into the same integrator.
- **Quantum:** RHF and UHF with STO-3G and 6-31G (H–Ne), analytic gradients,
  geometry optimisation, cube and Molden output.
- **Kinetics:** mass-action ODEs, Gillespie, tau-leaping, next-reaction SSA,
  Arrhenius and Eyring rates, parameter fitting.
- **Crystal:** Bravais lattices, general-position operations for all 230 space
  groups, Wigner–Seitz cells, Brillouin-zone grids, XRD simulation.
- **Aqueous:** Debye–Hückel, Davies, Pitzer, and an eNRTL skeleton, with
  calcium-phosphate speciation. Its constants are illustrative, not fitted.
- **I/O:** XYZ, PDB, MOL2, SDF, CIF, GRO, SMILES (read and canonical write),
  Gaussian input and output, and streaming XYZ trajectories.

Open items are listed in `todo.md` in the repository. Two that matter for
users: XTC trajectories and OPLS/UFF parameters. Kani proofs run on Linux CI
only.

## Licence

Every crate is `MIT OR Apache-2.0`. The dependency policy is stricter: no
Apache-2.0-only crate may enter the graph (see `deny.toml`).
