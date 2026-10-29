# Third-party software and data

This page lists what the project depends on, what it uses only as a
reference or generator, and what has been removed. The policy is that the
project uses no GPL-licensed code and no Apache-2.0-only code. Dependencies
must be MIT, BSD-style, ISC, Zlib, CC0, Unicode, MPL-2.0, or a dual licence
that includes one of these.

## Shipped Rust dependencies

The dependency graph is checked by `cargo deny check` against the allow-list
in `deny.toml`. The check fails if any crate's licence falls outside the list,
so a GPL or Apache-2.0-only crate cannot enter unnoticed. The allow-list is
MIT, BSD-2-Clause, BSD-3-Clause, ISC, Zlib, CC0-1.0, Unicode-3.0,
Unicode-DFS-2016, MPL-2.0, and `Apache-2.0 WITH LLVM-exception`. Dual
`MIT OR Apache-2.0` crates satisfy the list through their MIT alternative.

The main crates are `serde`, `serde_json`, `rayon`, `libm`, `proptest` and
`tpt-math` (path dependency, same licence as this workspace). `criterion`
(dev-dependency), `pyo3` and `wasm-bindgen` are used only by the benchmark,
Python and WebAssembly crates.

## Development tools and reference data

These are not shipped. Each one generates test fixtures, generates data
files, or provides an independent reference for a test.

| Tool or source | Licence | Used for |
| --- | --- | --- |
| OpenMM 8.x | MIT | Reference energies and forces (AMBER, OPLS, Ewald/PME, LJ), the DCD fixture (`tools/openmm_dcd_fixture.py`), the AMBER reference (`tools/amber_reference.py`) |
| amber14 force-field XML (distributed with OpenMM) | Parameters published by the AMBER project; the XML ships under OpenMM's MIT licence | Atom types, charges and bonded terms for the AMBER subset |
| foyer `oplsaa.xml` (mosdef-hub/foyer, 0.7.6) | MIT | OPLS-AA parameters for the butane reference (`tools/opls_reference.py`) |
| spglib 2.x | BSD-3-Clause | Symmetry operations for the 230 space groups (`tools/gen_spacegroups.py`) |
| basis_set_exchange 0.12 | BSD-3-Clause | cc-pVDZ basis data for H–Ar (`tools/gen_cc_pvdz.py`) |
| numpy | BSD-3-Clause | Reference scripts |
| pytest | MIT | Python tests |
| maturin | MIT or Apache-2.0 | Building the Python wheel |
| CRC Handbook of Chemistry and Physics, 97th ed. | Published data | Measured mean activity coefficients for NaCl, KCl, CaCl₂ and MgCl₂ (`activity/pitzer.rs`). The numbers are entered as facts |
| PHREEQC `pitzer.dat` | Published database; check PHREEQC's terms before redistributing | The 25 °C MgCl₂ Pitzer values, cited in `activity/pitzer.rs` |

## Optional, not a dependency

| Tool | Licence | Status |
| --- | --- | --- |
| ASE (Atomic Simulation Environment) | LGPL-2.1-or-later | `crates/tpt-chem-py/ase_calculator.py` is a calculator for ASE users. The Python test for it skips when ASE is not installed. ASE is not required by the project or its default tests |

## Removed

| Tool | Licence | What changed |
| --- | --- | --- |
| MDAnalysis | LGPL-3.0-or-later | The DCD cross-check scripts and fixture were replaced with OpenMM's `DCDFile` output (`tests/dcd_openmm.rs`). The DCD fixture from MDAnalysis is gone |
| pyEQL | LGPL-3.0-or-later | Its benchmark CSV had been read to find the CRC values. No file from it is in the repository now. The CRC numbers are cited to the Handbook |

## Provenance note on the DCD module

The first version of `tpt-chem-io/src/dcd.rs` was written after reading the
header code in MDAnalysis's bundled `readdcd.h`, to confirm the record layout.
That code was not copied. The module has since been checked against OpenMM's
independently written writer, and the frame records match byte for byte. The
reader and writer in `dcd.rs` are the project's own code. A future review that
wants to be strict about provenance can compare `dcd.rs` against the CHARMM and
VMD DCD format descriptions and the OpenMM output in `tests/data/`.
