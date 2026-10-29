# Units cheat sheet

Each crate uses the units listed below. Convert at the boundary, not inside
a calculation.

| Quantity | MD crates | Quantum | Conversion |
| --- | --- | --- | --- |
| Length | Å | Å in, Bohr inside | 1 Bohr = 0.529177210544 Å |
| Energy | kJ·mol⁻¹ | Hartree | 1 Eₕ = 2625.499639481 kJ·mol⁻¹ |
| Energy | kcal·mol⁻¹ | Hartree | 1 Eₕ = 627.5094740631 kcal·mol⁻¹ |
| Time | fs | atomic units | 1 a.u. of time = 0.024188843265864 fs |
| Pressure | bar (from kJ·mol⁻¹·Å⁻³) | — | `KJ_MOL_ANG3_TO_BAR` |

The constants and helpers are in `tpt_chem_core::units`, for example
`bohr_to_angstrom`, `hartree_to_kj_per_mol`, and `atomic_time_to_fs`.
