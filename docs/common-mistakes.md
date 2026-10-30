# Common mistakes

The workspace uses different unit systems in different crates on purpose.
Most wrong answers come from mixing them. This page lists the conventions and
the traps.

Written for: developers and researchers calling tpt-chemistry from their own code.

## Units by crate

| Crate | Lengths | Energies | Time | Notes |
| --- | --- | --- | --- | --- |
| `tpt-chem-core` / `tpt-chem-md` | Å | kJ·mol⁻¹ | fs | MD unit system; charges in e |
| `tpt-chem-quantum` | Å in input, Bohr internally | Hartree | atomic units in Molden | Gradients in Hartree/Bohr |
| `tpt-chem-io` (`.gro`) | nm in file | nm at the file boundary | — | Converted to Å on read; written back as nm |
| `tpt-chem-aqueous` | — | `log₁₀ K` at 25 °C | — | Concentrations are molalities (mol per kg water) |

Conversion helpers live in `tpt_chem_core::units`, for example
`angstrom_to_bohr`, `bohr_to_angstrom`, `hartree_to_kj_per_mol` and
`atomic_time_to_fs`.

## Traps

**Å vs Bohr.** The HF examples take geometries in Å but report minima in Bohr
(`h2_hf` prints both). Passing a Bohr value where the API expects Å moves
every atom by a factor of about 1.89. If an HF geometry looks absurdly
stretched, check this first.

**Hartree vs kJ·mol⁻¹.** One Hartree is about 2625.5 kJ·mol⁻¹. Quantum
energies never come out in kJ·mol⁻¹. Convert before comparing with MD or
thermochemistry values.

**Molality vs molarity.** `tpt-chem-aqueous` concentrations are mol per kg
of water. Molarity (mol per litre of solution) differs by the solution
density, which matters at high ionic strength. Mixing them shifts saturation
indices.

**Activity vs concentration.** `pH` in `tpt-chem-aqueous` is on the activity
scale, `a(H⁺) = 10^(−pH)`. Setting pH to a measured value that was taken on a
concentration scale gives a different equilibrium.

**Saturation index is not precipitation.** `SI > 0` only says a phase can
lower its free energy. Visible precipitation also needs nucleation.
Use `CalciumPhosphateSystem::kinetic_onset_ph` for a kinetic estimate, and
treat its default `NucleationKinetics` as illustrative until calibrated.

**Illustrative constants.** Phosphate, chelator and nucleation defaults in
`tpt-chem-aqueous` are placeholders. They are overridable and should be
replaced with measured values before any quantitative use.

**Basis set scope.** STO-3G and 6-31G are implemented for H–Ne. d-function
basis sets (cc-pVDZ and similar) and heavier elements are not yet supported.
