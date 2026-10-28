# Sample structures

Small geometries for trying the readers and examples. All coordinates are in Å.

| File | Contents | Source |
| --- | --- | --- |
| `water.xyz` / `.pdb` | H₂O, gas-phase geometry | Idealized: O–H 0.96 Å, angle ≈ 104.5° |
| `benzene.xyz` / `.pdb` | C₆H₆, D₆h | Idealized hexagon: C–C 1.39 Å, C–H 1.09 Å; PDB carries CONECT bonds |
| `nacl_diatomic.xyz` / `.pdb` | NaCl gas-phase diatomic | Bond length 2.36 Å |
| `argon_fcc.xyz` / `.cif` | Solid argon, 256 atoms (4×4×4 FCC conventional cells) | Ideal lattice, a = 5.26 Å; CIF is P1 with cell edge 21.04 Å |
| `silicon_diamond.xyz` / `.cif` | Diamond-cubic Si, 64 atoms (2×2×2 conventional cells) | Ideal lattice, a = 5.431 Å; CIF is P1 with cell edge 10.862 Å |

The CIF files list every atom in the cell (`P 1`), so they can go straight into
`tpt_chemistry::pipeline::xrd_from_cif` without symmetry expansion.

These are idealized starting structures, not experimental or optimized
geometries. Run a geometry optimization before using them for reference
energies.
