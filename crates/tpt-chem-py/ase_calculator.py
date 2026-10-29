"""ASE calculator backed by tpt-chemistry's RHF/STO-3G energies and forces.

    from ase.build import molecule
    from ase_calculator import TptRhf

    atoms = molecule("H2O")
    atoms.calc = TptRhf()
    atoms.get_potential_energy()   # eV

Positions are Å, energies are converted from Hartree to eV, and forces are
converted from Hartree/Å to eV/Å. Only elements supported by STO-3G (H–Ne) and
neutral closed-shell molecules are accepted.

This file is not part of the wheel. Put it on `sys.path` next to `tpt_chem`.
"""

import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from ase.units import Hartree

import tpt_chem


def _to_xyz(atoms) -> str:
    lines = [str(len(atoms)), "ase"]
    for symbol, (x, y, z) in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        lines.append(f"{symbol} {x:.12f} {y:.12f} {z:.12f}")
    return "\n".join(lines) + "\n"


class TptRhf(Calculator):
    """RHF/STO-3G energies and analytic forces."""

    implemented_properties = ["energy", "forces"]

    def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        energy, forces = tpt_chem.rhf_energy_forces_xyz(_to_xyz(self.atoms))
        self.results["energy"] = energy * Hartree
        self.results["forces"] = np.array(forces) * Hartree
