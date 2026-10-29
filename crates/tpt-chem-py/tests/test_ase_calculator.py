"""ASE calculator: energies and forces agree with the native solver and with finite differences.

ASE is LGPL-2.1-or-later and is not a dependency of this project. These tests
run only when ASE is installed and are skipped otherwise.
"""

import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import tpt_chem  # noqa: E402

try:
    from ase import Atoms
    from ase.units import Hartree

    from ase_calculator import TptRhf

    HAVE_ASE = True
except ImportError:
    HAVE_ASE = False

WATER = "3\nwater\nO 0 0 0.117\nH 0 0.757 -0.469\nH 0 -0.757 -0.469\n"


def water_atoms():
    atoms = Atoms(
        "OH2",
        positions=[(0.0, 0.0, 0.117), (0.0, 0.757, -0.469), (0.0, -0.757, -0.469)],
    )
    atoms.calc = TptRhf()
    return atoms


@unittest.skipUnless(HAVE_ASE, "ASE is not installed (optional, LGPL-2.1)")
class AseCalculator(unittest.TestCase):
    def test_energy_matches_native_solver(self):
        atoms = water_atoms()
        self.assertAlmostEqual(
            atoms.get_potential_energy(), tpt_chem.rhf_energy_xyz(WATER) * Hartree, places=6
        )

    def test_forces_match_finite_differences(self):
        h = 1e-4  # Å
        atoms = water_atoms()
        forces = atoms.get_forces()
        numeric = np.zeros_like(forces)
        for i in range(len(atoms)):
            for axis in range(3):
                plus = atoms.copy()
                minus = atoms.copy()
                plus.calc = TptRhf()
                minus.calc = TptRhf()
                plus.positions[i, axis] += h
                minus.positions[i, axis] -= h
                numeric[i, axis] = -(plus.get_potential_energy() - minus.get_potential_energy()) / (2 * h)
        np.testing.assert_allclose(forces, numeric, atol=2e-4)

    def test_net_force_vanishes(self):
        # Translation invariance: the forces on an isolated molecule sum to zero.
        forces = water_atoms().get_forces()
        np.testing.assert_allclose(forces.sum(axis=0), 0.0, atol=1e-6)


if __name__ == "__main__":
    unittest.main()
